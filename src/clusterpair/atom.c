/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <allocate.h>
#include <atom.h>
#include <force.h>
#include <parameter.h>
#include <util.h>

inline int get_ncj_from_nci(int nci) {
#ifdef CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
    return nci << 3;
#else
#if CLUSTER_M == CLUSTER_N
    return nci;
#elif CLUSTER_M < CLUSTER_N
    return nci >> 1;
#else
    return nci << 1;
#endif
#endif
}

int write_atoms_to_file(Atom* atom, char* name) {
    FILE *fp = fopen(name, "wb");
    if (fp == NULL) {
        perror("Error opening file");
        return -1;
    }

    for (int i = 0; i < atom->Nlocal; ++i) {
        fprintf(fp,
            "%lf %lf %lf %lf %lf %lf %d\n",
            atom_x(i),
            atom_y(i),
            atom_z(i),
            atom_vx(i),
            atom_vy(i),
            atom_vz(i),
            atom->type[i]);
    }
    fclose(fp);
    return 0;
}

void initAtom(Atom* atom) {
    atom->x               = NULL;
    atom->y               = NULL;
    atom->z               = NULL;
    atom->vx              = NULL;
    atom->vy              = NULL;
    atom->vz              = NULL;
    atom->cl_x            = NULL;
    atom->cl_v            = NULL;
    atom->cl_f            = NULL;
    atom->cl_t            = NULL;
    atom->cl_sqrt_epsilon = NULL;
    atom->cl_sigma3       = NULL;
    atom->Natoms          = 0;
    atom->Nlocal          = 0;
    atom->Nghost          = 0;
    atom->Nmax            = 0;
    atom->Nclusters_local = 0;
    atom->Nclusters_ghost = 0;
    atom->Nclusters_max   = 0;
    atom->type            = NULL;
    atom->ntypes          = 0;
    atom->epsilon         = NULL;
    atom->sigma6          = NULL;
    atom->cutforcesq      = NULL;
    atom->cutneighsq      = NULL;
    atom->sqrt_epsilon_per_type = NULL;
    atom->sigma3_per_type = NULL;
    atom->iclusters       = NULL;
    atom->jclusters       = NULL;
    atom->cluster_bin     = NULL;
    atom->siclusters      = NULL;

    initMasks(atom);
    // MPI New features
    Box* mybox      = &(atom->mybox);
    atom->NmaxGhost = 0;
    atom->PBCx      = NULL;
    atom->PBCy      = NULL;
    atom->PBCz      = NULL;
    mybox->xprd     = 0;
    mybox->yprd     = 0;
    mybox->zprd     = 0;
    mybox->lo[0]    = 0;
    mybox->lo[1]    = 0;
    mybox->lo[2]    = 0;
    mybox->hi[0]    = 0;
    mybox->hi[1]    = 0;
    mybox->hi[2]    = 0;
}

void createAtom(Atom* atom, Parameter* param) {
    int me = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    MD_FLOAT xlo  = 0.0;
    MD_FLOAT xhi  = param->xprd;
    MD_FLOAT ylo  = 0.0;
    MD_FLOAT yhi  = param->yprd;
    MD_FLOAT zlo  = 0.0;
    MD_FLOAT zhi  = param->zprd;
    atom->Natoms  = 4 * param->nx * param->ny * param->nz;
    atom->Nlocal  = 0;
    atom->ntypes  = param->ntypes;
    atom->epsilon = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma6  = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutforcesq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutneighsq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));

    // Compute type-pair LJ parameters
    atom->sqrt_epsilon_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma3_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    computePerTypeLJParameters(atom->ntypes, param, atom->sqrt_epsilon_per_type, atom->sigma3_per_type);
    computeTypePairLJParameters(atom->ntypes, atom->sqrt_epsilon_per_type, atom->sigma3_per_type, atom->epsilon, atom->sigma6);
    // Keep per-type arrays for populating cluster LJ parameters
    for (int i = 0; i < atom->ntypes * atom->ntypes; i++) {
        atom->cutneighsq[i] = param->cutneigh * param->cutneigh;
        atom->cutforcesq[i] = param->cutforce * param->cutforce;
    }

    MD_FLOAT alat = pow((4.0 / param->rho), (1.0 / 3.0));
    int ilo       = (int)(xlo / (0.5 * alat) - 1);
    int ihi       = (int)(xhi / (0.5 * alat) + 1);
    int jlo       = (int)(ylo / (0.5 * alat) - 1);
    int jhi       = (int)(yhi / (0.5 * alat) + 1);
    int klo       = (int)(zlo / (0.5 * alat) - 1);
    int khi       = (int)(zhi / (0.5 * alat) + 1);

    ilo = MAX(ilo, 0);
    ihi = MIN(ihi, 2 * param->nx - 1);
    jlo = MAX(jlo, 0);
    jhi = MIN(jhi, 2 * param->ny - 1);
    klo = MAX(klo, 0);
    khi = MIN(khi, 2 * param->nz - 1);

    MD_FLOAT xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp;
    int i, j, k, m, n;
    int sx        = 0;
    int sy        = 0;
    int sz        = 0;
    int ox        = 0;
    int oy        = 0;
    int oz        = 0;
    int subboxdim = 8;

    if(param->setup) {
        while (oz * subboxdim <= khi) {
            k = oz * subboxdim + sz;
            j = oy * subboxdim + sy;
            i = ox * subboxdim + sx;

            if (((i + j + k) % 2 == 0) && (i >= ilo) && (i <= ihi) && (j >= jlo) &&
                (j <= jhi) && (k >= klo) && (k <= khi)) {
                xtmp = 0.5 * alat * i;
                ytmp = 0.5 * alat * j;
                ztmp = 0.5 * alat * k;

                if (xtmp >= xlo && xtmp < xhi && ytmp >= ylo && ytmp < yhi &&
                    ztmp >= zlo && ztmp < zhi) {
                    n = k * (2 * param->ny) * (2 * param->nx) + j * (2 * param->nx) + i +
                        1;
                    for (m = 0; m < 5; m++) {
                        myrandom(&n);
                    }
                    vxtmp = myrandom(&n);
                    for (m = 0; m < 5; m++) {
                        myrandom(&n);
                    }
                    vytmp = myrandom(&n);
                    for (m = 0; m < 5; m++) {
                        myrandom(&n);
                    }
                    vztmp = myrandom(&n);

                    if (atom->Nlocal == atom->Nmax) {
                        growAtom(atom);
                    }

                    atom_x(atom->Nlocal)     = xtmp;                 
                    atom_y(atom->Nlocal)     = ytmp;
                    atom_z(atom->Nlocal)     = ztmp;
                    atom->vx[atom->Nlocal]   = vxtmp;
                    atom->vy[atom->Nlocal]   = vytmp;
                    atom->vz[atom->Nlocal]   = vztmp;
                    atom->type[atom->Nlocal] = rand() % atom->ntypes;
                    atom->Nlocal++;
                }
            }

            sx++;
            if (sx == subboxdim) {
                sx = 0;
                sy++;
            }
            if (sy == subboxdim) {
                sy = 0;
                sz++;
            }
            if (sz == subboxdim) {
                sz = 0;
                ox++;
            }
            if (ox * subboxdim > ihi) {
                ox = 0;
                oy++;
            }
            if (oy * subboxdim > jhi) {
                oy = 0;
                oz++;
            }
        }
    }
}

int typeStr2int(const char* type) {
    if (strncmp(type, "Ar", 2) == 0) {
        return 0;
    } // Argon
    fprintf(stderr, "Invalid atom type: %s\n", type);
    exit(-1);
    return -1;
}

int readAtom(Atom* atom, Parameter* param) {
    int me = 0;

#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    int len = strlen(param->input_file);

    if (strncmp(&param->input_file[len - 4], ".pdb", 4) == 0) {
        return readAtomPdb(atom, param);
    }

    if (strncmp(&param->input_file[len - 4], ".gro", 4) == 0) {
        return readAtomGro(atom, param);
    }

    if (strncmp(&param->input_file[len - 4], ".dmp", 4) == 0) {
        return readAtomDmp(atom, param);
    }

    if (me == 0) {
        fprintf(stderr,
            "Invalid input file extension: %s\nValid choices are: pdb, gro, dmp\n",
            param->input_file);
    }

    exit(-1);
    return -1;
}

int readAtomPdb(Atom* atom, Parameter* param) {
    int me = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me); // New
#endif

    FILE* fp = fopen(param->input_file, "r");
    char line[MAXLINE];
    int readAtoms = 0;

    if (!fp) {
        if (me == 0) {
            fprintf(stderr, "Could not open input file: %s\n", param->input_file);
        }

        exit(-1);
        return -1;
    }

    while (fgets(line, MAXLINE, fp) != NULL) {
        char* item = strtok(line, " ");

        if (strncmp(item, "CRYST1", 6) == 0) {
            param->xlo  = 0.0;
            param->xhi  = atof(strtok(NULL, " "));
            param->ylo  = 0.0;
            param->yhi  = atof(strtok(NULL, " "));
            param->zlo  = 0.0;
            param->zhi  = atof(strtok(NULL, " "));
            param->xprd = param->xhi - param->xlo;
            param->yprd = param->yhi - param->ylo;
            param->zprd = param->zhi - param->zlo;
            // alpha, beta, gamma, sGroup, z
        } else if (strncmp(item, "ATOM", 4) == 0) {
            char* label;
            int atomId, compId;
            MD_FLOAT occupancy, charge;
            atomId = atoi(strtok(NULL, " ")) - 1;

            while (atomId + 1 >= atom->Nmax) {
                growAtom(atom);
            }

            atom->type[atomId] = typeStr2int(strtok(NULL, " "));
            label              = strtok(NULL, " ");
            compId             = atoi(strtok(NULL, " "));
            atom_x(atomId)     = atof(strtok(NULL, " "));
            atom_y(atomId)     = atof(strtok(NULL, " "));
            atom_z(atomId)     = atof(strtok(NULL, " "));
            atom->vx[atomId]   = 0.0;
            atom->vy[atomId]   = 0.0;
            atom->vz[atomId]   = 0.0;
            occupancy          = atof(strtok(NULL, " "));
            charge             = atof(strtok(NULL, " "));
            atom->ntypes       = MAX(atom->type[atomId] + 1, atom->ntypes);
            atom->Natoms++;
            atom->Nlocal++;
            readAtoms++;
        } else if (strncmp(item, "HEADER", 6) == 0 || strncmp(item, "REMARK", 6) == 0 ||
                   strncmp(item, "MODEL", 5) == 0 || strncmp(item, "TER", 3) == 0 ||
                   strncmp(item, "ENDMDL", 6) == 0) {
            // Do nothing
        } else {
            if (me == 0) {
                fprintf(stderr, "Invalid item: %s\n", item);
            }

            exit(-1);
            return -1;
        }
    }

    if (!readAtoms) {
        if (me == 0) {
            fprintf(stderr, "Input error: No atoms read!\n");
        }
        exit(-1);
        return -1;
    }

    atom->epsilon = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma6  = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutforcesq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutneighsq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    // Compute type-pair LJ parameters
    atom->sqrt_epsilon_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma3_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    computePerTypeLJParameters(atom->ntypes, param, atom->sqrt_epsilon_per_type, atom->sigma3_per_type);
    computeTypePairLJParameters(atom->ntypes, atom->sqrt_epsilon_per_type, atom->sigma3_per_type, atom->epsilon, atom->sigma6);
    // Keep per-type arrays for populating cluster LJ parameters
    for (int i = 0; i < atom->ntypes * atom->ntypes; i++) {
        atom->cutneighsq[i] = param->cutneigh * param->cutneigh;
        atom->cutforcesq[i] = param->cutforce * param->cutforce;
    }

    if (me == 0) {
        fprintf(stdout, "Read %d atoms from %s\n", readAtoms, param->input_file);
    }
    fclose(fp);
    return readAtoms;
}

int readAtomGro(Atom* atom, Parameter* param) {
    int me = 0;

#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me); // New
#endif

    FILE* fp = fopen(param->input_file, "r");
    char line[MAXLINE];
    char desc[MAXLINE];
    int readAtoms   = 0;
    int atomsToRead = 0;
    int i           = 0;

    if (!fp) {
        if (me == 0) {
            fprintf(stderr, "Could not open input file: %s\n", param->input_file);
        }
        exit(-1);
        return -1;
    }

    if (fgets(desc, MAXLINE, fp) == NULL && ferror(fp) != 0) {
        fprintf(stderr, "Error in fgets function\n");
    }

    for (i = 0; desc[i] != '\n'; i++)
        ;
    desc[i] = '\0';

    if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
        fprintf(stderr, "Error in fgets function\n");
    }

    atomsToRead = atoi(strtok(line, " "));
    if (me == 0) {
        fprintf(stdout, "System: %s with %d atoms\n", desc, atomsToRead);
    }

    while (readAtoms < atomsToRead) {
        if (fgets(line, MAXLINE, fp) == NULL) {
            break;
        }

        char* label = strtok(line, " ");
        int type    = typeStr2int(strtok(NULL, " "));
        int atomId  = atoi(strtok(NULL, " ")) - 1;
        atomId      = readAtoms;
        while (atomId + 1 >= atom->Nmax) {
            growAtom(atom);
        }

        atom->type[atomId] = type;
        atom_x(atomId)     = atof(strtok(NULL, " "));
        atom_y(atomId)     = atof(strtok(NULL, " "));
        atom_z(atomId)     = atof(strtok(NULL, " "));
        atom->vx[atomId]   = atof(strtok(NULL, " "));
        atom->vy[atomId]   = atof(strtok(NULL, " "));
        atom->vz[atomId]   = atof(strtok(NULL, " "));
        atom->ntypes       = MAX(atom->type[atomId] + 1, atom->ntypes);
        atom->Natoms++;
        atom->Nlocal++;
        readAtoms++;
    }

    if (fgets(line, MAXLINE, fp) != NULL) {
        param->xlo  = 0.0;
        param->xhi  = atof(strtok(line, " "));
        param->ylo  = 0.0;
        param->yhi  = atof(strtok(NULL, " "));
        param->zlo  = 0.0;
        param->zhi  = atof(strtok(NULL, " "));
        param->xprd = param->xhi - param->xlo;
        param->yprd = param->yhi - param->ylo;
        param->zprd = param->zhi - param->zlo;
    }

    if (readAtoms != atomsToRead) {
        if (me == 0) {
            fprintf(stderr,
                "Input error: Number of atoms read do not match (%d/%d).\n",
                readAtoms,
                atomsToRead);
        }
        exit(-1);
        return -1;
    }

    atom->epsilon = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma6  = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutforcesq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutneighsq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    // Compute type-pair LJ parameters
    atom->sqrt_epsilon_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma3_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    computePerTypeLJParameters(atom->ntypes, param, atom->sqrt_epsilon_per_type, atom->sigma3_per_type);
    computeTypePairLJParameters(atom->ntypes, atom->sqrt_epsilon_per_type, atom->sigma3_per_type, atom->epsilon, atom->sigma6);
    // Keep per-type arrays for populating cluster LJ parameters
    for (int i = 0; i < atom->ntypes * atom->ntypes; i++) {
        atom->cutneighsq[i] = param->cutneigh * param->cutneigh;
        atom->cutforcesq[i] = param->cutforce * param->cutforce;
    }

    if (me == 0) {
        fprintf(stdout, "Read %d atoms from %s\n", readAtoms, param->input_file);
    }

    fclose(fp);
    return readAtoms;
}

int readAtomDmp(Atom* atom, Parameter* param) {
    int me = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    FILE* fp = fopen(param->input_file, "r");
    char line[MAXLINE];
    int natoms    = 0;
    int readAtoms = 0;
    int atomId    = -1;
    int ts        = -1;

    if (!fp) {
        if (me == 0) {
            fprintf(stderr, "Could not open input file: %s\n", param->input_file);
        }
        exit(-1);
        return -1;
    }

    while (fgets(line, MAXLINE, fp) != NULL && ts < 1 && !readAtoms) {
        if (strncmp(line, "ITEM: ", 6) == 0) {
            char* item = &line[6];

            if (strncmp(item, "TIMESTEP", 8) == 0) {
                if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
                    fprintf(stderr, "Error in fgets function\n");
                }

                ts = atoi(line);
            } else if (strncmp(item, "NUMBER OF ATOMS", 15) == 0) {
                if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
                    fprintf(stderr, "Error in fgets function\n");
                }
                natoms       = atoi(line);
                atom->Natoms = natoms;
                atom->Nlocal = natoms;
                while (atom->Nlocal >= atom->Nmax) {
                    growAtom(atom);
                }
            } else if (strncmp(item, "BOX BOUNDS pp pp pp", 19) == 0) {
                if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
                    fprintf(stderr, "Error in fgets function\n");
                }
                param->xlo  = atof(strtok(line, " "));
                param->xhi  = atof(strtok(NULL, " "));
                param->xprd = param->xhi - param->xlo;

                if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
                    fprintf(stderr, "Error in fgets function\n");
                }
                param->ylo  = atof(strtok(line, " "));
                param->yhi  = atof(strtok(NULL, " "));
                param->yprd = param->yhi - param->ylo;

                if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
                    fprintf(stderr, "Error in fgets function\n");
                }
                param->zlo  = atof(strtok(line, " "));
                param->zhi  = atof(strtok(NULL, " "));
                param->zprd = param->zhi - param->zlo;
            } else if (strncmp(item, "ATOMS id type x y z vx vy vz", 28) == 0) {
                for (int i = 0; i < natoms; i++) {
                    if (fgets(line, MAXLINE, fp) == NULL && ferror(fp) != 0) {
                        fprintf(stderr, "Error in fgets function\n");
                    }
                    atomId             = atoi(strtok(line, " ")) - 1;
                    atom->type[atomId] = atoi(strtok(NULL, " "));
                    atom_x(atomId)     = atof(strtok(NULL, " "));
                    atom_y(atomId)     = atof(strtok(NULL, " "));
                    atom_z(atomId)     = atof(strtok(NULL, " "));
                    atom->vx[atomId]   = atof(strtok(NULL, " "));
                    atom->vy[atomId]   = atof(strtok(NULL, " "));
                    atom->vz[atomId]   = atof(strtok(NULL, " "));
                    atom->ntypes       = MAX(atom->type[atomId] + 1, atom->ntypes);
                    readAtoms++;
                }
            } else {
                if (me == 0) {
                    fprintf(stderr, "Invalid item: %s\n", item);
                }
                exit(-1);
                return -1;
            }
        } else {
            if (me == 0) {
                fprintf(stderr,
                    "Invalid input from file, expected item reference but got:\n%s\n",
                    line);
            }
            exit(-1);
            return -1;
        }
    }

    if (ts < 0 || !natoms || !readAtoms) {
        if (me == 0) {
            fprintf(stderr, "Input error: atom data was not read!\n");
        }
        exit(-1);
        return -1;
    }

    atom->epsilon = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma6  = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutforcesq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutneighsq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    // Compute type-pair LJ parameters
    atom->sqrt_epsilon_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma3_per_type = allocate(ALIGNMENT, atom->ntypes * sizeof(MD_FLOAT));
    computePerTypeLJParameters(atom->ntypes, param, atom->sqrt_epsilon_per_type, atom->sigma3_per_type);
    computeTypePairLJParameters(atom->ntypes, atom->sqrt_epsilon_per_type, atom->sigma3_per_type, atom->epsilon, atom->sigma6);
    // Keep per-type arrays for populating cluster LJ parameters
    for (int i = 0; i < atom->ntypes * atom->ntypes; i++) {
        atom->cutneighsq[i] = param->cutneigh * param->cutneigh;
        atom->cutforcesq[i] = param->cutforce * param->cutforce;
    }
    if (me == 0) {
        fprintf(stdout, "Read %d atoms from %s\n", natoms, param->input_file);
    }
    fclose(fp);
    return natoms;
}

void initMasks(Atom* atom) {
    const unsigned int halfMaskBits = VECTOR_WIDTH >> 1;
    unsigned int mask0, mask1, mask2, mask3;

    atom->exclusion_filter        = allocate(ALIGNMENT,
        CLUSTER_M * VECTOR_WIDTH * sizeof(MD_UINT));
    atom->diagonal_4xn_j_minus_i  = allocate(ALIGNMENT,
        MAX(CLUSTER_M, VECTOR_WIDTH) * sizeof(MD_UINT));
    atom->diagonal_2xnn_j_minus_i = allocate(ALIGNMENT, VECTOR_WIDTH * sizeof(MD_UINT));
    // atom->masks_2xnn = allocate(ALIGNMENT, 8 * sizeof(unsigned int));

    for (int j = 0; j < MAX(CLUSTER_M, VECTOR_WIDTH); j++) {
        atom->diagonal_4xn_j_minus_i[j] = (MD_FLOAT)(j)-0.5;
    }

    for (int j = 0; j < VECTOR_WIDTH / 2; j++) {
        atom->diagonal_2xnn_j_minus_i[j]                    = (MD_FLOAT)(j)-0.5;
        atom->diagonal_2xnn_j_minus_i[VECTOR_WIDTH / 2 + j] = (MD_FLOAT)(j - 1) - 0.5;
    }

    for (int i = 0; i < CLUSTER_M * VECTOR_WIDTH; i++) {
        atom->exclusion_filter[i] = (1U << i);
    }

#if CLUSTER_M == CLUSTER_N
    for (unsigned int cond0 = 0; cond0 < 2; cond0++) {
        mask0                              = (unsigned int)(0xf - 0x1 * cond0);
        mask1                              = (unsigned int)(0xf - 0x3 * cond0);
        mask2                              = (unsigned int)(0xf - 0x7 * cond0);
        mask3                              = (unsigned int)(0xf - 0xf * cond0);
        atom->masks_2xnn_hn[cond0 * 2 + 0] = (mask1 << halfMaskBits) | mask0;
        atom->masks_2xnn_hn[cond0 * 2 + 1] = (mask3 << halfMaskBits) | mask2;

        mask0                              = (unsigned int)(0xf - 0x1 * cond0);
        mask1                              = (unsigned int)(0xf - 0x2 * cond0);
        mask2                              = (unsigned int)(0xf - 0x4 * cond0);
        mask3                              = (unsigned int)(0xf - 0x8 * cond0);
        atom->masks_2xnn_fn[cond0 * 2 + 0] = (mask1 << halfMaskBits) | mask0;
        atom->masks_2xnn_fn[cond0 * 2 + 1] = (mask3 << halfMaskBits) | mask2;

        atom->masks_4xn_hn[cond0 * 4 + 0] = (unsigned int)(0xf - 0x1 * cond0);
        atom->masks_4xn_hn[cond0 * 4 + 1] = (unsigned int)(0xf - 0x3 * cond0);
        atom->masks_4xn_hn[cond0 * 4 + 2] = (unsigned int)(0xf - 0x7 * cond0);
        atom->masks_4xn_hn[cond0 * 4 + 3] = (unsigned int)(0xf - 0xf * cond0);

        atom->masks_4xn_fn[cond0 * 4 + 0] = (unsigned int)(0xf - 0x1 * cond0);
        atom->masks_4xn_fn[cond0 * 4 + 1] = (unsigned int)(0xf - 0x2 * cond0);
        atom->masks_4xn_fn[cond0 * 4 + 2] = (unsigned int)(0xf - 0x4 * cond0);
        atom->masks_4xn_fn[cond0 * 4 + 3] = (unsigned int)(0xf - 0x8 * cond0);
    }
#else
    for (unsigned int cond0 = 0; cond0 < 2; cond0++) {
        for (unsigned int cond1 = 0; cond1 < 2; cond1++) {
#if CLUSTER_M < CLUSTER_N
            mask0 = (unsigned int)(0xff - 0x1 * cond0 - 0x1f * cond1);
            mask1 = (unsigned int)(0xff - 0x3 * cond0 - 0x3f * cond1);
            mask2 = (unsigned int)(0xff - 0x7 * cond0 - 0x7f * cond1);
            mask3 = (unsigned int)(0xff - 0xf * cond0 - 0xff * cond1);
#else
            mask0 = (unsigned int)(0x3 - 0x1 * cond0);
            mask1 = (unsigned int)(0x3 - 0x3 * cond0);
            mask2 = (unsigned int)(0x3 - cond0 * 0x3 - 0x1 * cond1);
            mask3 = (unsigned int)(0x3 - cond0 * 0x3 - 0x3 * cond1);
#endif

            atom->masks_2xnn_hn[cond0 * 4 + cond1 * 2 + 0] = (mask1 << halfMaskBits) |
                                                             mask0;
            atom->masks_2xnn_hn[cond0 * 4 + cond1 * 2 + 1] = (mask3 << halfMaskBits) |
                                                             mask2;

#if CLUSTER_M < CLUSTER_N
            mask0 = (unsigned int)(0xff - 0x1 * cond0 - 0x10 * cond1);
            mask1 = (unsigned int)(0xff - 0x2 * cond0 - 0x20 * cond1);
            mask2 = (unsigned int)(0xff - 0x4 * cond0 - 0x40 * cond1);
            mask3 = (unsigned int)(0xff - 0x8 * cond0 - 0x80 * cond1);
#else
            mask0 = (unsigned int)(0x3 - 0x1 * cond0);
            mask1 = (unsigned int)(0x3 - 0x2 * cond0);
            mask2 = (unsigned int)(0x3 - 0x1 * cond1);
            mask3 = (unsigned int)(0x3 - 0x2 * cond1);
#endif

            atom->masks_2xnn_fn[cond0 * 4 + cond1 * 2 + 0] = (mask1 << halfMaskBits) |
                                                             mask0;
            atom->masks_2xnn_fn[cond0 * 4 + cond1 * 2 + 1] = (mask3 << halfMaskBits) |
                                                             mask2;

#if CLUSTER_M < CLUSTER_N
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 0] = (unsigned int)(0xff -
                                                                           0x1 * cond0 -
                                                                           0x1f * cond1);
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 1] = (unsigned int)(0xff -
                                                                           0x3 * cond0 -
                                                                           0x3f * cond1);
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 2] = (unsigned int)(0xff -
                                                                           0x7 * cond0 -
                                                                           0x7f * cond1);
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 3] = (unsigned int)(0xff -
                                                                           0xf * cond0 -
                                                                           0xff * cond1);

            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 0] = (unsigned int)(0xff -
                                                                           0x1 * cond0 -
                                                                           0x10 * cond1);
            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 1] = (unsigned int)(0xff -
                                                                           0x2 * cond0 -
                                                                           0x20 * cond1);
            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 2] = (unsigned int)(0xff -
                                                                           0x4 * cond0 -
                                                                           0x40 * cond1);
            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 3] = (unsigned int)(0xff -
                                                                           0x8 * cond0 -
                                                                           0x80 * cond1);
#else
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 0] = (unsigned int)(0x3 -
                                                                           0x1 * cond0);
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 1] = (unsigned int)(0x3 -
                                                                           0x3 * cond0);
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 2] = (unsigned int)(0x3 -
                                                                           0x3 * cond0 -
                                                                           0x1 * cond1);
            atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 3] = (unsigned int)(0x3 -
                                                                           0x3 * cond0 -
                                                                           0x3 * cond1);

            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 0] = (unsigned int)(0x3 -
                                                                           0x1 * cond0);
            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 1] = (unsigned int)(0x3 -
                                                                           0x2 * cond0);
            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 2] = (unsigned int)(0x3 -
                                                                           0x1 * cond1);
            atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 3] = (unsigned int)(0x3 -
                                                                           0x2 * cond1);
#endif
        }
    }
#endif
}

void growAtom(Atom* atom) {
    int nold = atom->Nmax;
    atom->Nmax += DELTA;

#ifdef ATOM_POSITION_AOS
    atom->x = (MD_FLOAT*)reallocate(atom->x,
        ALIGNMENT,
        atom->Nmax * sizeof(MD_FLOAT) * 3,
        nold * sizeof(MD_FLOAT) * 3);
#else
    atom->x = (MD_FLOAT*)reallocate(atom->x,
        ALIGNMENT,
        atom->Nmax * sizeof(MD_FLOAT),
        nold * sizeof(MD_FLOAT));
    atom->y = (MD_FLOAT*)reallocate(atom->y,
        ALIGNMENT,
        atom->Nmax * sizeof(MD_FLOAT),
        nold * sizeof(MD_FLOAT));
    atom->z = (MD_FLOAT*)reallocate(atom->z,
        ALIGNMENT,
        atom->Nmax * sizeof(MD_FLOAT),
        nold * sizeof(MD_FLOAT));
#endif
    atom->vx   = (MD_FLOAT*)reallocate(atom->vx,
        ALIGNMENT,
        atom->Nmax * sizeof(MD_FLOAT),
        nold * sizeof(MD_FLOAT));
    atom->vy   = (MD_FLOAT*)reallocate(atom->vy,
        ALIGNMENT,
        atom->Nmax * sizeof(MD_FLOAT),
        nold * sizeof(MD_FLOAT));
    atom->vz   = (MD_FLOAT*)reallocate(atom->vz,
        ALIGNMENT,
        atom->Nmax * sizeof(MD_FLOAT),
        nold * sizeof(MD_FLOAT));
    atom->type = (int*)
        reallocate(atom->type, ALIGNMENT, atom->Nmax * sizeof(int), nold * sizeof(int));
}

void growClusters(Atom* atom, int super_clustering) {
    int nold  = atom->Nclusters_max;
    // If M<N, we need to allocate more j-clusters
    int jterm = MAX(1, CLUSTER_N / CLUSTER_M);

    atom->Nclusters_max += DELTA;
    atom->iclusters    = (Cluster*)reallocate(atom->iclusters,
        ALIGNMENT,
        atom->Nclusters_max * SCLUSTER_SIZE * sizeof(Cluster),
        nold * SCLUSTER_SIZE * sizeof(Cluster));
    atom->jclusters    = (Cluster*)reallocate(atom->jclusters,
        ALIGNMENT,
        atom->Nclusters_max * SCLUSTER_SIZE * jterm * sizeof(Cluster),
        nold * SCLUSTER_SIZE * jterm * sizeof(Cluster));
    atom->cluster_bin = (int*)reallocate(atom->cluster_bin,
        ALIGNMENT,
        atom->Nclusters_max * sizeof(int),
        nold * sizeof(int));
    atom->cl_x         = (MD_FLOAT*)reallocate(atom->cl_x,
        ALIGNMENT,
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT),
        nold * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
    atom->cl_f         = (MD_FLOAT*)reallocate(atom->cl_f,
        ALIGNMENT,
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE  * 3 * sizeof(MD_FLOAT),
        nold * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
    atom->cl_v         = (MD_FLOAT*)reallocate(atom->cl_v,
        ALIGNMENT,
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE  * 3 * sizeof(MD_FLOAT),
        nold * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
    atom->cl_t         = (int*)reallocate(atom->cl_t,
        ALIGNMENT,
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE  * sizeof(int),
        nold * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
    atom->cl_sqrt_epsilon = (MD_FLOAT*)reallocate(atom->cl_sqrt_epsilon,
        ALIGNMENT,
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT),
        nold * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
    atom->cl_sigma3       = (MD_FLOAT*)reallocate(atom->cl_sigma3,
        ALIGNMENT,
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT),
        nold * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));

    // NUMA-aware first-touch initialization for newly allocated memory
    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp for schedule(runtime)
        for (int ci = nold; ci < atom->Nclusters_max; ci++) {
            int ci_vec_base = ci * CLUSTER_M * SCLUSTER_SIZE * 3;
            for (int idx = 0; idx < CLUSTER_M * SCLUSTER_SIZE * 3; idx++) {
                atom->cl_x[ci_vec_base + idx] = 0.0;
                atom->cl_f[ci_vec_base + idx] = 0.0;
                atom->cl_v[ci_vec_base + idx] = 0.0;
            }
        }

        #pragma omp for schedule(runtime)
        for (int ci = nold; ci < atom->Nclusters_max; ci++) {
            int ci_sca_base = ci * CLUSTER_M * SCLUSTER_SIZE;
            for (int idx = 0; idx < CLUSTER_M * SCLUSTER_SIZE; idx++) {
                atom->cl_t[ci_sca_base + idx] = 0;
            }
        }

        #pragma omp for schedule(runtime)
        for (int ci = nold; ci < atom->Nclusters_max; ci++) {
            atom->cluster_bin[ci] = 0;
            atom->iclusters[ci * SCLUSTER_SIZE].natoms = 0;
            atom->iclusters[ci * SCLUSTER_SIZE].bbminx = 0.0;
            atom->iclusters[ci * SCLUSTER_SIZE].bbmaxx = 0.0;
            atom->iclusters[ci * SCLUSTER_SIZE].bbminy = 0.0;
            atom->iclusters[ci * SCLUSTER_SIZE].bbmaxy = 0.0;
            atom->iclusters[ci * SCLUSTER_SIZE].bbminz = 0.0;
            atom->iclusters[ci * SCLUSTER_SIZE].bbmaxz = 0.0;
        }
    }
    #endif

    if(super_clustering) {
        atom->siclusters    = (SuperCluster*)reallocate(atom->siclusters,
            ALIGNMENT,
            atom->Nclusters_max * sizeof(SuperCluster),
            nold * sizeof(SuperCluster));
    }

#ifdef CUDA_TARGET
    growClustersCUDA(atom);
#endif
}

/* MPI added*/

void freeAtom(Atom* atom) {
#ifdef ATOM_POSITION_AOS
    free(atom->x);
    atom->x = NULL;
#else
    free(atom->x);
    atom->x = NULL;
    free(atom->y);
    atom->y = NULL;
    free(atom->z);
    atom->z = NULL;
#endif
    free(atom->vx);
    atom->vx = NULL;
    free(atom->vy);
    atom->vy = NULL;
    free(atom->vz);
    atom->vz = NULL;
    free(atom->type);
    atom->type = NULL;
}

void growPbc(Atom* atom) {
    int nold = atom->NmaxGhost;
    atom->NmaxGhost += DELTA;

    if (atom->PBCx || atom->PBCy || atom->PBCz) {
        atom->PBCx = (int*)reallocate(atom->PBCx,
            ALIGNMENT,
            atom->NmaxGhost * sizeof(int),
            nold * sizeof(int));
        atom->PBCy = (int*)reallocate(atom->PBCy,
            ALIGNMENT,
            atom->NmaxGhost * sizeof(int),
            nold * sizeof(int));
        atom->PBCz = (int*)reallocate(atom->PBCz,
            ALIGNMENT,
            atom->NmaxGhost * sizeof(int),
            nold * sizeof(int));
    } else {
        atom->PBCx = (int*)malloc(atom->NmaxGhost * sizeof(int));
        atom->PBCy = (int*)malloc(atom->NmaxGhost * sizeof(int));
        atom->PBCz = (int*)malloc(atom->NmaxGhost * sizeof(int));
    }
}

void packForward(Atom* atom, int nc, int* list, MD_FLOAT* buf, int* pbc) {
    for (int i = 0; i < nc; i++) {
        int cj          = list[i];
        int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
        MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
        int displ       = i * CLUSTER_N;

        for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
            buf[3 * (displ + cjj) + 0] = cj_x[CL_X_INDEX(cjj)] +
                                         pbc[0] * atom->mybox.xprd;
            buf[3 * (displ + cjj) + 1] = cj_x[CL_Y_INDEX(cjj)] +
                                         pbc[1] * atom->mybox.yprd;
            buf[3 * (displ + cjj) + 2] = cj_x[CL_Z_INDEX(cjj)] +
                                         pbc[2] * atom->mybox.zprd;
        }

        for (int cjj = atom->jclusters[cj].natoms; cjj < CLUSTER_N; cjj++) {
            buf[3 * (displ + cjj) + 0] = -1; // x
            buf[3 * (displ + cjj) + 1] = -1; // y
            buf[3 * (displ + cjj) + 2] = -1; // z
        }
    }
}

void unpackForward(Atom* atom, int nc, int c0, MD_FLOAT* buf) {
    for (int i = 0; i < nc; i++) {
        int cj          = c0 + i;
        int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
        MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
        int displ       = i * CLUSTER_N;

        for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
            if (cj_x[CL_X_INDEX(cjj)] < INF)
                cj_x[CL_X_INDEX(cjj)] = buf[3 * (displ + cjj) + 0];
            if (cj_x[CL_Y_INDEX(cjj)] < INF)
                cj_x[CL_Y_INDEX(cjj)] = buf[3 * (displ + cjj) + 1];
            if (cj_x[CL_Z_INDEX(cjj)] < INF)
                cj_x[CL_Z_INDEX(cjj)] = buf[3 * (displ + cjj) + 2];
        }
    }
}

int packGhost(Atom* atom, int cj, MD_FLOAT* buf, int* pbc) {
    // # of elements per cluster natoms,x0,y0,z0,type_0, . .
    // ,xn,yn,zn,type_n,bbminx,bbmaxxy,bbminy,bbmaxy,bbminz,bbmaxz count = 4*N_CLUSTER+7,
    // if N_CLUSTER =4 => count = 23 value/cluster + trackpbc[x] + trackpbc[y] +
    // trackpbc[z]
    int m = 0;
    if (atom->jclusters[cj].natoms > 0) {
        int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
        int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
        MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
        MD_FLOAT bbminx = INF, bbmaxx = -INF;
        MD_FLOAT bbminy = INF, bbmaxy = -INF;
        MD_FLOAT bbminz = INF, bbmaxz = -INF;

        buf[m++] = (MD_FLOAT) atom->jclusters[cj].natoms;

        for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
            MD_FLOAT xtmp = cj_x[CL_X_INDEX(cjj)] + pbc[0] * atom->mybox.xprd;
            MD_FLOAT ytmp = cj_x[CL_Y_INDEX(cjj)] + pbc[1] * atom->mybox.yprd;
            MD_FLOAT ztmp = cj_x[CL_Z_INDEX(cjj)] + pbc[2] * atom->mybox.zprd;

            buf[m++] = xtmp;
            buf[m++] = ytmp;
            buf[m++] = ztmp;
            buf[m++] = (MD_FLOAT)atom->cl_t[cj_sca_base + cjj];

            if (bbminx > xtmp) {
                bbminx = xtmp;
            }
            if (bbmaxx < xtmp) {
                bbmaxx = xtmp;
            }
            if (bbminy > ytmp) {
                bbminy = ytmp;
            }
            if (bbmaxy < ytmp) {
                bbmaxy = ytmp;
            }
            if (bbminz > ztmp) {
                bbminz = ztmp;
            }
            if (bbmaxz < ztmp) {
                bbmaxz = ztmp;
            }
        }

        for (int cjj = atom->jclusters[cj].natoms; cjj < CLUSTER_N; cjj++) {
            buf[m++] = -1.0; // x
            buf[m++] = -1.0; // y
            buf[m++] = -1.0; // z
            buf[m++] = -1.0; // type
        }

        buf[m++] = bbminx;
        buf[m++] = bbmaxx;
        buf[m++] = bbminy;
        buf[m++] = bbmaxy;
        buf[m++] = bbminz;
        buf[m++] = bbmaxz;
        // TODO: check atom->ncj
        int ghostId = cj - atom->ncj;
        // check for ghost particles
        buf[m++] = (MD_FLOAT)(cj - atom->ncj >= 0) ? pbc[0] + atom->PBCx[ghostId]
                                                   : pbc[0];
        buf[m++] = (MD_FLOAT)(cj - atom->ncj >= 0) ? pbc[1] + atom->PBCy[ghostId]
                                                   : pbc[1];
        buf[m++] = (MD_FLOAT)(cj - atom->ncj >= 0) ? pbc[2] + atom->PBCz[ghostId]
                                                   : pbc[2];
    }

    return m;
}

int unpackGhost(Parameter *param, Atom* atom, int cj, MD_FLOAT* buf) {
    int m    = 0;
    int jfac = MAX(1, CLUSTER_N / CLUSTER_M);
    if (cj * jfac >= atom->Nclusters_max) {
        growClusters(atom, param->super_clustering);
    }

    if (atom->Nclusters_ghost >= atom->NmaxGhost) {
        growPbc(atom);
    }

    int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
    int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
    MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];

    atom->jclusters[cj].natoms = (int)buf[m++];
    for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
        cj_x[CL_X_INDEX(cjj)]       = buf[m++];
        cj_x[CL_Y_INDEX(cjj)]       = buf[m++];
        cj_x[CL_Z_INDEX(cjj)]       = buf[m++];
        atom->cl_t[cj_sca_base + cjj] = (int)buf[m++];
        atom->Nghost++;
    }

    for (int cjj = atom->jclusters[cj].natoms; cjj < CLUSTER_N; cjj++) {
        cj_x[CL_X_INDEX(cjj)]       = INF;
        cj_x[CL_Y_INDEX(cjj)]       = INF;
        cj_x[CL_Z_INDEX(cjj)]       = INF;
        atom->cl_t[cj_sca_base + cjj] = 0;
        m += 4;
    }

    atom->jclusters[cj].bbminx        = buf[m++];
    atom->jclusters[cj].bbmaxx        = buf[m++];
    atom->jclusters[cj].bbminy        = buf[m++];
    atom->jclusters[cj].bbmaxy        = buf[m++];
    atom->jclusters[cj].bbminz        = buf[m++];
    atom->jclusters[cj].bbmaxz        = buf[m++];
    atom->PBCx[atom->Nclusters_ghost] = (int)buf[m++];
    atom->PBCy[atom->Nclusters_ghost] = (int)buf[m++];
    atom->PBCz[atom->Nclusters_ghost] = (int)buf[m++];
    atom->Nclusters_ghost++;
    return m;
}

void packReverse(Atom* atom, int nc, int c0, MD_FLOAT* buf) {
    for (int i = 0; i < nc; i++) {
        int cj          = c0 + i;
        int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
        MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];
        int displ       = i * CLUSTER_N;

        for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
            buf[3 * (displ + cjj) + 0] = cj_f[CL_X_INDEX_3D(cjj)];
            buf[3 * (displ + cjj) + 1] = cj_f[CL_Y_INDEX_3D(cjj)];
            buf[3 * (displ + cjj) + 2] = cj_f[CL_Z_INDEX_3D(cjj)];
        }

        for (int cjj = atom->jclusters[cj].natoms; cjj < CLUSTER_N; cjj++) {
            buf[3 * (displ + cjj) + 0] = -1; // x
            buf[3 * (displ + cjj) + 1] = -1; // y
            buf[3 * (displ + cjj) + 2] = -1; // z
        }
    }
}

void unpackReverse(Atom* atom, int nc, int* list, MD_FLOAT* buf) {
    for (int i = 0; i < nc; i++) {
        int cj          = list[i];
        int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
        MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];
        int displ       = i * CLUSTER_N;

        for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
            cj_f[CL_X_INDEX_3D(cjj)] += buf[3 * (displ + cjj) + 0];
            cj_f[CL_Y_INDEX_3D(cjj)] += buf[3 * (displ + cjj) + 1];
            cj_f[CL_Z_INDEX_3D(cjj)] += buf[3 * (displ + cjj) + 2];
        }
    }
}

int packExchange(Atom* atom, int i, MD_FLOAT* buf) {
    int m    = 0;
    buf[m++] = atom_x(i);
    buf[m++] = atom_y(i);
    buf[m++] = atom_z(i);
    buf[m++] = atom_vx(i);
    buf[m++] = atom_vy(i);
    buf[m++] = atom_vz(i);
    buf[m++] = atom->type[i];
    return m;
}

int unpackExchange(Atom* atom, int i, MD_FLOAT* buf) {
    while (i >= atom->Nmax) {
        growAtom(atom);
    }

    int m         = 0;
    atom_x(i)     = buf[m++];
    atom_y(i)     = buf[m++];
    atom_z(i)     = buf[m++];
    atom_vx(i)    = buf[m++];
    atom_vy(i)    = buf[m++];
    atom_vz(i)    = buf[m++];
    atom->type[i] = buf[m++];
    return m;
}

void pbc(Atom* atom)
{
    for (int i = 0; i < atom->Nlocal; i++) {
        MD_FLOAT xprd = atom->mybox.xprd;
        MD_FLOAT yprd = atom->mybox.yprd;
        MD_FLOAT zprd = atom->mybox.zprd;
        if (atom_x(i) < 0.0) atom_x(i) += xprd;
        if (atom_y(i) < 0.0) atom_y(i) += yprd;
        if (atom_z(i) < 0.0) atom_z(i) += zprd;
        if (atom_x(i) >= xprd) atom_x(i) -= xprd;
        if (atom_y(i) >= yprd) atom_y(i) -= yprd;
        if (atom_z(i) >= zprd) atom_z(i) -= zprd;
    }
}

void copy(Atom* atom, int i, int j)
{
    atom_x(i)     = atom_x(j);
    atom_y(i)     = atom_y(j);
    atom_z(i)     = atom_z(j);
    atom_vx(i)    = atom_vx(j);
    atom_vy(i)    = atom_vy(j);
    atom_vz(i)    = atom_vz(j);
    atom->type[i] = atom->type[j];
}
