/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>

#include <atom.h>
#include <force.h>
#include <util.h>
#include <string.h>
#include <vtk.h>

#ifdef _MPI
#include <mpi.h>
static MPI_File _fh;
static inline void flushBuffer(char*);
#endif

void write_data_to_vtk_file(const char* filename, Atom* atom, Parameter* param, int timestep) {
    write_local_atoms_to_vtk_file(filename, atom, param, timestep);
    write_ghost_atoms_to_vtk_file(filename, atom, timestep);
    write_local_cluster_edges_to_vtk_file(filename, atom, timestep);
    write_ghost_cluster_edges_to_vtk_file(filename, atom, timestep);
    write_super_clusters_to_vtk_file(filename, atom, timestep);
}

int write_super_clusters_to_vtk_file(const char* filename, Atom* atom, int timestep) {
    char timestep_filename[128];
    snprintf(timestep_filename,
        sizeof timestep_filename,
        "%s_sup_%d.vtk",
        filename,
        timestep);
    FILE* fp = fopen(timestep_filename, "wb");

    if (fp == NULL) {
        fprintf(stderr, "Could not open VTK file for writing!\n");
        return -1;
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Particle data\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d double\n", atom->Nclusters_local * SCLUSTER_SIZE * CLUSTER_M);

    for (int sci = 0; sci < atom->Nclusters_local; sci++) {
        int sci_vec_base = SCI_VECTOR_BASE_INDEX(sci);
        MD_FLOAT* sci_x  = &atom->cl_x[sci_vec_base];

        for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
            for (int cii = 0; cii < CLUSTER_M; cii++) {
                fprintf(fp,
                    "%.4f %.4f %.4f\n",
                    sci_x[CL_X_INDEX(sci_ci * CLUSTER_M + cii)],
                    sci_x[CL_Y_INDEX(sci_ci * CLUSTER_M + cii)],
                    sci_x[CL_Z_INDEX(sci_ci * CLUSTER_M + cii)]);
            }
        }
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "CELLS %d %d\n", atom->Nlocal, atom->Nlocal * 2);
    for (int i = 0; i < atom->Nlocal; ++i) {
        fprintf(fp, "1 %d\n", i);
    }
    fprintf(fp, "\n\n");
    fprintf(fp, "CELL_TYPES %d\n", atom->Nlocal);
    for (int i = 0; i < atom->Nlocal; ++i) {
        fprintf(fp, "1\n");
    }
    fprintf(fp, "\n\n");
    fprintf(fp, "POINT_DATA %d\n", atom->Nlocal);
    fprintf(fp, "SCALARS mass double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < atom->Nlocal; i++) {
        fprintf(fp, "1.0\n");
    }
    fprintf(fp, "\n\n");
    fclose(fp);
    return 0;
}

int write_local_atoms_to_vtk_file(const char* filename, Atom* atom, Parameter* param, int timestep) {
    char timestepFilename[128];
    snprintf(timestepFilename,
        sizeof timestepFilename,
        "%s_local_%d.vtk",
        filename,
        timestep);
    FILE* fp = fopen(timestepFilename, "wb");

    if (fp == NULL) {
        fprintf(stderr, "Could not open VTK file for writing!\n");
        return -1;
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Particle data\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d double\n", atom->Nlocal);

    if(param->super_clustering) {
        for (int sci = 0; sci < atom->Nclusters_local; ++sci) {
            int sciVecBase = SCI_VECTOR_BASE_INDEX(sci);
            MD_FLOAT* sciX = &atom->cl_x[sciVecBase];

            for (int sci_ci = 0; sci_ci < atom->siclusters[sci].nclusters; ++sci_ci) {
                const int ci = sci * SCLUSTER_SIZE + sci_ci;

                for (int cii = 0; cii < atom->iclusters[ci].natoms; ++cii) {
                    fprintf(fp,
                        "%.4f %.4f %.4f\n",
                        sciX[CL_X_INDEX(sci_ci * CLUSTER_M + cii)],
                        sciX[CL_Y_INDEX(sci_ci * CLUSTER_M + cii)],
                        sciX[CL_Z_INDEX(sci_ci * CLUSTER_M + cii)]);
                }
            }
        }
    } else {
        for (int ci = 0; ci < atom->Nclusters_local; ++ci) {
            int ciVecBase = CI_VECTOR_BASE_INDEX(ci);
            MD_FLOAT* ciX = &atom->cl_x[ciVecBase];
            for (int cii = 0; cii < atom->iclusters[ci].natoms; ++cii) {
                fprintf(fp,
                    "%.4f %.4f %.4f\n",
                    ciX[CL_X_INDEX(cii)],
                    ciX[CL_Y_INDEX(cii)],
                    ciX[CL_Z_INDEX(cii)]);
            }
        }
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "CELLS %d %d\n", atom->Nlocal, atom->Nlocal * 2);
    for (int i = 0; i < atom->Nlocal; ++i) {
        fprintf(fp, "1 %d\n", i);
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "CELL_TYPES %d\n", atom->Nlocal);
    for (int i = 0; i < atom->Nlocal; ++i) {
        fprintf(fp, "1\n");
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "POINT_DATA %d\n", atom->Nlocal);
    fprintf(fp, "SCALARS mass double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    for (int i = 0; i < atom->Nlocal; i++) {
        fprintf(fp, "1.0\n");
    }

    fprintf(fp, "\n\n");
    fclose(fp);
    return 0;
}

int write_ghost_atoms_to_vtk_file(const char* filename, Atom* atom, int timestep) {
    char timestepFilename[128];
    snprintf(timestepFilename,
        sizeof timestepFilename,
        "%s_ghost_%d.vtk",
        filename,
        timestep);
    FILE* fp = fopen(timestepFilename, "wb");

    if (fp == NULL) {
        fprintf(stderr, "Could not open VTK file for writing!\n");
        return -1;
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Particle data\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fp, "POINTS %d double\n", atom->Nghost);

    for (int ci = atom->Nclusters_local;
         ci < atom->Nclusters_local + atom->Nclusters_ghost;
         ++ci) {
        int ciVecBase = CI_VECTOR_BASE_INDEX(ci);
        MD_FLOAT* ciX = &atom->cl_x[ciVecBase];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; ++cii) {
            fprintf(fp,
                "%.4f %.4f %.4f\n",
                ciX[CL_X_INDEX(cii)],
                ciX[CL_Y_INDEX(cii)],
                ciX[CL_Z_INDEX(cii)]);
        }
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "CELLS %d %d\n", atom->Nghost, atom->Nghost * 2);

    for (int i = 0; i < atom->Nghost; ++i) {
        fprintf(fp, "1 %d\n", i);
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "CELL_TYPES %d\n", atom->Nghost);

    for (int i = 0; i < atom->Nghost; ++i) {
        fprintf(fp, "1\n");
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "POINT_DATA %d\n", atom->Nghost);
    fprintf(fp, "SCALARS mass double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    for (int i = 0; i < atom->Nghost; i++) {
        fprintf(fp, "1.0\n");
    }

    fprintf(fp, "\n\n");
    fclose(fp);
    return 0;
}

int write_local_cluster_edges_to_vtk_file(const char* filename, Atom* atom, int timestep) {
    char timestepFilename[128];
    snprintf(timestepFilename,
        sizeof timestepFilename,
        "%s_local_edges_%d.vtk",
        filename,
        timestep);
    FILE* fp     = fopen(timestepFilename, "wb");
    int N        = atom->Nclusters_local;
    int totLines = 0;
    int i        = 0;

    if (fp == NULL) {
        fprintf(stderr, "Could not open VTK file for writing!\n");
        return -1;
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Particle data\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");
    fprintf(fp, "POINTS %d double\n", atom->Nlocal);

    for (int ci = 0; ci < N; ++ci) {
        int ciVecBase = CI_VECTOR_BASE_INDEX(ci);
        MD_FLOAT* ciX = &atom->cl_x[ciVecBase];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; ++cii) {
            fprintf(fp,
                "%.4f %.4f %.4f\n",
                ciX[CL_X_INDEX(cii)],
                ciX[CL_Y_INDEX(cii)],
                ciX[CL_Z_INDEX(cii)]);
        }

        totLines += atom->iclusters[ci].natoms;
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "LINES %d %d\n", N, N + totLines);

    for (int ci = 0; ci < N; ++ci) {
        fprintf(fp, "%d ", atom->iclusters[ci].natoms);
        for (int cii = 0; cii < atom->iclusters[ci].natoms; ++cii) {
            fprintf(fp, "%d ", i++);
        }

        fprintf(fp, "\n");
    }

    fprintf(fp, "\n\n");
    fclose(fp);
    return 0;
}

int write_ghost_cluster_edges_to_vtk_file(const char* filename, Atom* atom, int timestep) {
    char timestepFilename[128];
    snprintf(timestepFilename,
        sizeof timestepFilename,
        "%s_ghost_edges_%d.vtk",
        filename,
        timestep);
    FILE* fp     = fopen(timestepFilename, "wb");
    int N        = atom->Nclusters_local + atom->Nclusters_ghost;
    int totLines = 0;
    int i        = 0;

    if (fp == NULL) {
        fprintf(stderr, "Could not open VTK file for writing!\n");
        return -1;
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Particle data\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");
    fprintf(fp, "POINTS %d double\n", atom->Nghost);

    for (int ci = atom->Nclusters_local; ci < N; ++ci) {
        int ciVecBase = CI_VECTOR_BASE_INDEX(ci);
        MD_FLOAT* ciX = &atom->cl_x[ciVecBase];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; ++cii) {
            fprintf(fp,
                "%.4f %.4f %.4f\n",
                ciX[CL_X_INDEX(cii)],
                ciX[CL_Y_INDEX(cii)],
                ciX[CL_Z_INDEX(cii)]);
        }

        totLines += atom->iclusters[ci].natoms;
    }

    fprintf(fp, "\n\n");
    fprintf(fp, "LINES %d %d\n", atom->Nclusters_ghost, atom->Nclusters_ghost + totLines);

    for (int ci = atom->Nclusters_local; ci < N; ++ci) {
        fprintf(fp, "%d ", atom->iclusters[ci].natoms);
        for (int cii = 0; cii < atom->iclusters[ci].natoms; ++cii) {
            fprintf(fp, "%d ", i++);
        }

        fprintf(fp, "\n");
    }

    fprintf(fp, "\n\n");
    fclose(fp);
    return 0;
}

#ifdef _MPI
int vtkOpen(const char* filename, Comm* comm, Atom* atom, int timestep) {
    char msg[256];
    char timestep_filename[128];

    snprintf(timestep_filename,
        sizeof timestep_filename,
        "%s_%d.vtk",
        filename,
        timestep);

    MPI_File_open(MPI_COMM_WORLD,
        timestep_filename,
        MPI_MODE_WRONLY | MPI_MODE_CREATE,
        MPI_INFO_NULL,
        &_fh);

    if (_fh == MPI_FILE_NULL) {
        fprintf_once(comm->myproc, stderr, "Could not open VTK file for writing!\n");
        return -1;
    }

    if (comm->myproc == 0) {
        sprintf(msg, "# vtk DataFile Version 2.0\n");
        sprintf(msg, "%sParticle data\n", msg);
        sprintf(msg, "%sASCII\n", msg);
        sprintf(msg, "%sDATASET UNSTRUCTURED_GRID\n", msg);
        sprintf(msg, "%sPOINTS %d double\n", msg, atom->Natoms);
        flushBuffer(msg);
    }
}

int vtkVector(Comm* comm, Atom* atom, Parameter* param) {
    if (_fh == MPI_FILE_NULL) {
        fprintf_once(comm->myproc, stderr, "VTK not initialize, call vtkOpen first!\n");
        return -1;
    }

    int sizeline  = 25; // #initial guess of characters in "%.4f %.4f %.4f\n"
    int extrabuff = 100;
    int sizebuff  = sizeline * atom->Nlocal + extrabuff;
    int mysize    = 0;
    char* msg     = (char*)malloc(sizebuff);

    sprintf(msg, "");

    for (int i = 0; i < atom->Nlocal; i++) {
        if (mysize + extrabuff >= sizebuff) {
            sizebuff *= 1.5;
            msg = (char*)realloc(msg, sizebuff);
        }
        // TODO: do not forget to add param->xlo, param->ylo, param->zlo
        sprintf(msg, "%s%.4f %.4f %.4f\n", msg, atom_x(i), atom_y(i), atom_z(i));
        mysize = strlen(msg);
    }

    int gatherSize[comm->numproc];
    int offset     = 0;
    int globalSize = 0;

    MPI_Allgather(&mysize, 1, MPI_INT, gatherSize, 1, MPI_INT, MPI_COMM_WORLD);

    for (int i = 0; i < comm->myproc; i++)
        offset += gatherSize[i];

    for (int i = 0; i < comm->numproc; i++)
        globalSize += gatherSize[i];

    MPI_Offset displ;
    MPI_Datatype FileType;
    int GlobalSize[] = { globalSize };
    int LocalSize[]  = { mysize };
    int Start[]      = { offset };

    if (LocalSize[0] > 0) {
        MPI_Type_create_subarray(1,
            GlobalSize,
            LocalSize,
            Start,
            MPI_ORDER_C,
            MPI_CHAR,
            &FileType);
    } else {
        MPI_Type_vector(0, 0, 0, MPI_CHAR, &FileType);
    }

    MPI_Type_commit(&FileType);
    MPI_File_get_size(_fh, &displ);
    MPI_File_set_view(_fh, displ, MPI_CHAR, FileType, "native", MPI_INFO_NULL);
    MPI_File_write_all(_fh, msg, mysize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_set_view(_fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

    if (comm->myproc == 0) {
        sprintf(msg, "\n\n");
        sprintf(msg, "%sCELLS %d %d\n", msg, atom->Natoms, atom->Natoms * 2);

        for (int i = 0; i < atom->Natoms; i++)
            sprintf(msg, "%s1 %d\n", msg, i);

        flushBuffer(msg);

        sprintf(msg, "\n\n");
        sprintf(msg, "%sCELL_TYPES %d\n", msg, atom->Natoms);

        for (int i = 0; i < atom->Natoms; i++)
            sprintf(msg, "%s1\n", msg);

        flushBuffer(msg);

        sprintf(msg, "\n\n");
        sprintf(msg, "%sPOINT_DATA %d\n", msg, atom->Natoms);
        sprintf(msg, "%sSCALARS mass double\n", msg);
        sprintf(msg, "%sLOOKUP_TABLE default\n", msg);

        for (int i = 0; i < atom->Natoms; i++)
            sprintf(msg, "%s1.0\n", msg);

        sprintf(msg, "%s\n\n", msg);
        flushBuffer(msg);
    }
}

void vtkClose() {
    MPI_File_close(&_fh);
    _fh = MPI_FILE_NULL;
}

static inline void flushBuffer(char* msg) {
    MPI_Offset displ;
    MPI_File_get_size(_fh, &displ);
    MPI_File_write_at(_fh, displ, msg, strlen(msg), MPI_CHAR, MPI_STATUS_IGNORE);
}

#endif

// TODO: print ghost and cluster using MPI
void printvtk(const char* filename, Comm* comm, Atom* atom, Parameter* param, int timestep) { 
#ifdef _MPI
    vtkOpen(filename, comm, atom, timestep);
    vtkVector(comm, atom, param);
    vtkClose();
#else  
    write_data_to_vtk_file(filename, atom, param, timestep);
#endif
}
