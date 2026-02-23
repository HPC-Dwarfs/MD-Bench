/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <string.h>
//---
#include <likwid-marker.h>
//---
#include <allocate.h>
#include <atom.h>
#include <eam.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <pbc.h>
#include <stats.h>
#include <thermo.h>
#include <timers.h>
#include <timing.h>
#include <util.h>

#define HLINE                                                                            \
    "----------------------------------------------------------------------------\n"

// Patterns
#define P_SEQ  0
#define P_FIX  1
#define P_RAND 2

void init(Parameter* param)
{
    param->input_file  = NULL;
    param->force_field = FF_LJ;
    param->epsilon     = 1.0;
    param->sigma6      = 1.0;
    param->rho         = 0.8442;
    param->ntypes      = 4;
    param->ntimes      = 200;
    param->nx          = 1;
    param->ny          = 1;
    param->nz          = 1;
    param->lattice     = 1.0;
    param->cutforce    = 1000000.0;
    param->cutneigh    = param->cutforce;
    param->mass        = 1.0;
    param->half_neigh  = 0;
    // Unused
    param->dt            = 0.005;
    param->dtforce       = 0.5 * param->dt;
    param->nstat         = 100;
    param->temp          = 1.44;
    param->reneigh_every = 20;
    param->proc_freq     = 2.4;
    param->eam_file      = NULL;
}

void createNeighbors(Atom* atom, Neighbor* neighbor, int pattern, int nneighs, int nreps)
{
    const int maxneighs = nneighs * nreps;
    neighbor->numneigh  = (int*)malloc(atom->Nmax * sizeof(int));
    neighbor->neighbors = (int*)malloc(atom->Nmax * maxneighs * sizeof(int));

    if (pattern == P_RAND && atom->Nlocal <= nneighs) {
        fprintf(stderr,
            "Error: When using random pattern, number of atoms should be higher than "
            "number of neighbors per atom!\n");
        exit(-1);
    }

    for (int i = 0; i < atom->Nlocal; i++) {
        int j         = (pattern == P_SEQ) ? (i + 1) : 0;
        int m         = (pattern == P_SEQ) ? atom->Nlocal : nneighs;

        for (int k = 0; k < nneighs; k++) {
            if (pattern == P_RAND) {
                int found = 0;
                do {
                    j           = rand() % atom->Nlocal;
                    neighptr[k] = j;
                    found       = (int)(i == j);
                    for (int l = 0; l < k; l++) {
                        if (neighptr[l] == j) {
                            found = 1;
                        }
                    }
                } while (found == 1);
            } else {
		neighs(neighbor->neighbors, i, k, atom->Nlocal, neighbor->maxneighs) = j;
                j           = (j + 1) % m;
            }
        }

        for (int r = 1; r < nreps; r++) {
            for (int k = 0; k < nneighs; k++) {
		neighs(neighbor->neighbors, i, r * nneighs + k, atom->Nlocal, neighbor->maxneighs) = 
		    neighs(neighbor->neighbors, i, k, atom->Nlocal, neighbor->maxneighs);
            }
        }

        neighbor->numneigh[i] = nneighs * nreps;
    }
}

int main(int argc, const char* argv[])
{
    Eam eam;
    Atom atom_data;
    Atom* atom = (Atom*)(&atom_data);
    Neighbor neighbor;
    Stats stats;
    Parameter param;
    char* pattern_str = NULL;
    int pattern       = P_SEQ;
    int natoms        = 256;
    int nneighs       = 76;
    int nreps         = 1;
    int csv           = 0;

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_REGISTER("force");
    DEBUG_MESSAGE("Initializing parameters...\n");
    init(&param);

    for (int i = 0; i < argc; i++) {
        if ((strcmp(argv[i], "-f") == 0)) {
            if ((param.force_field = str2ff(argv[++i])) < 0) {
                fprintf(stderr, "Invalid force field!\n");
                exit(-1);
            }
            continue;
        }
        if ((strcmp(argv[i], "-p") == 0)) {
            pattern_str = strdup(argv[++i]);
            if (strncmp(pattern_str, "seq", 3) == 0) {
                pattern = P_SEQ;
            } else if (strncmp(pattern_str, "fix", 3) == 0) {
                pattern = P_FIX;
            } else if (strncmp(pattern_str, "rand", 3) == 0) {
                pattern = P_RAND;
            } else {
                fprintf(stderr, "Invalid pattern!\n");
                exit(-1);
            }
            continue;
        }
        if ((strcmp(argv[i], "-e") == 0)) {
            param.eam_file = strdup(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-n") == 0) || (strcmp(argv[i], "--nsteps") == 0)) {
            param.ntimes = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-na") == 0)) {
            natoms = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-nn") == 0)) {
            nneighs = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-nr") == 0)) {
            nreps = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "--freq") == 0)) {
            param.proc_freq = atof(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "--csv") == 0)) {
            csv = 1;
            continue;
        }
        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            printf("MD Bench: A minimalistic re-implementation of miniMD\n");
            printf(HLINE);
            printf("-f <string>:          force field (lj or eam), default lj\n");
            printf(
                "-p <string>:          pattern for data accesses (seq, fix or rand)\n");
            printf("-n / --nsteps <int>:  number of timesteps for simulation\n");
            printf("-na <int>:            number of atoms (default 256)\n");
            printf("-nn <int>:            number of neighbors per atom (default 76)\n");
            printf("-nr <int>:            number of times neighbor lists should be "
                   "replicated (default 1)\n");
            printf("--freq <real>:        set CPU frequency (GHz) and display average "
                   "cycles per atom and neighbors\n");
            printf("--csv:                set output as CSV style\n");
            printf(HLINE);
            exit(EXIT_SUCCESS);
        }
    }

    if (pattern_str == NULL) {
        pattern_str = strdup("seq\0");
    }

    if (param.force_field == FF_EAM) {
        DEBUG_MESSAGE("Initializing EAM parameters...\n");
        initEam(&param);
    }

    DEBUG_MESSAGE("Initializing atoms...\n");
    initAtom(atom);
    initStats(&stats);

    atom->ntypes  = param.ntypes;
    atom->epsilon = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->sigma6  = allocate(ALIGNMENT, atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutforcesq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    atom->cutneighsq = allocate(ALIGNMENT,
        atom->ntypes * atom->ntypes * sizeof(MD_FLOAT));
    for (int i = 0; i < atom->ntypes * atom->ntypes; i++) {
        atom->epsilon[i]    = param.epsilon;
        atom->sigma6[i]     = param.sigma6;
        atom->cutneighsq[i] = param.cutneigh * param.cutneigh;
        atom->cutforcesq[i] = param.cutforce * param.cutforce;
    }

    DEBUG_MESSAGE("Creating atoms...\n");
    for (int i = 0; i < natoms; ++i) {
        while (atom->Nlocal > atom->Nmax - natoms) {
            growAtom(atom);
        }

        atom->type[atom->Nlocal] = rand() % atom->ntypes;
        atom_x(atom->Nlocal)     = (MD_FLOAT)(i)*0.00001;
        atom_y(atom->Nlocal)     = (MD_FLOAT)(i)*0.00001;
        atom_z(atom->Nlocal)     = (MD_FLOAT)(i)*0.00001;
        atom_vx(atom->Nlocal)    = 0.0;
        atom_vy(atom->Nlocal)    = 0.0;
        atom_vz(atom->Nlocal)    = 0.0;
        atom->Nlocal++;
    }

    const double estim_atom_volume      = (double)(atom->Nlocal * 3 * sizeof(MD_FLOAT));
    const double estim_neighbors_volume = (double)(atom->Nlocal * (nneighs + 2) *
                                                   sizeof(int));
    const double estim_volume           = (double)(atom->Nlocal * 6 * sizeof(MD_FLOAT) +
                                         estim_neighbors_volume);

    if (!csv) {
        printf("Pattern: %s\n", pattern_str);
        printf("Number of timesteps: %d\n", param.ntimes);
        printf("Number of atoms: %d\n", natoms);
        printf("Number of neighbors per atom: %d\n", nneighs);
        printf("Number of times to replicate neighbor lists: %d\n", nreps);
        printf("Estimated total data volume (kB): %.4f\n", estim_volume / 1000.0);
        printf("Estimated atom data volume (kB): %.4f\n", estim_atom_volume / 1000.0);
        printf("Estimated neighborlist data volume (kB): %.4f\n",
            estim_neighbors_volume / 1000.0);
    }

    DEBUG_MESSAGE("Initializing neighbor lists...\n");
    initNeighbor(&neighbor, &param);
    DEBUG_MESSAGE("Creating neighbor lists...\n");
    createNeighbors(atom, &neighbor, pattern, nneighs, nreps);
    DEBUG_MESSAGE("Computing forces...\n");

    double T_accum = 0.0;
    for (int i = 0; i < param.ntimes; i++) {
#if defined(MEM_TRACER) || defined(INDEX_TRACER)
        traceAddresses(&param, atom, &neighbor, i + 1);
#endif

        if (param.force_field == FF_EAM) {
            computeForceEam(&param, atom, &neighbor, &stats);
        } else {
            if (param.half_neigh) {
                T_accum += computeForceLJHalfNeigh(&param, atom, &neighbor, &stats);
            } else {
                T_accum += computeForceLJFullNeigh(&param, atom, &neighbor, &stats);
            }
        }
    }

    double freq_hz                     = param.proc_freq * 1.e9;
    const double atoms_updates_per_sec = (double)(atom->Nlocal) / T_accum *
                                         (double)(param.ntimes);
    const double cycles_per_atom = T_accum / (double)(atom->Nlocal) /
                                   (double)(param.ntimes) * freq_hz;
    const double cycles_per_neigh = cycles_per_atom / (double)(nneighs);

    if (!csv) {
        printf("Total time: %.4f, Mega atom updates/s: %.4f\n",
            T_accum,
            atoms_updates_per_sec / 1.e6);
        if (param.proc_freq > 0.0) {
            printf("Cycles per atom: %.4f, Cycles per neighbor: %.4f\n",
                cycles_per_atom,
                cycles_per_neigh);
        }
    } else {
        printf("steps,pattern,natoms,nneighs,nreps,total vol.(kB),atoms vol.(kB),neigh "
               "vol.(kB),time(s),atom upds/s(M)");
        if (param.proc_freq > 0.0) {
            printf(",cy/atom,cy/neigh");
        }
        printf("\n");

        printf("%d,%s,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f",
            param.ntimes,
            pattern_str,
            natoms,
            nneighs,
            nreps,
            estim_volume / 1.e3,
            estim_atom_volume / 1.e3,
            estim_neighbors_volume / 1.e3,
            T_accum,
            atoms_updates_per_sec / 1.e6);

        if (param.proc_freq > 0.0) {
            printf(",%.4f,%.4f", cycles_per_atom, cycles_per_neigh);
        }
        printf("\n");
    }

    double timer[NUMTIMER];
    timer[FORCE] = T_accum;
    displayStatistics(atom, &param, &stats, timer);
    LIKWID_MARKER_CLOSE;
    return EXIT_SUCCESS;
}
