/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>

#include <atom.h>
#include <force.h>
#include <parameter.h>
#include <stats.h>
#include <timers.h>

#ifdef _MPI
#include <mpi.h>
#endif

void initStats(Stats* s) {
    s->calculated_forces       = 0;
    s->num_neighs              = 0;
    s->force_iters             = 0;
    s->atoms_within_cutoff     = 0;
    s->atoms_outside_cutoff    = 0;
    s->clusters_within_cutoff  = 0;
    s->clusters_outside_cutoff = 0;
}

void displayStatistics(Atom* atom, Parameter* param, Stats* stats, double* timer)
{
    int me = 0;
    long long int neigh_sum = 0;
    long long int forces_sum = 0;
    long long int iters_sum = 0;
    long long int Ncluster_sum = 0;

#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Reduce(&stats->num_neighs, &neigh_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats->calculated_forces, &forces_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats->force_iters, &iters_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&atom->Nclusters_local, &Ncluster_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    stats->num_neighs = neigh_sum;
    stats->calculated_forces = forces_sum;
    stats->force_iters = iters_sum;
    atom->Nclusters_local = Ncluster_sum;

#endif

#ifdef COMPUTE_STATS

    const int MxN            = CLUSTER_M * CLUSTER_N;
    double avgAtomsCluster   = (double)(atom->Natoms) / (double)(atom->Nclusters_local);
    double forceUsefulVolume = 1e-9 * ((double)(atom->Natoms * (param->ntimes + 1)) *
                                              (sizeof(MD_FLOAT) * 6 + sizeof(int)) +
                                          (double)(stats->num_neighs) *
                                              (sizeof(MD_FLOAT) * ATOM_DIM + sizeof(int)));
    double avgNeighAtom      = (stats->num_neighs * CLUSTER_N) /
                          (double)(atom->Natoms * (param->ntimes + 1));
    double avgNeighCluster = (double)(stats->num_neighs) /
                             (double)(stats->calculated_forces);
    double avgSimd = stats->force_iters / (double)(atom->Natoms * (param->ntimes + 1));

#ifndef ONE_ATOM_TYPE
    forceUsefulVolume += 1e-9 *
                         (double)((atom->Natoms * (param->ntimes + 1)) +
                                  stats->num_neighs) *
                         sizeof(int);
#endif
    if (me == 0) {
        printf("Statistics:\n");
        printf("\tVector width: %d, Processor frequency: %.4f GHz\n",
            VECTOR_WIDTH,
            param->proc_freq);
        printf("\tAverage atoms per cluster: %.4f\n", avgAtomsCluster);
        printf("\tAverage neighbors per atom: %.4f\n", avgNeighAtom);
        printf("\tAverage neighbors per cluster: %.4f\n", avgNeighCluster);
        printf("\tAverage SIMD iterations per atom: %.4f\n", avgSimd);
        printf("\tTotal number of computed pair interactions: %lld\n",
            stats->num_neighs * MxN);
        printf("\tTotal number of SIMD iterations: %lld\n", stats->force_iters);
        printf("\tUseful read data volume for force computation: %.2fGB\n",
            forceUsefulVolume);
        printf("\tCycles/SIMD iteration: %.4f\n",
            timer[FORCE] * param->proc_freq * 1e9 / stats->force_iters);
    }

#ifdef USE_REFERENCE_KERNEL

    long long int within_cutoff_sum = 0;
    long long int outside_cutoff_sum = 0;

#ifdef _MPI
    MPI_Reduce(&stats->atoms_within_cutoff, &within_cutoff_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats->atoms_outside_cutoff, &outside_cutoff_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    stats->atoms_within_cutoff = within_cutoff_sum;
    stats->atoms_outside_cutoff = outside_cutoff_sum;
#endif

    const double atoms_eff = (double)stats->atoms_within_cutoff /
                             (double)(stats->atoms_within_cutoff +
                                      stats->atoms_outside_cutoff) *
                             100.0;
    if (me == 0) {
        printf("\tAtoms within/outside cutoff radius: %lld/%lld (%.2f%%)\n",
            stats->atoms_within_cutoff,
            stats->atoms_outside_cutoff,
            atoms_eff);
    }
    const double clusters_eff = (double)stats->clusters_within_cutoff /
                                (double)(stats->clusters_within_cutoff +
                                         stats->clusters_outside_cutoff) *
                                100.0;
    if (me == 0) {
        printf("\tClusters within/outside cutoff radius: %lld/%lld (%.2f%%)\n",
            stats->clusters_within_cutoff,
            stats->clusters_outside_cutoff,
            clusters_eff);
    }
#endif

#endif
}
