/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>

#include <atom.h>
#include <parameter.h>
#include <stats.h>
#include <timers.h>

#ifdef _MPI
#include <mpi.h>
#endif  

void initStats(Stats* s)
{
    s->total_force_neighs   = 0;
    s->total_force_iters    = 0;
    s->atoms_within_cutoff  = 0;
    s->atoms_outside_cutoff = 0;
}

void displayStatistics(Atom* atom, Parameter* param, Stats* stats, double* timer)
{
    int me = 0;
    long long int neigh_sum = 0;
    long long int forces_sum = 0;
    long long int iters_sum = 0;

#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Reduce(&stats->total_force_neighs, &neigh_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&stats->total_force_iters, &iters_sum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    stats->total_force_neighs = neigh_sum;
    stats->total_force_iters = iters_sum;
#endif

#ifdef COMPUTE_STATS

    double force_useful_volume = 1e-9 * ((double)(atom->Natoms * (param->ntimes + 1)) *
                                                (sizeof(MD_FLOAT) * 6 + sizeof(int)) +
                                            (double)(stats->total_force_neighs) *
                                                (sizeof(MD_FLOAT) * 3 + sizeof(int)));
    double avg_neigh           = stats->total_force_neighs /
                       (double)(atom->Natoms * (param->ntimes + 1));
    double avg_simd = stats->total_force_iters /
                      (double)(atom->Natoms * (param->ntimes + 1));

#ifndef ONE_ATOM_TYPE
    force_useful_volume += 1e-9 *
                           (double)((atom->Natoms * (param->ntimes + 1)) +
                                    stats->total_force_neighs) *
                           sizeof(int);
#endif
    if (me == 0) {
        printf("Statistics:\n");
        printf("\tVector width: %d, Processor frequency: %.4f GHz\n",
            VECTOR_WIDTH,
            param->proc_freq);
        printf("\tAverage neighbors per atom: %.4f\n", avg_neigh);
        printf("\tAverage SIMD iterations per atom: %.4f\n", avg_simd);
        printf("\tTotal number of computed pair interactions: %lld\n",
            stats->total_force_neighs);
        printf("\tTotal number of SIMD iterations: %lld\n", stats->total_force_iters);
        printf("\tUseful read data volume for force computation: %.2fGB\n",
            force_useful_volume);
        printf("\tCycles/SIMD iteration: %.4f\n",
            timer[FORCE] * param->proc_freq * 1e9 / stats->total_force_iters);
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

    const double eff_pct = (double)stats->atoms_within_cutoff /
                           (double)(stats->atoms_within_cutoff +
                                    stats->atoms_outside_cutoff) *
                           100.0;
    if (me == 0) {
        printf("\tAtoms within/outside cutoff radius: %lld/%lld (%.2f%%)\n",
            stats->atoms_within_cutoff,
            stats->atoms_outside_cutoff,
            eff_pct);
    }
#endif

#endif
}
