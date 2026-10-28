/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

/* Per-thread force output of the half-list kernels.
 *
 * A half list writes the reaction force of a pair into the j-cluster, which another
 * thread may own as an i-cluster or also write as a j-cluster (local or ghost). Like the
 * per-thread output buffers of GROMACS' NBNxM, every thread but thread 0 accumulates into
 * its own force array and the arrays are summed into cl_f at the end of the kernel.
 * Thread 0 writes cl_f directly, a single thread therefore runs without extra work.
 *
 * Usage: threadForcesReserve() before the parallel region, threadForces() at its
 * start (all threads) and threadForcesReduce() at its end (all threads).
 */

#ifndef __THREAD_FORCES_H_
#define __THREAD_FORCES_H_

#include <allocate.h>
#include <atom.h>
#include <force.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* End of the part of cl_f a kernel writes: the local i-clusters, the ghost j-clusters and
 * the dummy j-cluster that pads the neighbor lists */
static inline int threadForcesEnd(const Atom* atom)
{
    int ncj_all = atom->ncj + atom->Nclusters_ghost;
    if (atom->dummy_cj >= ncj_all) {
        ncj_all = atom->dummy_cj + 1;
    }
#if CLUSTER_M == CLUSTER_N * 2
    // Two j-clusters share the coordinate blocks of one i-cluster
    ncj_all = (ncj_all + 1) & ~0x1;
#endif
    int iend = atom->Nclusters_local;
#if CLUSTER_M == CLUSTER_N / 2
    // Two i-clusters share the coordinate blocks of one j-cluster
    iend = (iend + 1) & ~0x1;
#endif
    const int fi = CI_VECTOR3_BASE_INDEX(iend);
    const int fj = CJ_VECTOR3_BASE_INDEX(ncj_all);
    return fi > fj ? fi : fj;
}

static inline int threadForcesMaxThreads(void)
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/* Make room for the force arrays of threads 1..n-1, returns their length */
static inline int threadForcesReserve(Atom* atom)
{
    const int nf        = threadForcesEnd(atom);
    const size_t needed = (size_t)(threadForcesMaxThreads() - 1) * (size_t)nf;
    if (needed > atom->cl_f_threads_size) {
        free(atom->cl_f_threads);
        atom->cl_f_threads = (MD_FLOAT*)allocate(ALIGNMENT, needed * sizeof(MD_FLOAT));
        atom->cl_f_threads_size = needed;
    }
    return nf;
}

/* Force array of the calling thread, zeroed (cl_f for thread 0, which zeroForces()
 * clears before the parallel region) */
static inline MD_FLOAT* threadForces(Atom* atom, const int nf)
{
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
    if (tid > 0) {
        MD_FLOAT* f = &atom->cl_f_threads[(size_t)(tid - 1) * nf];
        memset(f, 0, (size_t)nf * sizeof(MD_FLOAT));
        return f;
    }
#endif
    (void)nf;
    return atom->cl_f;
}

/* Sum the force arrays of the other threads into cl_f */
static inline void threadForcesReduce(Atom* atom, const int nf)
{
#ifdef _OPENMP
    const int nthreads = omp_get_num_threads();
    if (nthreads > 1) {
#pragma omp barrier
        MD_FLOAT* restrict f        = atom->cl_f;
        const MD_FLOAT* restrict tf = atom->cl_f_threads;

#pragma omp for schedule(static)
        for (int i = 0; i < nf; i++) {
            MD_FLOAT sum = f[i];

            for (int t = 0; t < nthreads - 1; t++) {
                sum += tf[(size_t)t * nf + i];
            }

            f[i] = sum;
        }
    }
#else
    (void)atom;
    (void)nf;
#endif
}

#endif
