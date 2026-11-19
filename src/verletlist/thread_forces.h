/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

/* Per-thread force output of the half-list kernels, as in clusterpair/thread_forces.h.
 *
 * A half list writes the reaction force of a pair into atom j, which another thread may
 * own as an i-atom or also write as a j-atom (local or ghost). Every thread but thread 0
 * accumulates into its own force array and the arrays are summed into the atom forces at
 * the end of the kernel. Thread 0 writes the atom forces directly, a single thread
 * therefore runs without extra work.
 *
 * Usage: threadForcesReserve() before the parallel region, threadForces() at its
 * start (all threads) and threadForcesReduce() at its end (all threads).
 */

#ifndef __THREAD_FORCES_H_
#define __THREAD_FORCES_H_

#include <allocate.h>
#include <atom.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    MD_FLOAT* fx;
    MD_FLOAT* fy;
    MD_FLOAT* fz;
} ThreadForces;

static inline int threadForcesMaxThreads(void)
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

/* Make room for the force arrays of threads 1..n-1, returns the number of values of one
 * (3 per local and ghost atom) */
static inline int threadForcesReserve(Atom* atom)
{
    const int nf        = 3 * (atom->Nlocal + atom->Nghost);
    const size_t needed = (size_t)(threadForcesMaxThreads() - 1) * (size_t)nf;
    if (needed > atom->f_threads_size) {
        free(atom->f_threads);
        atom->f_threads      = (MD_FLOAT*)allocate(ALIGNMENT, needed * sizeof(MD_FLOAT));
        atom->f_threads_size = needed;
    }
    return nf;
}

/* Force arrays of the calling thread, zeroed (the atom forces for thread 0, which the
 * kernel clears before the parallel region). They are indexed like the atom forces, with
 * atom_access_x/y/z() */
static inline ThreadForces threadForces(Atom* atom, const int nf)
{
    ThreadForces tf = { atom->fx, atom->fy, atom->fz };
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
    if (tid > 0) {
        MD_FLOAT* f = &atom->f_threads[(size_t)(tid - 1) * nf];
        memset(f, 0, (size_t)nf * sizeof(MD_FLOAT));
#ifdef ATOM_POSITION_AOS
        tf.fx = tf.fy = tf.fz = f;
#else
        tf.fx = f;
        tf.fy = f + nf / 3;
        tf.fz = f + 2 * (nf / 3);
#endif
    }
#endif
    (void)nf;
    return tf;
}

/* Sum the force arrays of the other threads into the atom forces */
static inline void threadForcesReduce(Atom* atom, const int nf)
{
#ifdef _OPENMP
    const int nthreads = omp_get_num_threads();
    if (nthreads > 1) {
#pragma omp barrier
        const MD_FLOAT* restrict tf = atom->f_threads;

#ifdef ATOM_POSITION_AOS
        MD_FLOAT* restrict f = atom->fx;
#pragma omp for schedule(static)
        for (int i = 0; i < nf; i++) {
            MD_FLOAT sum = f[i];
            for (int t = 0; t < nthreads - 1; t++) {
                sum += tf[(size_t)t * nf + i];
            }
            f[i] = sum;
        }
#else
        const int n           = nf / 3;
        MD_FLOAT* restrict fx = atom->fx;
        MD_FLOAT* restrict fy = atom->fy;
        MD_FLOAT* restrict fz = atom->fz;
#pragma omp for schedule(static)
        for (int i = 0; i < n; i++) {
            MD_FLOAT sx = fx[i], sy = fy[i], sz = fz[i];
            for (int t = 0; t < nthreads - 1; t++) {
                const MD_FLOAT* f = &tf[(size_t)t * nf];
                sx += f[i];
                sy += f[n + i];
                sz += f[2 * n + i];
            }
            fx[i] = sx;
            fy[i] = sy;
            fz[i] = sz;
        }
#endif
    }
#else
    (void)atom;
    (void)nf;
#endif
}

#endif
