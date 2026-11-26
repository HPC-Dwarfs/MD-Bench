/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// Verlet-list neighbor-list builds, a template on LJ_COMB_RULE: the Makefile compiles it
// once per rule and buildNeighborCPU() in neighbor.c calls the variant of the rule
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <allocate.h>
#include <atom.h>
#include <device.h>
#include <kernel_variant.h>
#include <neighbor.h>
#include <neighbor_internal.h>
#include <parameter.h>
#include <util.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// buildNeighborLists is the plain C build, buildNeighborListsSimd the intrinsics one
// (--neighbor-simd intrinsics). Under AOS/SOA the latter writes surviving neighbor
// indices through the neighs() macro directly into the padded per-atom slot; under CSR it
// appends them to that atom's growable per-thread buffer instead (see buildNeighborCPU)
// -- the vectorized bin-scan itself (gather candidate positions, build a lane mask,
// scalar bit-scan over survivors) is identical either way, only the storage step differs.
#ifdef HAVE_SIMD_NEIGHBOR
#include <simd.h>
#endif

#ifndef NBLIST_CSR

/* AOS/SOA, scalar bin scan. Layout only affects the neighs() macro used to
 * store each surviving index, so this one function covers both. */
void RULE_VARIANT_NAME(buildNeighborLists)(Atom* atom, Neighbor* neighbor)
{
    int resize = 1;

    while (resize) {
        int new_maxneighs = neighbor->maxneighs;
        int resize_local  = 0;

#pragma omp parallel for schedule(runtime) reduction(max : new_maxneighs)                \
    reduction(| : resize_local)
        for (int i = 0; i < atom->Nlocal; i++) {
            int n         = 0;
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            int ibin      = coord2bin(xtmp, ytmp, ztmp);
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int type_i = atom->type[i];
#endif
            for (int k = 0; k < nstencil; k++) {
                int jbin     = ibin + stencil[k];
                int* loc_bin = &bins[jbin * atoms_per_bin];

                for (int m = 0; m < bincount[jbin]; m++) {
                    int j = loc_bin[m];
                    if (i == j) continue;
                    if (neighbor->half_neigh && j < i) continue;
                    if (half_stencil && ibin == jbin && skipNeigh(atom, i, j)) continue;

                    MD_FLOAT delx = xtmp - atom_x(j);
                    MD_FLOAT dely = ytmp - atom_y(j);
                    MD_FLOAT delz = ztmp - atom_z(j);
                    MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                    int type_j = atom->type[j];
                    const MD_FLOAT cutoff =
                        atom->cutneighsq[type_i * atom->ntypes + type_j];
#else
                    const MD_FLOAT cutoff = cutneighsq;
#endif
                    if (rsq <= cutoff) {
                        // Past maxneighs only count, the list is resized and rebuilt
                        if (n < neighbor->maxneighs) {
                            neighs(neighbor->neighbors, i, n, atom->Nlocal, neighbor) = j;
                        }
                        n++;
                    }
                }
            }

            neighbor->numneigh[i] = n;
            if (n >= neighbor->maxneighs) {
                resize_local = 1;
                if (n > new_maxneighs) new_maxneighs = n;
            }
        }

        resize = resize_local;
        if (resize) {
            printf("RESIZE %d, PROC %d\n", neighbor->maxneighs, me);
            neighbor->maxneighs = new_maxneighs * 1.2;
            free(neighbor->neighbors);
            free(is_inner_buf);
            neighbor->neighbors = (int*)allocate(ALIGNMENT,
                atom->Nmax * neighbor->maxneighs * sizeof(int));
            memset(neighbor->neighbors,
                0,
                atom->Nmax * neighbor->maxneighs * sizeof(int));
            is_inner_buf = (int*)allocate(ALIGNMENT, neighbor->maxneighs * sizeof(int));
        }
    }
}

#endif // !NBLIST_CSR

#if !defined(NBLIST_CSR) && defined(HAVE_SIMD_NEIGHBOR)

/* AOS/SOA, vectorized bin scan: gather candidate positions, build a lane
 * mask from (i != j) & (half_neigh ? j >= i : true) & (rsq < cutoff), then
 * compact survivors via a scalar bit-scan. */
void RULE_VARIANT_NAME(buildNeighborListsSimd)(Atom* atom, Neighbor* neighbor)
{
    int resize = 1;

    while (resize) {
        int new_maxneighs = neighbor->maxneighs;
        int resize_local  = 0;

#pragma omp parallel for schedule(runtime) reduction(max : new_maxneighs)                \
    reduction(| : resize_local)
        for (int i = 0; i < atom->Nlocal; i++) {
            int n         = 0;
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            int ibin      = coord2bin(xtmp, ytmp, ztmp);
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int type_i = atom->type[i];
#endif
            MD_SIMD_INT ivec       = simd_i32_broadcast(i);
            MD_SIMD_FLOAT xtmp_vec = simd_real_broadcast(xtmp);
            MD_SIMD_FLOAT ytmp_vec = simd_real_broadcast(ytmp);
            MD_SIMD_FLOAT ztmp_vec = simd_real_broadcast(ztmp);
#if LJ_COMB_RULE == LJ_COMB_SINGLE
            MD_SIMD_FLOAT cutoff_vec_single = simd_real_broadcast(cutneighsq);
#else
            MD_SIMD_INT tbase_i = simd_i32_broadcast(type_i * atom->ntypes);
#endif
            for (int k = 0; k < nstencil; k++) {
                int jbin       = ibin + stencil[k];
                int* loc_bin   = &bins[jbin * atoms_per_bin];
                int bincount_j = bincount[jbin];

                /* skipNeigh() is an MPI half-stencil self-bin tie-break; fall
                 * back to scalar there rather than vectorizing that rare path. */
                if (half_stencil && ibin == jbin) {
                    for (int m = 0; m < bincount_j; m++) {
                        int j = loc_bin[m];
                        if (i == j) continue;
                        if (neighbor->half_neigh && j < i) continue;
                        if (skipNeigh(atom, i, j)) continue;

                        MD_FLOAT delx = xtmp - atom_x(j);
                        MD_FLOAT dely = ytmp - atom_y(j);
                        MD_FLOAT delz = ztmp - atom_z(j);
                        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                        int type_j = atom->type[j];
                        const MD_FLOAT cutoff =
                            atom->cutneighsq[type_i * atom->ntypes + type_j];
#else
                        const MD_FLOAT cutoff = cutneighsq;
#endif
                        if (rsq <= cutoff) {
                            if (n < neighbor->maxneighs) {
                                neighs(neighbor->neighbors,
                                    i,
                                    n,
                                    atom->Nlocal,
                                    neighbor) = j;
                            }
                            n++;
                        }
                    }
                    continue;
                }

                MD_SIMD_INT bincount_vec = simd_i32_broadcast(bincount_j);

                for (int m = 0; m < bincount_j; m += VECTOR_WIDTH) {
                    MD_SIMD_MASK bound_mask = simd_mask_i32_cond_lt(
                        simd_i32_add(simd_i32_broadcast(m), simd_i32_seq()),
                        bincount_vec);
                    MD_SIMD_INT jvec = simd_i32_mask_load(&loc_bin[m], bound_mask);

#ifdef ATOM_POSITION_AOS
                    MD_SIMD_INT j3   = simd_i32_add(simd_i32_add(jvec, jvec), jvec);
                    MD_SIMD_FLOAT xj = simd_real_gather(j3,
                        &(atom->x[0]),
                        sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT yj = simd_real_gather(j3,
                        &(atom->x[1]),
                        sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT zj = simd_real_gather(j3,
                        &(atom->x[2]),
                        sizeof(MD_FLOAT));
#else
                    MD_SIMD_FLOAT xj = simd_real_gather(jvec, atom->x, sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT yj = simd_real_gather(jvec, atom->y, sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT zj = simd_real_gather(jvec, atom->z, sizeof(MD_FLOAT));
#endif
                    MD_SIMD_FLOAT delx = simd_real_sub(xtmp_vec, xj);
                    MD_SIMD_FLOAT dely = simd_real_sub(ytmp_vec, yj);
                    MD_SIMD_FLOAT delz = simd_real_sub(ztmp_vec, zj);
                    MD_SIMD_FLOAT rsq  = simd_real_fma(delx,
                        delx,
                        simd_real_fma(dely, dely, simd_real_mul(delz, delz)));

                    MD_SIMD_MASK keep = simd_mask_and(bound_mask,
                        simd_mask_not(simd_mask_i32_cond_eq(jvec, ivec)));
                    if (neighbor->half_neigh) {
                        keep = simd_mask_and(keep,
                            simd_mask_not(simd_mask_i32_cond_lt(jvec, ivec)));
                    }

#if LJ_COMB_RULE == LJ_COMB_SINGLE
                    MD_SIMD_FLOAT cutoff_vec = cutoff_vec_single;
#else
                    MD_SIMD_INT type_j_vec   = simd_i32_gather(jvec,
                        atom->type,
                        sizeof(int));
                    MD_SIMD_INT tij          = simd_i32_add(tbase_i, type_j_vec);
                    MD_SIMD_FLOAT cutoff_vec = simd_real_gather(tij,
                        atom->cutneighsq,
                        sizeof(MD_FLOAT));
#endif
                    keep = simd_mask_and(keep, simd_mask_cond_lt(rsq, cutoff_vec));

                    unsigned int bits = simd_mask_to_u32(keep);
                    if (bits) {
                        int tmp_j[VECTOR_WIDTH] __attribute__((aligned(64)));
                        simd_i32_store(tmp_j, jvec);
                        while (bits) {
                            int lane = __builtin_ctz(bits);
                            bits &= bits - 1;
                            if (n < neighbor->maxneighs) {
                                neighs(neighbor->neighbors,
                                    i,
                                    n,
                                    atom->Nlocal,
                                    neighbor) = tmp_j[lane];
                            }
                            n++;
                        }
                    }
                }
            }

            neighbor->numneigh[i] = n;
            if (n >= neighbor->maxneighs) {
                resize_local = 1;
                if (n > new_maxneighs) new_maxneighs = n;
            }
        }

        resize = resize_local;
        if (resize) {
            printf("RESIZE %d, PROC %d\n", neighbor->maxneighs, me);
            neighbor->maxneighs = new_maxneighs * 1.2;
            free(neighbor->neighbors);
            free(is_inner_buf);
            neighbor->neighbors = (int*)allocate(ALIGNMENT,
                atom->Nmax * neighbor->maxneighs * sizeof(int));
            memset(neighbor->neighbors,
                0,
                atom->Nmax * neighbor->maxneighs * sizeof(int));
            is_inner_buf = (int*)allocate(ALIGNMENT, neighbor->maxneighs * sizeof(int));
        }
    }
}

#endif // !NBLIST_CSR && HAVE_SIMD_NEIGHBOR

#ifdef NBLIST_CSR

/* CSR, scalar bin scan: each thread fills its own growable buffer, then a
 * prefix sum over per-thread totals produces global offsets and buffers are
 * merged into the final compact array. */
void RULE_VARIANT_NAME(buildNeighborLists)(Atom* atom, Neighbor* neighbor)
{
    int Nlocal          = atom->Nlocal;
    int max_threads     = 1;
    int actual_nthreads = 1;
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#endif
    int** t_buf  = (int**)malloc(max_threads * sizeof(int*));
    int* t_cap   = (int*)malloc(max_threads * sizeof(int));
    int* t_pos   = (int*)malloc(max_threads * sizeof(int));
    int* t_owner = (int*)allocate(ALIGNMENT, Nlocal * sizeof(int));
    int init_cap = MAX(1, Nlocal / max_threads) * MAX(1, neighbor->maxneighs);
    for (int t = 0; t < max_threads; t++) {
        t_cap[t] = init_cap;
        t_pos[t] = 0;
        t_buf[t] = (int*)malloc(init_cap * sizeof(int));
    }

#pragma omp parallel
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#pragma omp single
        actual_nthreads = omp_get_num_threads();
#else
        int tid = 0;
#endif

#define CSR_APPEND(jval)                                                                 \
    do {                                                                                 \
        if (t_pos[tid] + n >= t_cap[tid]) {                                              \
            t_cap[tid] *= 2;                                                             \
            t_buf[tid] = (int*)realloc(t_buf[tid], t_cap[tid] * sizeof(int));            \
        }                                                                                \
        t_buf[tid][t_pos[tid] + n] = (jval);                                             \
        n++;                                                                             \
    } while (0)

#pragma omp for schedule(static)
        for (int i = 0; i < Nlocal; i++) {
            t_owner[i]               = tid;
            neighbor->neigh_start[i] = t_pos[tid];

            int n         = 0;
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            int ibin      = coord2bin(xtmp, ytmp, ztmp);
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int type_i = atom->type[i];
#endif
            for (int k = 0; k < nstencil; k++) {
                int jbin     = ibin + stencil[k];
                int* loc_bin = &bins[jbin * atoms_per_bin];
                for (int m = 0; m < bincount[jbin]; m++) {
                    int j = loc_bin[m];
                    if (i == j) continue;
                    if (neighbor->half_neigh && j < i) continue;
                    if (half_stencil && ibin == jbin && skipNeigh(atom, i, j)) continue;
                    MD_FLOAT delx = xtmp - atom_x(j);
                    MD_FLOAT dely = ytmp - atom_y(j);
                    MD_FLOAT delz = ztmp - atom_z(j);
                    MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                    int type_j = atom->type[j];
                    const MD_FLOAT cutoff =
                        atom->cutneighsq[type_i * atom->ntypes + type_j];
#else
                    const MD_FLOAT cutoff = cutneighsq;
#endif
                    if (rsq <= cutoff) {
                        CSR_APPEND(j);
                    }
                }
            }
            neighbor->numneigh[i] = n;
            t_pos[tid] += n;
        }
#undef CSR_APPEND
    }

    int* t_offsets = (int*)allocate(ALIGNMENT, (actual_nthreads + 1) * sizeof(int));
    t_offsets[0]   = 0;
    for (int t = 0; t < actual_nthreads; t++)
        t_offsets[t + 1] = t_offsets[t] + t_pos[t];
    int total = t_offsets[actual_nthreads];

    if (neighbor->neighbors) free(neighbor->neighbors);
    neighbor->neighbors = (int*)allocate(ALIGNMENT,
        (MAX(1, total) + VECTOR_WIDTH) * sizeof(int));
    memset(neighbor->neighbors + total, 0, VECTOR_WIDTH * sizeof(int));

    for (int i = 0; i < Nlocal; i++)
        neighbor->neigh_start[i] += t_offsets[t_owner[i]];
    neighbor->neigh_start[Nlocal] = total;
    for (int t = 0; t < actual_nthreads; t++) {
        memcpy(neighbor->neighbors + t_offsets[t], t_buf[t], t_pos[t] * sizeof(int));
        free(t_buf[t]);
    }

    int new_max = 0;
    for (int i = 0; i < Nlocal; i++)
        if (neighbor->numneigh[i] > new_max) new_max = neighbor->numneigh[i];
    if (new_max != neighbor->maxneighs) {
        neighbor->maxneighs = new_max;
        free(is_inner_buf);
        is_inner_buf = (int*)allocate(ALIGNMENT,
            MAX(1, neighbor->maxneighs) * sizeof(int));
    }

    free(t_owner);
    free(t_buf);
    free(t_cap);
    free(t_pos);
    free(t_offsets);
}

#endif // NBLIST_CSR

#if defined(NBLIST_CSR) && defined(HAVE_SIMD_NEIGHBOR)

/* CSR, vectorized bin scan: same per-thread growable-buffer/prefix-sum
 * structure as the scalar CSR build above, with the same vectorized scan
 * technique as the AOS/SOA SIMD build -- only the storage step differs
 * (CSR_APPEND instead of a neighs() macro write). */
void RULE_VARIANT_NAME(buildNeighborListsSimd)(Atom* atom, Neighbor* neighbor)
{
    int Nlocal          = atom->Nlocal;
    int max_threads     = 1;
    int actual_nthreads = 1;
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#endif
    int** t_buf  = (int**)malloc(max_threads * sizeof(int*));
    int* t_cap   = (int*)malloc(max_threads * sizeof(int));
    int* t_pos   = (int*)malloc(max_threads * sizeof(int));
    int* t_owner = (int*)allocate(ALIGNMENT, Nlocal * sizeof(int));
    int init_cap = MAX(1, Nlocal / max_threads) * MAX(1, neighbor->maxneighs);
    for (int t = 0; t < max_threads; t++) {
        t_cap[t] = init_cap;
        t_pos[t] = 0;
        t_buf[t] = (int*)malloc(init_cap * sizeof(int));
    }

#pragma omp parallel
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#pragma omp single
        actual_nthreads = omp_get_num_threads();
#else
        int tid = 0;
#endif

#define CSR_APPEND(jval)                                                                 \
    do {                                                                                 \
        if (t_pos[tid] + n >= t_cap[tid]) {                                              \
            t_cap[tid] *= 2;                                                             \
            t_buf[tid] = (int*)realloc(t_buf[tid], t_cap[tid] * sizeof(int));            \
        }                                                                                \
        t_buf[tid][t_pos[tid] + n] = (jval);                                             \
        n++;                                                                             \
    } while (0)

#pragma omp for schedule(static)
        for (int i = 0; i < Nlocal; i++) {
            t_owner[i]               = tid;
            neighbor->neigh_start[i] = t_pos[tid];

            int n         = 0;
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            int ibin      = coord2bin(xtmp, ytmp, ztmp);
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int type_i = atom->type[i];
#endif
            MD_SIMD_INT ivec       = simd_i32_broadcast(i);
            MD_SIMD_FLOAT xtmp_vec = simd_real_broadcast(xtmp);
            MD_SIMD_FLOAT ytmp_vec = simd_real_broadcast(ytmp);
            MD_SIMD_FLOAT ztmp_vec = simd_real_broadcast(ztmp);
#if LJ_COMB_RULE == LJ_COMB_SINGLE
            MD_SIMD_FLOAT cutoff_vec_single = simd_real_broadcast(cutneighsq);
#else
            MD_SIMD_INT tbase_i = simd_i32_broadcast(type_i * atom->ntypes);
#endif
            for (int k = 0; k < nstencil; k++) {
                int jbin       = ibin + stencil[k];
                int* loc_bin   = &bins[jbin * atoms_per_bin];
                int bincount_j = bincount[jbin];

                if (half_stencil && ibin == jbin) {
                    for (int m = 0; m < bincount_j; m++) {
                        int j = loc_bin[m];
                        if (i == j) continue;
                        if (neighbor->half_neigh && j < i) continue;
                        if (skipNeigh(atom, i, j)) continue;

                        MD_FLOAT delx = xtmp - atom_x(j);
                        MD_FLOAT dely = ytmp - atom_y(j);
                        MD_FLOAT delz = ztmp - atom_z(j);
                        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                        int type_j = atom->type[j];
                        const MD_FLOAT cutoff =
                            atom->cutneighsq[type_i * atom->ntypes + type_j];
#else
                        const MD_FLOAT cutoff = cutneighsq;
#endif
                        if (rsq <= cutoff) {
                            CSR_APPEND(j);
                        }
                    }
                    continue;
                }

                MD_SIMD_INT bincount_vec = simd_i32_broadcast(bincount_j);

                for (int m = 0; m < bincount_j; m += VECTOR_WIDTH) {
                    MD_SIMD_MASK bound_mask = simd_mask_i32_cond_lt(
                        simd_i32_add(simd_i32_broadcast(m), simd_i32_seq()),
                        bincount_vec);
                    MD_SIMD_INT jvec = simd_i32_mask_load(&loc_bin[m], bound_mask);

#ifdef ATOM_POSITION_AOS
                    MD_SIMD_INT j3   = simd_i32_add(simd_i32_add(jvec, jvec), jvec);
                    MD_SIMD_FLOAT xj = simd_real_gather(j3,
                        &(atom->x[0]),
                        sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT yj = simd_real_gather(j3,
                        &(atom->x[1]),
                        sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT zj = simd_real_gather(j3,
                        &(atom->x[2]),
                        sizeof(MD_FLOAT));
#else
                    MD_SIMD_FLOAT xj = simd_real_gather(jvec, atom->x, sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT yj = simd_real_gather(jvec, atom->y, sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT zj = simd_real_gather(jvec, atom->z, sizeof(MD_FLOAT));
#endif
                    MD_SIMD_FLOAT delx = simd_real_sub(xtmp_vec, xj);
                    MD_SIMD_FLOAT dely = simd_real_sub(ytmp_vec, yj);
                    MD_SIMD_FLOAT delz = simd_real_sub(ztmp_vec, zj);
                    MD_SIMD_FLOAT rsq  = simd_real_fma(delx,
                        delx,
                        simd_real_fma(dely, dely, simd_real_mul(delz, delz)));

                    MD_SIMD_MASK keep = simd_mask_and(bound_mask,
                        simd_mask_not(simd_mask_i32_cond_eq(jvec, ivec)));
                    if (neighbor->half_neigh) {
                        keep = simd_mask_and(keep,
                            simd_mask_not(simd_mask_i32_cond_lt(jvec, ivec)));
                    }

#if LJ_COMB_RULE == LJ_COMB_SINGLE
                    MD_SIMD_FLOAT cutoff_vec = cutoff_vec_single;
#else
                    MD_SIMD_INT type_j_vec   = simd_i32_gather(jvec,
                        atom->type,
                        sizeof(int));
                    MD_SIMD_INT tij          = simd_i32_add(tbase_i, type_j_vec);
                    MD_SIMD_FLOAT cutoff_vec = simd_real_gather(tij,
                        atom->cutneighsq,
                        sizeof(MD_FLOAT));
#endif
                    keep = simd_mask_and(keep, simd_mask_cond_lt(rsq, cutoff_vec));

                    unsigned int bits = simd_mask_to_u32(keep);
                    if (bits) {
                        int tmp_j[VECTOR_WIDTH] __attribute__((aligned(64)));
                        simd_i32_store(tmp_j, jvec);
                        while (bits) {
                            int lane = __builtin_ctz(bits);
                            bits &= bits - 1;
                            CSR_APPEND(tmp_j[lane]);
                        }
                    }
                }
            }
            neighbor->numneigh[i] = n;
            t_pos[tid] += n;
        }
#undef CSR_APPEND
    }

    int* t_offsets = (int*)allocate(ALIGNMENT, (actual_nthreads + 1) * sizeof(int));
    t_offsets[0]   = 0;
    for (int t = 0; t < actual_nthreads; t++)
        t_offsets[t + 1] = t_offsets[t] + t_pos[t];
    int total = t_offsets[actual_nthreads];

    if (neighbor->neighbors) free(neighbor->neighbors);
    neighbor->neighbors = (int*)allocate(ALIGNMENT,
        (MAX(1, total) + VECTOR_WIDTH) * sizeof(int));
    memset(neighbor->neighbors + total, 0, VECTOR_WIDTH * sizeof(int));

    for (int i = 0; i < Nlocal; i++)
        neighbor->neigh_start[i] += t_offsets[t_owner[i]];
    neighbor->neigh_start[Nlocal] = total;
    for (int t = 0; t < actual_nthreads; t++) {
        memcpy(neighbor->neighbors + t_offsets[t], t_buf[t], t_pos[t] * sizeof(int));
        free(t_buf[t]);
    }

    int new_max = 0;
    for (int i = 0; i < Nlocal; i++)
        if (neighbor->numneigh[i] > new_max) new_max = neighbor->numneigh[i];
    if (new_max != neighbor->maxneighs) {
        neighbor->maxneighs = new_max;
        free(is_inner_buf);
        is_inner_buf = (int*)allocate(ALIGNMENT,
            MAX(1, neighbor->maxneighs) * sizeof(int));
    }

    free(t_owner);
    free(t_buf);
    free(t_cap);
    free(t_pos);
    free(t_offsets);
}

#endif

void RULE_VARIANT_NAME(neighborGhost)(Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("neighborGhost begin\n");

    int Nshell = 0;
    int Nlocal = atom->Nlocal;
    int Nghost = atom->Nghost;
    if (neighbor->listshell) free(neighbor->listshell);
    neighbor->listshell = (int*)allocate(ALIGNMENT, Nghost * sizeof(int));
    int* listzone       = (int*)allocate(ALIGNMENT, 8 * Nghost * sizeof(int));
    int countAtoms[8]   = { 0, 0, 0, 0, 0, 0, 0, 0 };

    // Selecting ghost atoms for interaction
    for (int i = Nlocal; i < Nlocal + Nghost; i++) {
        int izone = ghostZone(atom, i);
        int* list = &listzone[Nghost * izone];
        int n     = countAtoms[izone];
        list[n]   = i;
        countAtoms[izone]++;
    }

    for (int zone = 1; zone <= 3; zone++) {
        int* list = &listzone[Nghost * zone];
        for (int n = 0; n < countAtoms[zone]; n++)
            neighbor->listshell[Nshell++] = list[n];
    }

    neighbor->Nshell = Nshell;
    if (neighbor->numNeighShell) free(neighbor->numNeighShell);
    if (neighbor->neighshell) free(neighbor->neighshell);
    neighbor->neighshell    = (int*)allocate(ALIGNMENT,
        Nshell * neighbor->maxneighs * sizeof(int));
    neighbor->numNeighShell = (int*)allocate(ALIGNMENT, Nshell * sizeof(int));
    int resize              = 1;

    while (resize) {
        resize = 0;
        for (int i = 0; i < Nshell; i++) {
            int n         = 0;
            int iatom     = neighbor->listshell[i];
            int izone     = ghostZone(atom, iatom);
            MD_FLOAT xtmp = atom_x(iatom);
            MD_FLOAT ytmp = atom_y(iatom);
            MD_FLOAT ztmp = atom_z(iatom);
            int ibin      = coord2bin(xtmp, ytmp, ztmp);

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int type_i = atom->type[iatom];
#endif

            for (int k = 0; k < nstencil; k++) {
                int jbin     = ibin + stencil[k];
                int* loc_bin = &bins[jbin * atoms_per_bin];
                for (int m = 0; m < bincount[jbin]; m++) {
                    int jatom = loc_bin[m];

                    int jzone = ghostZone(atom, jatom);

                    if (jzone <= izone) continue;
                    if (izone == 1 && (jzone == 5 || jzone == 6 || jzone == 7)) continue;
                    if (izone == 2 && (jzone == 4 || jzone == 6 || jzone == 7)) continue;
                    if (izone == 3 && (jzone == 4 || jzone == 5 || jzone == 7)) continue;

                    MD_FLOAT delx = xtmp - atom_x(jatom);
                    MD_FLOAT dely = ytmp - atom_y(jatom);
                    MD_FLOAT delz = ztmp - atom_z(jatom);
                    MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                    int type_j = atom->type[jatom];
                    const MD_FLOAT cutoff =
                        atom->cutneighsq[type_i * atom->ntypes + type_j];
#else
                    const MD_FLOAT cutoff = cutneighsq;
#endif
                    if (rsq <= cutoff) {
                        neighshell(neighbor->neighshell,
                            i,
                            n,
                            neighbor->maxneighs) = jatom;
                        n++;
                    }
                }
            }

            neighbor->numNeighShell[i] = n;
            if (n >= neighbor->maxneighs) {
                resize              = 1;
                neighbor->maxneighs = n * 1.2;
                break;
            }
        }

        if (resize) {
            free(neighbor->neighshell);
            neighbor->neighshell = (int*)allocate(ALIGNMENT,
                Nshell * neighbor->maxneighs * sizeof(int));
        }
    }

    free(listzone);
    DEBUG_MESSAGE("neighborGhost end\n");
}
