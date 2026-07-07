/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 *
 * SIMD-optimized Lennard-Jones force kernels for Verlet Lists
 * Supports: AVX2, AVX512, NEON, SVE (double precision)
 * Requires: __SIMD_KERNEL__ flag and NBLIST_AOS layout
 *
 * LJ combination rules (compile-time via -DLJ_COMB_RULE=<value>):
 *   LJ_COMB_SINGLE (0): Single atom type - broadcast global epsilon/sigma
 *   LJ_COMB_GEOM   (1): Geometric - sqrt(eps_i*eps_j), sigma3_i*sigma3_j
 *   LJ_COMB_NONE   (2): Full type-pair matrix lookup via gather
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//---
#include <atom.h>
#include <likwid-marker.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>

#ifdef __SIMD_KERNEL__
#include <simd.h>
#endif

#if defined(NBLIST_PAIRLIST) && defined(__SIMD_KERNEL__)
#include <allocate.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

// Compile-time guards for unsupported configurations.
// Only error out when the SIMD kernel is actually selected for this build;
// on CUDA/HIP builds, SIMD is forced to NONE and these functions are never
// called (force.c dispatches to computeForceLJCUDA), so the file just needs
// to compile to stubs.
#if defined(NBLIST_SOA) && defined(__SIMD_KERNEL__)
#error "SIMD kernel not implemented when NBLIST_DATA_LAYOUT is SOA"
#endif
#if defined(NBLIST_CSR) && defined(__SIMD_KERNEL__)
#error "SIMD kernel not implemented when NBLIST_DATA_LAYOUT is CSR"
#endif

double computeForceLJFullNeigh_simd(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    int Nlocal = atom->Nlocal;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#endif

    for (int i = 0; i < Nlocal; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    double S = getTimeStamp();

#ifndef __SIMD_KERNEL__
    fprintf(stderr, "Error: SIMD kernel not implemented for specified instruction set!");
    exit(-1);
#else
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
#else
    // Cutoff is uniform for all types, broadcast it
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(param->cutforce * param->cutforce);
#endif
    MD_SIMD_FLOAT c48_vec        = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec        = simd_real_broadcast(0.5);

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < Nlocal; i++) {
#ifdef NBLIST_CSR
            int* neighs               = &neighbor->neighbors[neighbor->neigh_start[i]];
#else
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
            MD_SIMD_FLOAT xtmp        = simd_real_broadcast(atom_x(i));
            MD_SIMD_FLOAT ytmp        = simd_real_broadcast(atom_y(i));
            MD_SIMD_FLOAT ztmp        = simd_real_broadcast(atom_z(i));
            MD_SIMD_FLOAT fix         = simd_real_zero();
            MD_SIMD_FLOAT fiy         = simd_real_zero();
            MD_SIMD_FLOAT fiz         = simd_real_zero();

#if LJ_COMB_RULE == LJ_COMB_GEOM
            // Broadcast per-atom LJ params for atom i (geometric combination)
            MD_SIMD_FLOAT sqrt_eps_i = simd_real_broadcast(atom->sqrt_epsilon[i]);
            MD_SIMD_FLOAT sigma3_i   = simd_real_broadcast(atom->sigma3[i]);
#elif LJ_COMB_RULE == LJ_COMB_NONE
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs; k += VECTOR_WIDTH) {
                // If the last iteration of this loop is separated from the rest, this
                // mask can be set only there
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);

#if LJ_COMB_RULE == LJ_COMB_GEOM
                // Direct gather of per-atom LJ params (avoids type index computation)
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j,
                    atom->sqrt_epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j,
                    atom->sigma3,
                    sizeof(MD_FLOAT));
                // Geometric combination: eps_ij = sqrt(eps_i) * sqrt(eps_j), sigma6_ij =
                // sigma3_i * sigma3_j
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_NONE
                MD_SIMD_INT tj           = simd_i32_gather(j, atom->type, sizeof(int));
                MD_SIMD_INT tij          = simd_i32_add(tbase_i, tj);
                MD_SIMD_FLOAT sigma6_vec = simd_real_gather(tij,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_gather(tij,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#endif

#ifdef ATOM_POSITION_AOS
                MD_SIMD_INT j3           = simd_i32_add(simd_i32_add(j, j), j); // j * 3
                MD_SIMD_FLOAT delx       = xtmp - simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely       = ytmp - simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz       = ztmp - simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT delx       = xtmp -
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely = ytmp -
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz = ztmp -
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT));
#endif
                MD_SIMD_FLOAT rsq        = simd_real_fma(delx,
                    delx,
                    simd_real_fma(dely, dely, simd_real_mul(delz, delz)));
                MD_SIMD_MASK cutoff_mask = simd_mask_and(mask_numneighs,
                    simd_mask_cond_lt(rsq, cutforcesq_vec));
                MD_SIMD_FLOAT sr2        = simd_real_reciprocal(rsq);
                MD_SIMD_FLOAT sr6        = simd_real_mul(sr2,
                    simd_real_mul(sr2, simd_real_mul(sr2, sigma6_vec)));
                MD_SIMD_FLOAT force      = simd_real_mul(c48_vec,
                    simd_real_mul(sr6,
                        simd_real_mul(simd_real_sub(sr6, c05_vec),
                            simd_real_mul(sr2, eps_vec))));

                fix = simd_real_masked_add(fix, simd_real_mul(delx, force), cutoff_mask);
                fiy = simd_real_masked_add(fiy, simd_real_mul(dely, force), cutoff_mask);
                fiz = simd_real_masked_add(fiz, simd_real_mul(delz, force), cutoff_mask);
            }

            atom_fx(i) += simd_real_h_reduce_sum(fix);
            atom_fy(i) += simd_real_h_reduce_sum(fiy);
            atom_fz(i) += simd_real_h_reduce_sum(fiz);
        }

        LIKWID_MARKER_STOP("force");
    }
#endif

    double E = getTimeStamp();
    return E - S;
}

double computeForceLJHalfNeigh_simd(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    int Nlocal = atom->Nlocal;
    int Nghost = atom->Nghost;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#endif

    for (int i = 0; i < Nlocal + Nghost; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    double S = getTimeStamp();

#ifndef __SIMD_KERNEL__
    fprintf(stderr, "Error: SIMD kernel not implemented for specified instruction set!");
    exit(-1);
#else
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
#else
    // Cutoff is uniform for all types, broadcast it
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(param->cutforce * param->cutforce);
#endif
    MD_SIMD_FLOAT c48_vec        = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec        = simd_real_broadcast(0.5);
    MD_SIMD_INT nlocal_vec       = simd_i32_broadcast(Nlocal);

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < Nlocal; i++) {
#ifdef NBLIST_CSR
            int* neighs               = &neighbor->neighbors[neighbor->neigh_start[i]];
#else
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
            MD_SIMD_FLOAT xtmp        = simd_real_broadcast(atom_x(i));
            MD_SIMD_FLOAT ytmp        = simd_real_broadcast(atom_y(i));
            MD_SIMD_FLOAT ztmp        = simd_real_broadcast(atom_z(i));
            MD_SIMD_FLOAT fix         = simd_real_zero();
            MD_SIMD_FLOAT fiy         = simd_real_zero();
            MD_SIMD_FLOAT fiz         = simd_real_zero();

#if LJ_COMB_RULE == LJ_COMB_GEOM
            // Broadcast per-atom LJ params for atom i (geometric combination)
            MD_SIMD_FLOAT sqrt_eps_i = simd_real_broadcast(atom->sqrt_epsilon[i]);
            MD_SIMD_FLOAT sigma3_i   = simd_real_broadcast(atom->sigma3[i]);
#elif LJ_COMB_RULE == LJ_COMB_NONE
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs; k += VECTOR_WIDTH) {
                // Mask for valid neighbors in this iteration
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);

#if LJ_COMB_RULE == LJ_COMB_GEOM
                // Direct gather of per-atom LJ params (avoids type index computation)
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j,
                    atom->sqrt_epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j,
                    atom->sigma3,
                    sizeof(MD_FLOAT));
                // Geometric combination: eps_ij = sqrt(eps_i) * sqrt(eps_j), sigma6_ij =
                // sigma3_i * sigma3_j
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_NONE
                MD_SIMD_INT tj           = simd_i32_gather(j, atom->type, sizeof(int));
                MD_SIMD_INT tij          = simd_i32_add(tbase_i, tj);
                MD_SIMD_FLOAT sigma6_vec = simd_real_gather(tij,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_gather(tij,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#endif

#ifdef ATOM_POSITION_AOS
                MD_SIMD_INT j3           = simd_i32_add(simd_i32_add(j, j), j); // j * 3
                MD_SIMD_FLOAT delx       = xtmp - simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely       = ytmp - simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz       = ztmp - simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT delx       = xtmp -
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely = ytmp -
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz = ztmp -
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT));
#endif
                MD_SIMD_FLOAT rsq        = simd_real_fma(delx,
                    delx,
                    simd_real_fma(dely, dely, simd_real_mul(delz, delz)));
                MD_SIMD_MASK cutoff_mask = simd_mask_and(mask_numneighs,
                    simd_mask_cond_lt(rsq, cutforcesq_vec));
                MD_SIMD_FLOAT sr2        = simd_real_reciprocal(rsq);
                MD_SIMD_FLOAT sr6        = simd_real_mul(sr2,
                    simd_real_mul(sr2, simd_real_mul(sr2, sigma6_vec)));
                MD_SIMD_FLOAT force      = simd_real_mul(c48_vec,
                    simd_real_mul(sr6,
                        simd_real_mul(simd_real_sub(sr6, c05_vec),
                            simd_real_mul(sr2, eps_vec))));

                // Compute force components and accumulate for atom i
                MD_SIMD_FLOAT fx_tmp = simd_real_mul(delx, force);
                MD_SIMD_FLOAT fy_tmp = simd_real_mul(dely, force);
                MD_SIMD_FLOAT fz_tmp = simd_real_mul(delz, force);

                fix = simd_real_masked_add(fix, fx_tmp, cutoff_mask);
                fiy = simd_real_masked_add(fiy, fy_tmp, cutoff_mask);
                fiz = simd_real_masked_add(fiz, fz_tmp, cutoff_mask);

                // Apply Newton's third law using vectorized scatter
                // Note: not thread-safe under OpenMP (no atomic scatter support)
                MD_SIMD_MASK j_local_mask  = simd_mask_i32_cond_lt(j, nlocal_vec);
                MD_SIMD_MASK j_update_mask = simd_mask_and(cutoff_mask,
                    param->method ? cutoff_mask : j_local_mask);
#ifdef ATOM_POSITION_AOS
                simd_real_masked_scatter_sub(&atom->fx[0], j3, fx_tmp, j_update_mask);
                simd_real_masked_scatter_sub(&atom->fx[1], j3, fy_tmp, j_update_mask);
                simd_real_masked_scatter_sub(&atom->fx[2], j3, fz_tmp, j_update_mask);
#else
                simd_real_masked_scatter_sub(atom->fx, j, fx_tmp, j_update_mask);
                simd_real_masked_scatter_sub(atom->fy, j, fy_tmp, j_update_mask);
                simd_real_masked_scatter_sub(atom->fz, j, fz_tmp, j_update_mask);
#endif
            }

            atom_fx(i) += simd_real_h_reduce_sum(fix);
            atom_fy(i) += simd_real_h_reduce_sum(fiy);
            atom_fz(i) += simd_real_h_reduce_sum(fiz);
        }

        LIKWID_MARKER_STOP("force");
    }
#endif

    double E = getTimeStamp();
    return E - S;
}

#ifdef __SIMD_COMPRESS__

// Extra per-i (and, for LJ_COMB_NONE, per-pair) state that ljforce_compressed() needs
// on top of the buffered (j, delx, dely, delz) triples, depending on LJ_COMB_RULE.
#if LJ_COMB_RULE == LJ_COMB_GEOM
#define LJ_EXTRA_PARAMS MD_SIMD_FLOAT sqrt_eps_i, MD_SIMD_FLOAT sigma3_i
#define LJ_EXTRA_ARGS   sqrt_eps_i, sigma3_i
#elif LJ_COMB_RULE == LJ_COMB_NONE
#define LJ_EXTRA_PARAMS MD_SIMD_INT tbase_i
#define LJ_EXTRA_ARGS   tbase_i
#else
#define LJ_EXTRA_PARAMS MD_SIMD_FLOAT eps_vec, MD_SIMD_FLOAT sigma6_vec
#define LJ_EXTRA_ARGS   eps_vec, sigma6_vec
#endif

// Computes the LJ force magnitude/r factor for a fully-compacted vector of pairs
// that have already passed the cutoff check, so (unlike the masked kernels above)
// no per-lane masking is needed here.
static inline MD_SIMD_FLOAT ljforce_compressed(Atom* atom, MD_SIMD_INT jv,
    MD_SIMD_FLOAT dvx, MD_SIMD_FLOAT dvy, MD_SIMD_FLOAT dvz, MD_SIMD_FLOAT c48_vec,
    MD_SIMD_FLOAT c05_vec, LJ_EXTRA_PARAMS)
{
    MD_SIMD_FLOAT rsq = simd_real_fma(
        dvx, dvx, simd_real_fma(dvy, dvy, simd_real_mul(dvz, dvz)));

#if LJ_COMB_RULE == LJ_COMB_GEOM
    MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(jv, atom->sqrt_epsilon, sizeof(MD_FLOAT));
    MD_SIMD_FLOAT sigma3_j   = simd_real_gather(jv, atom->sigma3, sizeof(MD_FLOAT));
    MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
    MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_NONE
    MD_SIMD_INT tj           = simd_i32_gather(jv, atom->type, sizeof(int));
    MD_SIMD_INT tij          = simd_i32_add(tbase_i, tj);
    MD_SIMD_FLOAT sigma6_vec = simd_real_gather(tij, atom->sigma6, sizeof(MD_FLOAT));
    MD_SIMD_FLOAT eps_vec    = simd_real_gather(tij, atom->epsilon, sizeof(MD_FLOAT));
#endif

    MD_SIMD_FLOAT sr2 = simd_real_reciprocal(rsq);
    MD_SIMD_FLOAT sr6 = simd_real_mul(
        sr2, simd_real_mul(sr2, simd_real_mul(sr2, sigma6_vec)));
    return simd_real_mul(c48_vec,
        simd_real_mul(sr6,
            simd_real_mul(simd_real_sub(sr6, c05_vec), simd_real_mul(sr2, eps_vec))));
}

// SIMD lane-compression variants of the kernels above: instead of running the
// expensive LJ stage (parameter gather, reciprocal, force calc) on every raw
// VECTOR_WIDTH-wide block of the neighbor list and masking out lanes that fail
// the cutoff check, only the (index, delx, dely, delz) of pairs that already
// passed the cutoff are appended to a small staging buffer. The expensive stage
// only runs once that buffer holds a full, unmasked vector's worth of valid
// pairs, or once as a masked "epilogue" for whatever is left over after atom
// i's neighbor list is exhausted. See TODO.md and the design plan for details.
double computeForceLJFullNeigh_simd_compress(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    int Nlocal = atom->Nlocal;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#endif

    for (int i = 0; i < Nlocal; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    double S = getTimeStamp();

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
#else
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(param->cutforce * param->cutforce);
#endif
    MD_SIMD_FLOAT c48_vec = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec = simd_real_broadcast(0.5);
    MD_SIMD_MASK all_true = simd_mask_from_u32((1u << VECTOR_WIDTH) - 1u);

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

        // Per-thread staging buffers: up to 2*VECTOR_WIDTH-1 cutoff-passing pairs can
        // be buffered at once (one raw block's worth on top of a not-yet-flushed
        // full vector), then compacted back down below VECTOR_WIDTH after a flush.
        int buf_j[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT buf_delx[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT buf_dely[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT buf_delz[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        int tmp_j[VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT tmp_delx[VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT tmp_dely[VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT tmp_delz[VECTOR_WIDTH] __attribute__((aligned(64)));

#pragma omp for schedule(runtime)
        for (int i = 0; i < Nlocal; i++) {
#ifdef NBLIST_CSR
            int* neighs               = &neighbor->neighbors[neighbor->neigh_start[i]];
#else
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
            MD_SIMD_FLOAT xtmp        = simd_real_broadcast(atom_x(i));
            MD_SIMD_FLOAT ytmp        = simd_real_broadcast(atom_y(i));
            MD_SIMD_FLOAT ztmp        = simd_real_broadcast(atom_z(i));
            MD_SIMD_FLOAT fix         = simd_real_zero();
            MD_SIMD_FLOAT fiy         = simd_real_zero();
            MD_SIMD_FLOAT fiz         = simd_real_zero();
            int cnt                   = 0;

            // Zero buf_j so that any not-yet-written slots read by a partial
            // epilogue flush are safe gather indices (0), not stack garbage.
            memset(buf_j, 0, sizeof(buf_j));

#if LJ_COMB_RULE == LJ_COMB_GEOM
            MD_SIMD_FLOAT sqrt_eps_i = simd_real_broadcast(atom->sqrt_epsilon[i]);
            MD_SIMD_FLOAT sigma3_i   = simd_real_broadcast(atom->sigma3[i]);
#elif LJ_COMB_RULE == LJ_COMB_NONE
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs; k += VECTOR_WIDTH) {
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);

#ifdef ATOM_POSITION_AOS
                MD_SIMD_INT j3           = simd_i32_add(simd_i32_add(j, j), j);
                MD_SIMD_FLOAT delx       = xtmp - simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely       = ytmp - simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz       = ztmp - simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT delx       = xtmp -
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely = ytmp -
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz = ztmp -
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT));
#endif
                MD_SIMD_FLOAT rsq        = simd_real_fma(delx,
                    delx,
                    simd_real_fma(dely, dely, simd_real_mul(delz, delz)));
                MD_SIMD_MASK cutoff_mask = simd_mask_and(mask_numneighs,
                    simd_mask_cond_lt(rsq, cutforcesq_vec));

                unsigned int bits = simd_mask_to_u32(cutoff_mask);
                if (bits == 0) {
                    continue;
                }

                simd_i32_store(tmp_j, j);
                simd_real_store(tmp_delx, delx);
                simd_real_store(tmp_dely, dely);
                simd_real_store(tmp_delz, delz);

                while (bits) {
                    int lane = __builtin_ctz(bits);
                    bits &= bits - 1;
                    buf_j[cnt]    = tmp_j[lane];
                    buf_delx[cnt] = tmp_delx[lane];
                    buf_dely[cnt] = tmp_dely[lane];
                    buf_delz[cnt] = tmp_delz[lane];
                    cnt++;
                }

                while (cnt >= VECTOR_WIDTH) {
                    MD_SIMD_INT jv    = simd_i32_load(buf_j);
                    MD_SIMD_FLOAT dvx = simd_real_load(buf_delx);
                    MD_SIMD_FLOAT dvy = simd_real_load(buf_dely);
                    MD_SIMD_FLOAT dvz = simd_real_load(buf_delz);
                    MD_SIMD_FLOAT force = ljforce_compressed(
                        atom, jv, dvx, dvy, dvz, c48_vec, c05_vec, LJ_EXTRA_ARGS);

                    fix = simd_real_masked_add(fix, simd_real_mul(dvx, force), all_true);
                    fiy = simd_real_masked_add(fiy, simd_real_mul(dvy, force), all_true);
                    fiz = simd_real_masked_add(fiz, simd_real_mul(dvz, force), all_true);

                    int remaining = cnt - VECTOR_WIDTH;
                    if (remaining > 0) {
                        memmove(buf_j, buf_j + VECTOR_WIDTH, remaining * sizeof(int));
                        memmove(buf_delx,
                            buf_delx + VECTOR_WIDTH,
                            remaining * sizeof(MD_FLOAT));
                        memmove(buf_dely,
                            buf_dely + VECTOR_WIDTH,
                            remaining * sizeof(MD_FLOAT));
                        memmove(buf_delz,
                            buf_delz + VECTOR_WIDTH,
                            remaining * sizeof(MD_FLOAT));
                    }
                    cnt = remaining;
                }
            }

            // Epilogue: whatever is left in the buffer once atom i's neighbor list
            // is exhausted is processed as one final masked (partial) vector.
            if (cnt > 0) {
                MD_SIMD_MASK epilogue_mask = simd_mask_from_u32((1u << cnt) - 1u);
                MD_SIMD_INT jv    = simd_i32_load(buf_j);
                MD_SIMD_FLOAT dvx = simd_real_load(buf_delx);
                MD_SIMD_FLOAT dvy = simd_real_load(buf_dely);
                MD_SIMD_FLOAT dvz = simd_real_load(buf_delz);
                MD_SIMD_FLOAT force = ljforce_compressed(
                    atom, jv, dvx, dvy, dvz, c48_vec, c05_vec, LJ_EXTRA_ARGS);

                fix = simd_real_masked_add(fix, simd_real_mul(dvx, force), epilogue_mask);
                fiy = simd_real_masked_add(fiy, simd_real_mul(dvy, force), epilogue_mask);
                fiz = simd_real_masked_add(fiz, simd_real_mul(dvz, force), epilogue_mask);
            }

            atom_fx(i) += simd_real_h_reduce_sum(fix);
            atom_fy(i) += simd_real_h_reduce_sum(fiy);
            atom_fz(i) += simd_real_h_reduce_sum(fiz);
        }

        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    return E - S;
}

double computeForceLJHalfNeigh_simd_compress(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    int Nlocal = atom->Nlocal;
    int Nghost = atom->Nghost;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#endif

    for (int i = 0; i < Nlocal + Nghost; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    double S = getTimeStamp();

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
#else
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(param->cutforce * param->cutforce);
#endif
    MD_SIMD_FLOAT c48_vec  = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec  = simd_real_broadcast(0.5);
    MD_SIMD_INT nlocal_vec = simd_i32_broadcast(Nlocal);
    MD_SIMD_MASK all_true  = simd_mask_from_u32((1u << VECTOR_WIDTH) - 1u);

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

        int buf_j[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT buf_delx[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT buf_dely[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT buf_delz[2 * VECTOR_WIDTH] __attribute__((aligned(64)));
        int tmp_j[VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT tmp_delx[VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT tmp_dely[VECTOR_WIDTH] __attribute__((aligned(64)));
        MD_FLOAT tmp_delz[VECTOR_WIDTH] __attribute__((aligned(64)));

#pragma omp for schedule(runtime)
        for (int i = 0; i < Nlocal; i++) {
#ifdef NBLIST_CSR
            int* neighs               = &neighbor->neighbors[neighbor->neigh_start[i]];
#else
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
            MD_SIMD_FLOAT xtmp        = simd_real_broadcast(atom_x(i));
            MD_SIMD_FLOAT ytmp        = simd_real_broadcast(atom_y(i));
            MD_SIMD_FLOAT ztmp        = simd_real_broadcast(atom_z(i));
            MD_SIMD_FLOAT fix         = simd_real_zero();
            MD_SIMD_FLOAT fiy         = simd_real_zero();
            MD_SIMD_FLOAT fiz         = simd_real_zero();
            int cnt                   = 0;

            memset(buf_j, 0, sizeof(buf_j));

#if LJ_COMB_RULE == LJ_COMB_GEOM
            MD_SIMD_FLOAT sqrt_eps_i = simd_real_broadcast(atom->sqrt_epsilon[i]);
            MD_SIMD_FLOAT sigma3_i   = simd_real_broadcast(atom->sigma3[i]);
#elif LJ_COMB_RULE == LJ_COMB_NONE
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs; k += VECTOR_WIDTH) {
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);

#ifdef ATOM_POSITION_AOS
                MD_SIMD_INT j3           = simd_i32_add(simd_i32_add(j, j), j);
                MD_SIMD_FLOAT delx       = xtmp - simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely       = ytmp - simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz       = ztmp - simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT delx       = xtmp -
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT dely = ytmp -
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT delz = ztmp -
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT));
#endif
                MD_SIMD_FLOAT rsq        = simd_real_fma(delx,
                    delx,
                    simd_real_fma(dely, dely, simd_real_mul(delz, delz)));
                MD_SIMD_MASK cutoff_mask = simd_mask_and(mask_numneighs,
                    simd_mask_cond_lt(rsq, cutforcesq_vec));

                unsigned int bits = simd_mask_to_u32(cutoff_mask);
                if (bits == 0) {
                    continue;
                }

                simd_i32_store(tmp_j, j);
                simd_real_store(tmp_delx, delx);
                simd_real_store(tmp_dely, dely);
                simd_real_store(tmp_delz, delz);

                while (bits) {
                    int lane = __builtin_ctz(bits);
                    bits &= bits - 1;
                    buf_j[cnt]    = tmp_j[lane];
                    buf_delx[cnt] = tmp_delx[lane];
                    buf_dely[cnt] = tmp_dely[lane];
                    buf_delz[cnt] = tmp_delz[lane];
                    cnt++;
                }

                while (cnt >= VECTOR_WIDTH) {
                    MD_SIMD_INT jv    = simd_i32_load(buf_j);
                    MD_SIMD_FLOAT dvx = simd_real_load(buf_delx);
                    MD_SIMD_FLOAT dvy = simd_real_load(buf_dely);
                    MD_SIMD_FLOAT dvz = simd_real_load(buf_delz);
                    MD_SIMD_FLOAT force = ljforce_compressed(
                        atom, jv, dvx, dvy, dvz, c48_vec, c05_vec, LJ_EXTRA_ARGS);

                    MD_SIMD_FLOAT fx_tmp = simd_real_mul(dvx, force);
                    MD_SIMD_FLOAT fy_tmp = simd_real_mul(dvy, force);
                    MD_SIMD_FLOAT fz_tmp = simd_real_mul(dvz, force);

                    fix = simd_real_masked_add(fix, fx_tmp, all_true);
                    fiy = simd_real_masked_add(fiy, fy_tmp, all_true);
                    fiz = simd_real_masked_add(fiz, fz_tmp, all_true);

                    // All lanes already passed the cutoff, so the update mask only
                    // needs to exclude ghost atoms (unless param->method wants them too).
                    MD_SIMD_MASK j_local_mask  = simd_mask_i32_cond_lt(jv, nlocal_vec);
                    MD_SIMD_MASK j_update_mask = simd_mask_and(all_true,
                        param->method ? all_true : j_local_mask);
#ifdef ATOM_POSITION_AOS
                    MD_SIMD_INT jv3 = simd_i32_add(simd_i32_add(jv, jv), jv);
                    simd_real_masked_scatter_sub(&atom->fx[0], jv3, fx_tmp, j_update_mask);
                    simd_real_masked_scatter_sub(&atom->fx[1], jv3, fy_tmp, j_update_mask);
                    simd_real_masked_scatter_sub(&atom->fx[2], jv3, fz_tmp, j_update_mask);
#else
                    simd_real_masked_scatter_sub(atom->fx, jv, fx_tmp, j_update_mask);
                    simd_real_masked_scatter_sub(atom->fy, jv, fy_tmp, j_update_mask);
                    simd_real_masked_scatter_sub(atom->fz, jv, fz_tmp, j_update_mask);
#endif

                    int remaining = cnt - VECTOR_WIDTH;
                    if (remaining > 0) {
                        memmove(buf_j, buf_j + VECTOR_WIDTH, remaining * sizeof(int));
                        memmove(buf_delx,
                            buf_delx + VECTOR_WIDTH,
                            remaining * sizeof(MD_FLOAT));
                        memmove(buf_dely,
                            buf_dely + VECTOR_WIDTH,
                            remaining * sizeof(MD_FLOAT));
                        memmove(buf_delz,
                            buf_delz + VECTOR_WIDTH,
                            remaining * sizeof(MD_FLOAT));
                    }
                    cnt = remaining;
                }
            }

            if (cnt > 0) {
                MD_SIMD_MASK epilogue_mask = simd_mask_from_u32((1u << cnt) - 1u);
                MD_SIMD_INT jv    = simd_i32_load(buf_j);
                MD_SIMD_FLOAT dvx = simd_real_load(buf_delx);
                MD_SIMD_FLOAT dvy = simd_real_load(buf_dely);
                MD_SIMD_FLOAT dvz = simd_real_load(buf_delz);
                MD_SIMD_FLOAT force = ljforce_compressed(
                    atom, jv, dvx, dvy, dvz, c48_vec, c05_vec, LJ_EXTRA_ARGS);

                MD_SIMD_FLOAT fx_tmp = simd_real_mul(dvx, force);
                MD_SIMD_FLOAT fy_tmp = simd_real_mul(dvy, force);
                MD_SIMD_FLOAT fz_tmp = simd_real_mul(dvz, force);

                fix = simd_real_masked_add(fix, fx_tmp, epilogue_mask);
                fiy = simd_real_masked_add(fiy, fy_tmp, epilogue_mask);
                fiz = simd_real_masked_add(fiz, fz_tmp, epilogue_mask);

                MD_SIMD_MASK j_local_mask  = simd_mask_i32_cond_lt(jv, nlocal_vec);
                MD_SIMD_MASK j_update_mask = simd_mask_and(epilogue_mask,
                    param->method ? epilogue_mask : j_local_mask);
#ifdef ATOM_POSITION_AOS
                MD_SIMD_INT jv3 = simd_i32_add(simd_i32_add(jv, jv), jv);
                simd_real_masked_scatter_sub(&atom->fx[0], jv3, fx_tmp, j_update_mask);
                simd_real_masked_scatter_sub(&atom->fx[1], jv3, fy_tmp, j_update_mask);
                simd_real_masked_scatter_sub(&atom->fx[2], jv3, fz_tmp, j_update_mask);
#else
                simd_real_masked_scatter_sub(atom->fx, jv, fx_tmp, j_update_mask);
                simd_real_masked_scatter_sub(atom->fy, jv, fy_tmp, j_update_mask);
                simd_real_masked_scatter_sub(atom->fz, jv, fz_tmp, j_update_mask);
#endif
            }

            atom_fx(i) += simd_real_h_reduce_sum(fix);
            atom_fy(i) += simd_real_h_reduce_sum(fiy);
            atom_fz(i) += simd_real_h_reduce_sum(fiz);
        }

        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    return E - S;
}

#undef LJ_EXTRA_PARAMS
#undef LJ_EXTRA_ARGS

#endif // __SIMD_COMPRESS__

#if defined(NBLIST_PAIRLIST) && defined(__SIMD_KERNEL__)
// Per-thread scratch force buffer, same shape/purpose as the one in
// force_lj.c's computeForceLJPairListRef (duplicated rather than shared,
// matching this file's existing self-contained-translation-unit style).
static MD_FLOAT* pairlist_simd_fbuf   = NULL;
static size_t pairlist_simd_fbuf_cap  = 0;

static MD_FLOAT* getPairlistSimdForceBuffer(int nall, int nthreads)
{
    size_t needed = (size_t)nthreads * (size_t)nall * 3;
    if (needed > pairlist_simd_fbuf_cap) {
        if (pairlist_simd_fbuf) free(pairlist_simd_fbuf);
        size_t alloc_size    = needed > 0 ? needed : 1;
        pairlist_simd_fbuf     = (MD_FLOAT*)allocate(ALIGNMENT, alloc_size * sizeof(MD_FLOAT));
        pairlist_simd_fbuf_cap = needed;
    }
    return pairlist_simd_fbuf;
}

// SIMD pair-list kernel: unlike the per-i kernels above, which broadcast one
// atom i's position/params across all lanes and only gather j, this gathers
// BOTH i and j per lane, so a single vector instruction mixes interactions
// belonging to different i atoms (flattenPairList() interleaves the pair
// list round-robin across i's specifically to make this land in practice).
//
// The expensive math (distance, reciprocal, LJ formula) stays fully
// vectorized. Only the force write-back is scalarized: two lanes in the
// same chunk could target the same atom (as either i or j), so a vector
// scatter-store would race within the instruction; unpacking to scalars and
// adding sequentially is always correct regardless of duplicate indices,
// at the same asymptotic cost the scalar reference kernel already pays per
// pair for its own write-back.
double computeForceLJPairList_simd(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJPairList_simd begin\n");

    int nlocal = atom->Nlocal;
    int nghost = atom->Nghost;
    int nall   = nlocal + nghost;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#endif

    for (int i = 0; i < nall; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    MD_FLOAT* fbuf = getPairlistSimdForceBuffer(MAX(1, nall), nthreads);
    int npairs     = neighbor->npairs;

    double S = getTimeStamp();

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
#else
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(param->cutforce * param->cutforce);
#endif
#if LJ_COMB_RULE == LJ_COMB_NONE
    MD_SIMD_INT ntypes_vec = simd_i32_broadcast(atom->ntypes);
#endif
    MD_SIMD_FLOAT c48_vec  = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec  = simd_real_broadcast(0.5);
    MD_SIMD_INT nlocal_vec = simd_i32_broadcast(nlocal);
    MD_SIMD_INT npairs_vec = simd_i32_broadcast(npairs);

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        MD_FLOAT* tf = fbuf + (size_t)tid * nall * 3;
        for (int k = 0; k < nall * 3; k++) {
            tf[k] = 0.0;
        }

#pragma omp for schedule(static)
        for (int p = 0; p < npairs; p += VECTOR_WIDTH) {
            MD_SIMD_MASK valid_mask = simd_mask_i32_cond_lt(
                simd_i32_add(simd_i32_broadcast(p), simd_i32_seq()), npairs_vec);

            MD_SIMD_INT ivec = simd_i32_load(&neighbor->pair_i[p]);
            MD_SIMD_INT jvec = simd_i32_load(&neighbor->pair_j[p]);

#if LJ_COMB_RULE == LJ_COMB_GEOM
            MD_SIMD_FLOAT sqrt_eps_i = simd_real_gather(ivec,
                atom->sqrt_epsilon,
                sizeof(MD_FLOAT));
            MD_SIMD_FLOAT sigma3_i   = simd_real_gather(ivec, atom->sigma3, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(jvec,
                atom->sqrt_epsilon,
                sizeof(MD_FLOAT));
            MD_SIMD_FLOAT sigma3_j   = simd_real_gather(jvec, atom->sigma3, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
            MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_NONE
            MD_SIMD_INT ti           = simd_i32_gather(ivec, atom->type, sizeof(int));
            MD_SIMD_INT tj           = simd_i32_gather(jvec, atom->type, sizeof(int));
            MD_SIMD_INT tij          = simd_i32_add(simd_i32_mul(ti, ntypes_vec), tj);
            MD_SIMD_FLOAT sigma6_vec = simd_real_gather(tij, atom->sigma6, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT eps_vec    = simd_real_gather(tij, atom->epsilon, sizeof(MD_FLOAT));
#endif

#ifdef ATOM_POSITION_AOS
            MD_SIMD_INT i3     = simd_i32_add(simd_i32_add(ivec, ivec), ivec);
            MD_SIMD_INT j3     = simd_i32_add(simd_i32_add(jvec, jvec), jvec);
            MD_SIMD_FLOAT xi   = simd_real_gather(i3, &(atom->x[0]), sizeof(MD_FLOAT));
            MD_SIMD_FLOAT yi   = simd_real_gather(i3, &(atom->x[1]), sizeof(MD_FLOAT));
            MD_SIMD_FLOAT zi   = simd_real_gather(i3, &(atom->x[2]), sizeof(MD_FLOAT));
            MD_SIMD_FLOAT xj   = simd_real_gather(j3, &(atom->x[0]), sizeof(MD_FLOAT));
            MD_SIMD_FLOAT yj   = simd_real_gather(j3, &(atom->x[1]), sizeof(MD_FLOAT));
            MD_SIMD_FLOAT zj   = simd_real_gather(j3, &(atom->x[2]), sizeof(MD_FLOAT));
#else
            MD_SIMD_FLOAT xi = simd_real_gather(ivec, atom->x, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT yi = simd_real_gather(ivec, atom->y, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT zi = simd_real_gather(ivec, atom->z, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT xj = simd_real_gather(jvec, atom->x, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT yj = simd_real_gather(jvec, atom->y, sizeof(MD_FLOAT));
            MD_SIMD_FLOAT zj = simd_real_gather(jvec, atom->z, sizeof(MD_FLOAT));
#endif

            MD_SIMD_FLOAT delx = simd_real_sub(xi, xj);
            MD_SIMD_FLOAT dely = simd_real_sub(yi, yj);
            MD_SIMD_FLOAT delz = simd_real_sub(zi, zj);
            MD_SIMD_FLOAT rsq  = simd_real_fma(delx,
                delx,
                simd_real_fma(dely, dely, simd_real_mul(delz, delz)));
            MD_SIMD_MASK cutoff_mask = simd_mask_and(valid_mask,
                simd_mask_cond_lt(rsq, cutforcesq_vec));

            MD_SIMD_FLOAT sr2   = simd_real_reciprocal(rsq);
            MD_SIMD_FLOAT sr6   = simd_real_mul(sr2,
                simd_real_mul(sr2, simd_real_mul(sr2, sigma6_vec)));
            MD_SIMD_FLOAT force = simd_real_mul(c48_vec,
                simd_real_mul(sr6,
                    simd_real_mul(simd_real_sub(sr6, c05_vec), simd_real_mul(sr2, eps_vec))));

            // Zero out invalid/failed-cutoff lanes so the scalar write-back
            // below can add unconditionally without a per-lane branch.
            MD_SIMD_FLOAT fx_masked = simd_real_masked_add(simd_real_zero(),
                simd_real_mul(delx, force),
                cutoff_mask);
            MD_SIMD_FLOAT fy_masked = simd_real_masked_add(simd_real_zero(),
                simd_real_mul(dely, force),
                cutoff_mask);
            MD_SIMD_FLOAT fz_masked = simd_real_masked_add(simd_real_zero(),
                simd_real_mul(delz, force),
                cutoff_mask);

            MD_FLOAT fxa[VECTOR_WIDTH] __attribute__((aligned(64)));
            MD_FLOAT fya[VECTOR_WIDTH] __attribute__((aligned(64)));
            MD_FLOAT fza[VECTOR_WIDTH] __attribute__((aligned(64)));
            simd_real_store(fxa, fx_masked);
            simd_real_store(fya, fy_masked);
            simd_real_store(fza, fz_masked);

            for (int lane = 0; lane < VECTOR_WIDTH; lane++) {
                int li = neighbor->pair_i[p + lane];
                int lj = neighbor->pair_j[p + lane];

                tf[li * 3 + 0] += fxa[lane];
                tf[li * 3 + 1] += fya[lane];
                tf[li * 3 + 2] += fza[lane];

                if ((param->half_neigh && lj < nlocal) || param->method) {
                    tf[lj * 3 + 0] -= fxa[lane];
                    tf[lj * 3 + 1] -= fya[lane];
                    tf[lj * 3 + 2] -= fza[lane];
                }
            }
        }

#pragma omp for schedule(static)
        for (int i = 0; i < nall; i++) {
            MD_FLOAT sx = 0.0, sy = 0.0, sz = 0.0;
            for (int t = 0; t < nthreads; t++) {
                MD_FLOAT* tfi = fbuf + (size_t)t * nall * 3;
                sx += tfi[i * 3 + 0];
                sy += tfi[i * 3 + 1];
                sz += tfi[i * 3 + 2];
            }
            atom_fx(i) += sx;
            atom_fy(i) += sy;
            atom_fz(i) += sz;
        }

        LIKWID_MARKER_STOP("force");
    }

    addStat(stats->total_force_neighs, npairs);

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJPairList_simd end\n");
    return E - S;
}
#endif
