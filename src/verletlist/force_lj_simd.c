/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 *
 * SIMD-optimized Lennard-Jones force kernels for Verlet Lists
 * Supports: AVX2, AVX512, NEON, SVE (double precision)
 * Requires: __SIMD_KERNEL__ flag. Works under NBLIST_AOS, NBLIST_SOA, and
 * NBLIST_CSR neighbor-list layouts (see the per-layout neighbor-index load
 * below, and neighbor.c for the corresponding buffer-safety guarantees each
 * layout relies on).
 *
 * LJ combination rules (compile-time via -DLJ_COMB_RULE=<value>):
 *   LJ_COMB_SINGLE (0): Single atom type - broadcast global epsilon/sigma
 *   LJ_COMB_GEOM   (1): Geometric - sqrt(eps_i*eps_j), sigma3_i*sigma3_j
 *   LJ_COMB_FULL   (2): Full type-pair matrix lookup via gather
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

#ifdef __SIMD_KERNEL__
#include <simd.h>
#endif

// NBLIST_AOS/NBLIST_CSR are both contiguous per-atom (neighbors[i*maxneighs+j]
// and neighbors[neigh_start[i]+j] respectively), so they share the same
// neighs[k] pointer-walk below. NBLIST_SOA is strided (neighbors[j*Nlocal+i]),
// so it instead computes a per-lane gather index directly (clamped to
// maxneighs-1, which neighbor.c's shared AOS/SOA buffer sizing guarantees is
// always in-bounds regardless of ghost-atom count) and gathers with
// simd_i32_gather. Both CSR's buffer (padded by VECTOR_WIDTH zeroed sentinel
// elements) and the shared AOS/SOA buffer (fully zeroed on allocation) are
// safe for a full-width masked load/gather past an atom's actual neighbor
// count -- see neighbor.c.

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
#elif !defined(NBLIST_SOA)
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            int numneighs_aligned     = numneighs - (numneighs % VECTOR_WIDTH);
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
#ifdef NBLIST_SOA
            MD_SIMD_INT nlocal_vec_soa   = simd_i32_broadcast(Nlocal);
            MD_SIMD_INT i_vec_soa        = simd_i32_broadcast(i);
            MD_SIMD_INT maxneighs_m1_vec = simd_i32_broadcast(neighbor->maxneighs - 1);
#endif
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
#elif LJ_COMB_RULE == LJ_COMB_FULL
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            // Aligned main loop: k+VECTOR_WIDTH-1 < numneighs always holds here, so
            // no bounds mask is needed -- only the cutoff test gates accumulation.
            for (int k = 0; k < numneighs_aligned; k += VECTOR_WIDTH) {
#ifdef NBLIST_SOA
                MD_SIMD_INT k_vec   = simd_i32_add(simd_i32_broadcast(k), simd_i32_seq());
                MD_SIMD_INT soa_idx = simd_i32_add(
                    simd_i32_mul(k_vec, nlocal_vec_soa), i_vec_soa);
                MD_SIMD_INT j       = simd_i32_gather(
                    soa_idx, neighbor->neighbors, sizeof(int));
#else
                MD_SIMD_INT j = simd_i32_loadu(&neighs[k]);
#endif

#if LJ_COMB_RULE == LJ_COMB_GEOM
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j,
                    atom->sqrt_epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j,
                    atom->sigma3,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_FULL
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
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp, simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely       = simd_real_sub(ytmp, simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz       = simd_real_sub(ztmp, simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT)));
#else
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp,
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely = simd_real_sub(ytmp,
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz = simd_real_sub(ztmp,
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT)));
#endif
                MD_SIMD_FLOAT rsq        = simd_real_fma(delx,
                    delx,
                    simd_real_fma(dely, dely, simd_real_mul(delz, delz)));
                MD_SIMD_MASK cutoff_mask = simd_mask_cond_lt(rsq, cutforcesq_vec);
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

            // Tail: whatever is left below one full VECTOR_WIDTH, masked.
            if (numneighs_aligned < numneighs) {
                int k = numneighs_aligned;
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
#ifdef NBLIST_SOA
                MD_SIMD_INT k_safe  = simd_i32_min(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()), maxneighs_m1_vec);
                MD_SIMD_INT soa_idx = simd_i32_add(
                    simd_i32_mul(k_safe, nlocal_vec_soa), i_vec_soa);
                MD_SIMD_INT j       = simd_i32_gather(
                    soa_idx, neighbor->neighbors, sizeof(int));
#else
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);
#endif

#if LJ_COMB_RULE == LJ_COMB_GEOM
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j,
                    atom->sqrt_epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j,
                    atom->sigma3,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_FULL
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
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp, simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely       = simd_real_sub(ytmp, simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz       = simd_real_sub(ztmp, simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT)));
#else
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp,
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely = simd_real_sub(ytmp,
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz = simd_real_sub(ztmp,
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT)));
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
#elif !defined(NBLIST_SOA)
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            int numneighs_aligned     = numneighs - (numneighs % VECTOR_WIDTH);
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
#ifdef NBLIST_SOA
            MD_SIMD_INT nlocal_vec_soa   = simd_i32_broadcast(Nlocal);
            MD_SIMD_INT i_vec_soa        = simd_i32_broadcast(i);
            MD_SIMD_INT maxneighs_m1_vec = simd_i32_broadcast(neighbor->maxneighs - 1);
#endif
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
#elif LJ_COMB_RULE == LJ_COMB_FULL
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            // Aligned main loop: k+VECTOR_WIDTH-1 < numneighs always holds here, so
            // no bounds mask is needed -- only the cutoff test gates accumulation.
            for (int k = 0; k < numneighs_aligned; k += VECTOR_WIDTH) {
#ifdef NBLIST_SOA
                MD_SIMD_INT k_vec   = simd_i32_add(simd_i32_broadcast(k), simd_i32_seq());
                MD_SIMD_INT soa_idx = simd_i32_add(
                    simd_i32_mul(k_vec, nlocal_vec_soa), i_vec_soa);
                MD_SIMD_INT j       = simd_i32_gather(
                    soa_idx, neighbor->neighbors, sizeof(int));
#else
                MD_SIMD_INT j = simd_i32_loadu(&neighs[k]);
#endif

#if LJ_COMB_RULE == LJ_COMB_GEOM
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j,
                    atom->sqrt_epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j,
                    atom->sigma3,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_FULL
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
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp, simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely       = simd_real_sub(ytmp, simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz       = simd_real_sub(ztmp, simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT)));
#else
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp,
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely = simd_real_sub(ytmp,
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz = simd_real_sub(ztmp,
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT)));
#endif
                MD_SIMD_FLOAT rsq        = simd_real_fma(delx,
                    delx,
                    simd_real_fma(dely, dely, simd_real_mul(delz, delz)));
                MD_SIMD_MASK cutoff_mask = simd_mask_cond_lt(rsq, cutforcesq_vec);
                MD_SIMD_FLOAT sr2        = simd_real_reciprocal(rsq);
                MD_SIMD_FLOAT sr6        = simd_real_mul(sr2,
                    simd_real_mul(sr2, simd_real_mul(sr2, sigma6_vec)));
                MD_SIMD_FLOAT force      = simd_real_mul(c48_vec,
                    simd_real_mul(sr6,
                        simd_real_mul(simd_real_sub(sr6, c05_vec),
                            simd_real_mul(sr2, eps_vec))));

                MD_SIMD_FLOAT fx_tmp = simd_real_mul(delx, force);
                MD_SIMD_FLOAT fy_tmp = simd_real_mul(dely, force);
                MD_SIMD_FLOAT fz_tmp = simd_real_mul(delz, force);

                fix = simd_real_masked_add(fix, fx_tmp, cutoff_mask);
                fiy = simd_real_masked_add(fiy, fy_tmp, cutoff_mask);
                fiz = simd_real_masked_add(fiz, fz_tmp, cutoff_mask);

                // Newton's third law via vectorized scatter -- not thread-safe
                // under OpenMP (no atomic scatter support).
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

            // Tail: whatever is left below one full VECTOR_WIDTH, masked.
            if (numneighs_aligned < numneighs) {
                int k = numneighs_aligned;
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
#ifdef NBLIST_SOA
                MD_SIMD_INT k_safe  = simd_i32_min(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()), maxneighs_m1_vec);
                MD_SIMD_INT soa_idx = simd_i32_add(
                    simd_i32_mul(k_safe, nlocal_vec_soa), i_vec_soa);
                MD_SIMD_INT j       = simd_i32_gather(
                    soa_idx, neighbor->neighbors, sizeof(int));
#else
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);
#endif

#if LJ_COMB_RULE == LJ_COMB_GEOM
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j,
                    atom->sqrt_epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j,
                    atom->sigma3,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_FULL
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
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp, simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely       = simd_real_sub(ytmp, simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz       = simd_real_sub(ztmp, simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT)));
#else
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp,
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely = simd_real_sub(ytmp,
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz = simd_real_sub(ztmp,
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT)));
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

                MD_SIMD_FLOAT fx_tmp = simd_real_mul(delx, force);
                MD_SIMD_FLOAT fy_tmp = simd_real_mul(dely, force);
                MD_SIMD_FLOAT fz_tmp = simd_real_mul(delz, force);

                fix = simd_real_masked_add(fix, fx_tmp, cutoff_mask);
                fiy = simd_real_masked_add(fiy, fy_tmp, cutoff_mask);
                fiz = simd_real_masked_add(fiz, fz_tmp, cutoff_mask);

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

// Extra per-i (and, for LJ_COMB_FULL, per-pair) state that ljforce_compressed() needs
// on top of the buffered (j, delx, dely, delz) triples, depending on LJ_COMB_RULE.
#if LJ_COMB_RULE == LJ_COMB_GEOM
#define LJ_EXTRA_PARAMS MD_SIMD_FLOAT sqrt_eps_i, MD_SIMD_FLOAT sigma3_i
#define LJ_EXTRA_ARGS   sqrt_eps_i, sigma3_i
#elif LJ_COMB_RULE == LJ_COMB_FULL
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
#elif LJ_COMB_RULE == LJ_COMB_FULL
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
#elif !defined(NBLIST_SOA)
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
#ifdef NBLIST_SOA
            MD_SIMD_INT nlocal_vec_soa   = simd_i32_broadcast(Nlocal);
            MD_SIMD_INT i_vec_soa        = simd_i32_broadcast(i);
            MD_SIMD_INT maxneighs_m1_vec = simd_i32_broadcast(neighbor->maxneighs - 1);
#endif
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
#elif LJ_COMB_RULE == LJ_COMB_FULL
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs; k += VECTOR_WIDTH) {
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
#ifdef NBLIST_SOA
                MD_SIMD_INT k_safe  = simd_i32_min(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()), maxneighs_m1_vec);
                MD_SIMD_INT soa_idx = simd_i32_add(
                    simd_i32_mul(k_safe, nlocal_vec_soa), i_vec_soa);
                MD_SIMD_INT j       = simd_i32_gather(
                    soa_idx, neighbor->neighbors, sizeof(int));
#else
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);
#endif

#ifdef ATOM_POSITION_AOS
                MD_SIMD_INT j3           = simd_i32_add(simd_i32_add(j, j), j);
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp, simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely       = simd_real_sub(ytmp, simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz       = simd_real_sub(ztmp, simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT)));
#else
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp,
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely = simd_real_sub(ytmp,
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz = simd_real_sub(ztmp,
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT)));
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
#elif !defined(NBLIST_SOA)
            int* neighs               = &neighbor->neighbors[i * neighbor->maxneighs];
#endif
            int numneighs             = neighbor->numneigh_inner[i];
            MD_SIMD_INT numneighs_vec = simd_i32_broadcast(numneighs);
#ifdef NBLIST_SOA
            MD_SIMD_INT nlocal_vec_soa   = simd_i32_broadcast(Nlocal);
            MD_SIMD_INT i_vec_soa        = simd_i32_broadcast(i);
            MD_SIMD_INT maxneighs_m1_vec = simd_i32_broadcast(neighbor->maxneighs - 1);
#endif
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
#elif LJ_COMB_RULE == LJ_COMB_FULL
            MD_SIMD_INT tbase_i = simd_i32_broadcast(atom->type[i] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs; k += VECTOR_WIDTH) {
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
#ifdef NBLIST_SOA
                MD_SIMD_INT k_safe  = simd_i32_min(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()), maxneighs_m1_vec);
                MD_SIMD_INT soa_idx = simd_i32_add(
                    simd_i32_mul(k_safe, nlocal_vec_soa), i_vec_soa);
                MD_SIMD_INT j       = simd_i32_gather(
                    soa_idx, neighbor->neighbors, sizeof(int));
#else
                MD_SIMD_INT j = simd_i32_mask_load(&neighs[k], mask_numneighs);
#endif

#ifdef ATOM_POSITION_AOS
                MD_SIMD_INT j3           = simd_i32_add(simd_i32_add(j, j), j);
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp, simd_real_gather(j3,
                                                &(atom->x[0]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely       = simd_real_sub(ytmp, simd_real_gather(j3,
                                                &(atom->x[1]),
                                                sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz       = simd_real_sub(ztmp, simd_real_gather(j3,
                                                &(atom->x[2]),
                                                sizeof(MD_FLOAT)));
#else
                MD_SIMD_FLOAT delx       = simd_real_sub(xtmp,
                                     simd_real_gather(j, atom->x, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT dely = simd_real_sub(ytmp,
                                     simd_real_gather(j, atom->y, sizeof(MD_FLOAT)));
                MD_SIMD_FLOAT delz = simd_real_sub(ztmp,
                                     simd_real_gather(j, atom->z, sizeof(MD_FLOAT)));
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
