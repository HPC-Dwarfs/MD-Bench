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

// Compile-time guards for unsupported configurations
#ifdef NBLIST_SOA
#error "SIMD kernel not implemented when NBLIST_DATA_LAYOUT is SOA"
#endif

double computeForceLJFullNeigh_simd(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    int Nlocal = atom->Nlocal;
    int* neighs;
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
            neighs                    = &neighbor->neighbors[i * neighbor->maxneighs];
            int numneighs             = neighbor->numneigh[i];
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
                MD_SIMD_INT j            = simd_i32_mask_load(&neighs[k], mask_numneighs);

#if LJ_COMB_RULE == LJ_COMB_GEOM
                // Direct gather of per-atom LJ params (avoids type index computation)
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j, atom->sqrt_epsilon, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j, atom->sigma3, sizeof(MD_FLOAT));
                // Geometric combination: eps_ij = sqrt(eps_i) * sqrt(eps_j), sigma6_ij = sigma3_i * sigma3_j
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_NONE
                MD_SIMD_INT tj           = simd_i32_gather(j, atom->type, sizeof(int));
                MD_SIMD_INT tij          = simd_i32_add(tbase_i, tj);
                MD_SIMD_FLOAT sigma6_vec = simd_real_gather(tij, atom->sigma6, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_gather(tij, atom->epsilon, sizeof(MD_FLOAT));
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
                MD_SIMD_FLOAT delx = xtmp -
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
    int* neighs;
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
            neighs                    = &neighbor->neighbors[i * neighbor->maxneighs];
            int numneighs             = neighbor->numneigh[i];
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
                MD_SIMD_INT j            = simd_i32_mask_load(&neighs[k], mask_numneighs);

#if LJ_COMB_RULE == LJ_COMB_GEOM
                // Direct gather of per-atom LJ params (avoids type index computation)
                MD_SIMD_FLOAT sqrt_eps_j = simd_real_gather(j, atom->sqrt_epsilon, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma3_j   = simd_real_gather(j, atom->sigma3, sizeof(MD_FLOAT));
                // Geometric combination: eps_ij = sqrt(eps_i) * sqrt(eps_j), sigma6_ij = sigma3_i * sigma3_j
                MD_SIMD_FLOAT eps_vec    = simd_real_mul(sqrt_eps_i, sqrt_eps_j);
                MD_SIMD_FLOAT sigma6_vec = simd_real_mul(sigma3_i, sigma3_j);
#elif LJ_COMB_RULE == LJ_COMB_NONE
                MD_SIMD_INT tj           = simd_i32_gather(j, atom->type, sizeof(int));
                MD_SIMD_INT tij          = simd_i32_add(tbase_i, tj);
                MD_SIMD_FLOAT sigma6_vec = simd_real_gather(tij, atom->sigma6, sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps_vec    = simd_real_gather(tij, atom->epsilon, sizeof(MD_FLOAT));
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
                MD_SIMD_FLOAT delx = xtmp -
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
