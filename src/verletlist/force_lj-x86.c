/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
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

double computeForceLJFullNeigh_simd(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
#ifdef NBLIST_SOA
    fprintf(stderr, "Error: SIMD kernel not implemented when NBLIST_DATA_LAYOUT is SOA!");
    exit(-1);
#endif
    int Nlocal = atom->Nlocal;
    int* neighs;
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;

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
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
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

            for (int k = 0; k < numneighs; k += VECTOR_WIDTH) {
                // If the last iteration of this loop is separated from the rest, this
                // mask can be set only there
                MD_SIMD_MASK mask_numneighs = simd_mask_i32_cond_lt(
                    simd_i32_add(simd_i32_broadcast(k), simd_i32_seq()),
                    numneighs_vec);
                MD_SIMD_INT j            = simd_i32_mask_load(&neighs[k], mask_numneighs);
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
