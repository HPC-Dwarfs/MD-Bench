/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>

#include <atom.h>
#include <force.h>
#include <likwid-marker.h>
#include <math.h>
#include <neighbor.h>
#include <parameter.h>
#include <simd.h>
#include <stats.h>
#include <timing.h>
#include <util.h>

void computeForceGhostShell(Parameter*, Atom*, Neighbor*);

#ifdef USE_REFERENCE_KERNEL
double computeForceLJRef(Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ begin\n");
    const int nbM = atom->Nclusters_local;
    const int nbN = neighbor->maxneighs;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#endif

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ci_f  = &atom->cl_f[ci_vec_base];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ci_f[CL_X_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Y_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Z_INDEX_3D(cii)] = 0.0;
        }
    }

    for (int cg = atom->ncj; cg < atom->ncj + atom->Nclusters_ghost; cg++) {
        int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cg);
        MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];
        for (int cjj = 0; cjj < atom->jclusters[cg].natoms; cjj++) {
            cj_f[CL_X_INDEX_3D(cjj)] = 0.0;
            cj_f[CL_Y_INDEX_3D(cjj)] = 0.0;
            cj_f[CL_Z_INDEX_3D(cjj)] = 0.0;
        }
    }

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_cj0      = CJ0_FROM_CI(ci);
            int ci_cj1      = CJ1_FROM_CI(ci);
            int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x  = &atom->cl_x[ci_vec_base];
            MD_FLOAT* ci_f  = &atom->cl_f[ci_vec_base];
            int numneighs   = neighbor->numneigh[ci];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int ci_sca_base = CI_SCALAR_BASE_INDEX(ci);
            int* ci_t       = &atom->cl_t[ci_sca_base];
#endif

            for (int k = 0; k < numneighs; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                int any         = 0;
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
                MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                for (int cii = 0; cii < CLUSTER_M; cii++) {
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                    int type_i = ci_t[cii];
#endif
                    MD_FLOAT xtmp = ci_x[CL_X_INDEX_3D(cii)];
                    MD_FLOAT ytmp = ci_x[CL_Y_INDEX_3D(cii)];
                    MD_FLOAT ztmp = ci_x[CL_Z_INDEX_3D(cii)];
                    MD_FLOAT fix  = 0;
                    MD_FLOAT fiy  = 0;
                    MD_FLOAT fiz  = 0;

                    for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
                        int cond;
#if CLUSTER_M == CLUSTER_N
                        cond = neighbor->half_neigh ? (ci_cj0 != cj || cii < cjj)
                                                    : (ci_cj0 != cj || cii != cjj);
#elif CLUSTER_M < CLUSTER_N
                        cond = neighbor->half_neigh
                                   ? (ci_cj0 != cj || cii + CLUSTER_M * (ci & 0x1) < cjj)
                                   : (ci_cj0 != cj ||
                                         cii + CLUSTER_M * (ci & 0x1) != cjj);
#else
                        cond = neighbor->half_neigh
                                   ? (ci_cj0 != cj || cii < cjj) &&
                                         (ci_cj1 != cj || cii < cjj + CLUSTER_N)
                                   : (ci_cj0 != cj || cii != cjj) &&
                                         (ci_cj1 != cj || cii != cjj + CLUSTER_N);
#endif
                        if (cond) {
                            MD_FLOAT delx = xtmp - cj_x[CL_X_INDEX_3D(cjj)];
                            MD_FLOAT dely = ytmp - cj_x[CL_Y_INDEX_3D(cjj)];
                            MD_FLOAT delz = ztmp - cj_x[CL_Z_INDEX_3D(cjj)];
                            MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                            int type_j          = cj_t[cjj];
                            int type_index      = type_i * atom->ntypes + type_j;
                            MD_FLOAT cutforcesq = atom->cutforcesq[type_index];
                            MD_FLOAT sigma6     = atom->sigma6[type_index];
                            MD_FLOAT epsilon    = atom->epsilon[type_index];
#endif

                            if (rsq < cutforcesq) {
                                MD_FLOAT sr2   = 1.0 / rsq;
                                MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
                                MD_FLOAT force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;

                                if (neighbor->half_neigh || param->method) {
                                    cj_f[CL_X_INDEX_3D(cjj)] -= delx * force;
                                    cj_f[CL_Y_INDEX_3D(cjj)] -= dely * force;
                                    cj_f[CL_Z_INDEX_3D(cjj)] -= delz * force;
                                }

                                fix += delx * force;
                                fiy += dely * force;
                                fiz += delz * force;
                                any = 1;
                                addStat(stats->atoms_within_cutoff, 1);
                            } else {
                                addStat(stats->atoms_outside_cutoff, 1);
                            }
                        }
                    }

                    if (any != 0) {
                        addStat(stats->clusters_within_cutoff, 1);
                    } else {
                        addStat(stats->clusters_outside_cutoff, 1);
                    }

                    ci_f[CL_X_INDEX_3D(cii)] += fix;
                    ci_f[CL_Y_INDEX_3D(cii)] += fiy;
                    ci_f[CL_Z_INDEX_3D(cii)] += fiz;
                }
            }

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters,
                (long long int)((double)numneighs * CLUSTER_M / CLUSTER_N));
        }
        if (param->method == eightShell) {
            computeForceGhostShell(param, atom, neighbor);
        }
        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ end\n");
    return E - S;
}
#else
double computeForceLJ2xnnHalfNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_2xnn begin\n");
    const int nbM = atom->Nclusters_local;
    const int nbN = neighbor->maxneighs;
    MD_FLOAT cutforcesq          = param->cutforce * param->cutforce;
    MD_FLOAT sigma6              = param->sigma6;
    MD_FLOAT epsilon             = param->epsilon;
    MD_SIMD_FLOAT c48_vec        = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec        = simd_real_broadcast(0.5);

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
#endif
    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ci_f  = &atom->cl_f[ci_vec_base];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ci_f[CL_X_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Y_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Z_INDEX_3D(cii)] = 0.0;
        }
    }

    for (int cg = atom->ncj; cg < atom->ncj + atom->Nclusters_ghost; cg++) {
        int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cg);
        MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];
        for (int cjj = 0; cjj < atom->jclusters[cg].natoms; cjj++) {
            cj_f[CL_X_INDEX_3D(cjj)] = 0.0;
            cj_f[CL_Y_INDEX_3D(cjj)] = 0.0;
            cj_f[CL_Z_INDEX_3D(cjj)] = 0.0;
        }
    }

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_cj0           = CJ0_FROM_CI(ci);
#if CLUSTER_M > CLUSTER_N
            int ci_cj1           = CJ1_FROM_CI(ci);
#endif
            int ci_vec_base      = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x       = &atom->cl_x[ci_vec_base];
            MD_FLOAT* ci_f       = &atom->cl_f[ci_vec_base];
            int numneighs        = neighbor->numneigh[ci];
            int numneighs_masked = neighbor->numneigh_masked[ci];

            MD_SIMD_FLOAT xi0_tmp = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(0)]);
            MD_SIMD_FLOAT xi2_tmp = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(2)]);
            MD_SIMD_FLOAT yi0_tmp = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(0)]);
            MD_SIMD_FLOAT yi2_tmp = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(2)]);
            MD_SIMD_FLOAT zi0_tmp = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(0)]);
            MD_SIMD_FLOAT zi2_tmp = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(2)]);
            MD_SIMD_FLOAT fix0    = simd_real_zero();
            MD_SIMD_FLOAT fiy0    = simd_real_zero();
            MD_SIMD_FLOAT fiz0    = simd_real_zero();
            MD_SIMD_FLOAT fix2    = simd_real_zero();
            MD_SIMD_FLOAT fiy2    = simd_real_zero();
            MD_SIMD_FLOAT fiz2    = simd_real_zero();

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int ci_sca_base       = CI_SCALAR_BASE_INDEX(ci);
            int* ci_t             = &atom->cl_t[ci_sca_base];
            MD_SIMD_INT tbase0    = simd_i32_load_h_dual_scaled(&ci_t[0], atom->ntypes);
            MD_SIMD_INT tbase2    = simd_i32_load_h_dual_scaled(&ci_t[2], atom->ntypes);
#endif

            for (int k = 0; k < numneighs_masked; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                // int imask = neighs_imask[k];
                MD_FLOAT* cj_x = &atom->cl_x[cj_vec_base];
                MD_FLOAT* cj_f = &atom->cl_f[cj_vec_base];
                // MD_SIMD_MASK interact0;
                // MD_SIMD_MASK interact2;

                // gmx_load_simd_2xnn_interactions((int)imask, filter0, filter2,
                // &interact0, &interact2);

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp    = simd_real_load_h_duplicate(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp    = simd_real_load_h_duplicate(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp    = simd_real_load_h_duplicate(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0     = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0     = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0     = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2     = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2     = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2     = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT rsq0      = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq2      = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));

#if CLUSTER_M == CLUSTER_N
                unsigned int cond0      = (unsigned int)(cj == ci_cj0);
                MD_SIMD_MASK excl_mask0 = simd_mask_from_u32(
                    atom->masks_2xnn_hn[cond0 * 2 + 0]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_2xnn_hn[cond0 * 2 + 1]);
#else
#if CLUSTER_M < CLUSTER_N
                unsigned int cond0      = (unsigned int)((cj << 1) + 0 == ci);
                unsigned int cond1      = (unsigned int)((cj << 1) + 1 == ci);
#else
                unsigned int cond0 = (unsigned int)(cj == ci_cj0);
                unsigned int cond1 = (unsigned int)(cj == ci_cj1);
#endif
                MD_SIMD_MASK excl_mask0 = simd_mask_from_u32(
                    atom->masks_2xnn_hn[cond0 * 4 + cond1 * 2 + 0]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_2xnn_hn[cond0 * 4 + cond1 * 2 + 1]);
#endif

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp = simd_i32_load_h_duplicate(cj_t);
                MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_0    = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2    = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps0        = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2        = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT sigma6_0    = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2    = sigma6_vec;
                MD_SIMD_FLOAT eps0        = eps_vec;
                MD_SIMD_FLOAT eps2        = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0, cutforcesq0);
                MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2, cutforcesq2);
                cutoff_mask0              = simd_mask_and(cutoff_mask0, excl_mask0);
                cutoff_mask2              = simd_mask_and(cutoff_mask2, excl_mask2);

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));

                MD_SIMD_FLOAT tx0 = simd_real_select_by_mask(simd_real_mul(delx0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT ty0 = simd_real_select_by_mask(simd_real_mul(dely0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tz0 = simd_real_select_by_mask(simd_real_mul(delz0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tx2 = simd_real_select_by_mask(simd_real_mul(delx2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT ty2 = simd_real_select_by_mask(simd_real_mul(dely2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT tz2 = simd_real_select_by_mask(simd_real_mul(delz2, force2),
                    cutoff_mask2);

                fix0 = simd_real_add(fix0, tx0);
                fiy0 = simd_real_add(fiy0, ty0);
                fiz0 = simd_real_add(fiz0, tz0);
                fix2 = simd_real_add(fix2, tx2);
                fiy2 = simd_real_add(fiy2, ty2);
                fiz2 = simd_real_add(fiz2, tz2);

                if (cj < CJ0_FROM_CI(atom->Nclusters_local) || param->method) {
                    simd_real_h_decr3(cj_f,
                        simd_real_add(tx0, tx2),
                        simd_real_add(ty0, ty2),
                        simd_real_add(tz0, tz2));
                }
            }

            for (int k = numneighs_masked; k < numneighs; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
                MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp = simd_real_load_h_duplicate(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp = simd_real_load_h_duplicate(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp = simd_real_load_h_duplicate(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0  = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0  = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0  = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2  = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2  = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2  = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT rsq0   = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq2   = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp   = simd_i32_load_h_duplicate(cj_t);
                MD_SIMD_INT tvec0    = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec2    = simd_i32_add(tbase2, tj_tmp);

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_0    = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2    = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps0        = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2        = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT sigma6_0    = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2    = sigma6_vec;
                MD_SIMD_FLOAT eps0        = eps_vec;
                MD_SIMD_FLOAT eps2        = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0, cutforcesq0);
                MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2, cutforcesq2);

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));

                MD_SIMD_FLOAT tx0 = simd_real_select_by_mask(simd_real_mul(delx0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT ty0 = simd_real_select_by_mask(simd_real_mul(dely0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tz0 = simd_real_select_by_mask(simd_real_mul(delz0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tx2 = simd_real_select_by_mask(simd_real_mul(delx2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT ty2 = simd_real_select_by_mask(simd_real_mul(dely2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT tz2 = simd_real_select_by_mask(simd_real_mul(delz2, force2),
                    cutoff_mask2);

                fix0 = simd_real_add(fix0, tx0);
                fiy0 = simd_real_add(fiy0, ty0);
                fiz0 = simd_real_add(fiz0, tz0);
                fix2 = simd_real_add(fix2, tx2);
                fiy2 = simd_real_add(fiy2, ty2);
                fiz2 = simd_real_add(fiz2, tz2);

                if (cj < CJ0_FROM_CI(atom->Nclusters_local) || param->method) {
                    simd_real_h_decr3(cj_f,
                        simd_real_add(tx0, tx2),
                        simd_real_add(ty0, ty2),
                        simd_real_add(tz0, tz2));
                }
            }

            simd_real_h_dual_incr_reduced_sum(&ci_f[CL_X_INDEX_3D(0)], fix0, fix2);
            simd_real_h_dual_incr_reduced_sum(&ci_f[CL_Y_INDEX_3D(0)], fiy0, fiy2);
            simd_real_h_dual_incr_reduced_sum(&ci_f[CL_Z_INDEX_3D(0)], fiz0, fiz2);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters,
                (long long int)((double)numneighs * CLUSTER_M / CLUSTER_N));
        }
        if (param->method == eightShell) {
            computeForceGhostShell(param, atom, neighbor);
        }
        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_2xnn end\n");
    return E - S;
}

double computeForceLJ2xnnFullNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_2xnn begin\n");
    const int nbM                = atom->Nclusters_local;
    const int nbN                = neighbor->maxneighs;
    MD_FLOAT cutforcesq          = param->cutforce * param->cutforce;
    MD_FLOAT sigma6              = param->sigma6;
    MD_FLOAT epsilon             = param->epsilon;
    MD_SIMD_FLOAT c48_vec        = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec        = simd_real_broadcast(0.5);

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
#endif

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ci_f  = &atom->cl_f[ci_vec_base];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ci_f[CL_X_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Y_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Z_INDEX_3D(cii)] = 0.0;
        }
    }

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_cj0           = CJ0_FROM_CI(ci);
#if CLUSTER_M > CLUSTER_N
            int ci_cj1           = CJ1_FROM_CI(ci);
#endif
            int ci_vec_base      = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x       = &atom->cl_x[ci_vec_base];
            MD_FLOAT* ci_f       = &atom->cl_f[ci_vec_base];
            const int numneighs        = neighbor->numneigh[ci];
            const int numneighs_masked = neighbor->numneigh_masked[ci];

            MD_SIMD_FLOAT xi0_tmp = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(0)]);
            MD_SIMD_FLOAT xi2_tmp = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(2)]);
            MD_SIMD_FLOAT yi0_tmp = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(0)]);
            MD_SIMD_FLOAT yi2_tmp = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(2)]);
            MD_SIMD_FLOAT zi0_tmp = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(0)]);
            MD_SIMD_FLOAT zi2_tmp = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(2)]);
            MD_SIMD_FLOAT fix0    = simd_real_zero();
            MD_SIMD_FLOAT fiy0    = simd_real_zero();
            MD_SIMD_FLOAT fiz0    = simd_real_zero();
            MD_SIMD_FLOAT fix2    = simd_real_zero();
            MD_SIMD_FLOAT fiy2    = simd_real_zero();
            MD_SIMD_FLOAT fiz2    = simd_real_zero();
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int ci_sca_base       = CI_SCALAR_BASE_INDEX(ci);
            int* ci_t             = &atom->cl_t[ci_sca_base];
            MD_SIMD_INT tbase0    = simd_i32_load_h_dual_scaled(&ci_t[0], atom->ntypes);
            MD_SIMD_INT tbase2    = simd_i32_load_h_dual_scaled(&ci_t[2], atom->ntypes);
#endif

            for (int k = 0; k < numneighs_masked; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
                unsigned int mask0, mask1, mask2, mask3;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp    = simd_real_load_h_duplicate(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp    = simd_real_load_h_duplicate(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp    = simd_real_load_h_duplicate(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0     = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0     = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0     = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2     = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2     = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2     = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT rsq0      = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq2      = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));

#if CLUSTER_M == CLUSTER_N
                unsigned int cond0      = (unsigned int)(cj == ci_cj0);
                MD_SIMD_MASK excl_mask0 = simd_mask_from_u32(
                    atom->masks_2xnn_fn[cond0 * 2 + 0]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_2xnn_fn[cond0 * 2 + 1]);
#else
#if CLUSTER_M < CLUSTER_N
                unsigned int cond0        = (unsigned int)((cj << 1) + 0 == ci);
                unsigned int cond1        = (unsigned int)((cj << 1) + 1 == ci);
#else
                unsigned int cond0 = (unsigned int)(cj == ci_cj0);
                unsigned int cond1 = (unsigned int)(cj == ci_cj1);
#endif
                MD_SIMD_MASK excl_mask0   = simd_mask_from_u32(
                    atom->masks_2xnn_fn[cond0 * 4 + cond1 * 2 + 0]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_2xnn_fn[cond0 * 4 + cond1 * 2 + 1]);
#endif

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp = simd_i32_load_h_duplicate(cj_t);
                MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_0    = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2    = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps0        = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2        = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT sigma6_0    = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2    = sigma6_vec;
                MD_SIMD_FLOAT eps0        = eps_vec;
                MD_SIMD_FLOAT eps2        = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_and(excl_mask0,
                    simd_mask_cond_lt(rsq0, cutforcesq0));
                MD_SIMD_MASK cutoff_mask2 = simd_mask_and(excl_mask2,
                    simd_mask_cond_lt(rsq2, cutforcesq2));

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));

                fix0 = simd_real_masked_add(fix0,
                    simd_real_mul(delx0, force0),
                    cutoff_mask0);
                fiy0 = simd_real_masked_add(fiy0,
                    simd_real_mul(dely0, force0),
                    cutoff_mask0);
                fiz0 = simd_real_masked_add(fiz0,
                    simd_real_mul(delz0, force0),
                    cutoff_mask0);
                fix2 = simd_real_masked_add(fix2,
                    simd_real_mul(delx2, force2),
                    cutoff_mask2);
                fiy2 = simd_real_masked_add(fiy2,
                    simd_real_mul(dely2, force2),
                    cutoff_mask2);
                fiz2 = simd_real_masked_add(fiz2,
                    simd_real_mul(delz2, force2),
                    cutoff_mask2);
            }

            for (int k = numneighs_masked; k < numneighs; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp = simd_real_load_h_duplicate(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp = simd_real_load_h_duplicate(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp = simd_real_load_h_duplicate(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0  = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0  = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0  = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2  = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2  = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2  = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT rsq0   = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq2   = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp   = simd_i32_load_h_duplicate(cj_t);
                MD_SIMD_INT tvec0    = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec2    = simd_i32_add(tbase2, tj_tmp);

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_0    = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2    = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps0        = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2        = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT sigma6_0    = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2    = sigma6_vec;
                MD_SIMD_FLOAT eps0        = eps_vec;
                MD_SIMD_FLOAT eps2        = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0, cutforcesq0);
                MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2, cutforcesq2);

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));

                fix0 = simd_real_masked_add(fix0,
                    simd_real_mul(delx0, force0),
                    cutoff_mask0);
                fiy0 = simd_real_masked_add(fiy0,
                    simd_real_mul(dely0, force0),
                    cutoff_mask0);
                fiz0 = simd_real_masked_add(fiz0,
                    simd_real_mul(delz0, force0),
                    cutoff_mask0);
                fix2 = simd_real_masked_add(fix2,
                    simd_real_mul(delx2, force2),
                    cutoff_mask2);
                fiy2 = simd_real_masked_add(fiy2,
                    simd_real_mul(dely2, force2),
                    cutoff_mask2);
                fiz2 = simd_real_masked_add(fiz2,
                    simd_real_mul(delz2, force2),
                    cutoff_mask2);
            }

            simd_real_h_dual_incr_reduced_sum(&ci_f[CL_X_INDEX_3D(0)], fix0, fix2);
            simd_real_h_dual_incr_reduced_sum(&ci_f[CL_Y_INDEX_3D(0)], fiy0, fiy2);
            simd_real_h_dual_incr_reduced_sum(&ci_f[CL_Z_INDEX_3D(0)], fiz0, fiz2);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters, (long long int)((double)numneighs));
            // addStat(stats->force_iters, (long long int)((double)numneighs * CLUSTER_M /
            // CLUSTER_N));
        }

        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_2xnn end\n");
    return E - S;
}

double computeForceLJ4xnHalfNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_4xn begin\n");
    const int nbM                = atom->Nclusters_local;
    const int nbN                = neighbor->maxneighs;
    MD_FLOAT cutforcesq          = param->cutforce * param->cutforce;
    MD_FLOAT sigma6              = param->sigma6;
    MD_FLOAT epsilon             = param->epsilon;
    MD_SIMD_FLOAT c48_vec        = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec        = simd_real_broadcast(0.5);

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
#endif

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ci_f  = &atom->cl_f[ci_vec_base];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ci_f[CL_X_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Y_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Z_INDEX_3D(cii)] = 0.0;
        }
    }

    for (int cg = atom->ncj; cg < atom->ncj + atom->Nclusters_ghost; cg++) {
        int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cg);
        MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];
        for (int cjj = 0; cjj < atom->jclusters[cg].natoms; cjj++) {
            cj_f[CL_X_INDEX_3D(cjj)] = 0.0;
            cj_f[CL_Y_INDEX_3D(cjj)] = 0.0;
            cj_f[CL_Z_INDEX_3D(cjj)] = 0.0;
        }
    }

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_cj0           = CJ0_FROM_CI(ci);
#if CLUSTER_M > CLUSTER_N
            int ci_cj1           = CJ1_FROM_CI(ci);
#endif
            int ci_vec_base      = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x       = &atom->cl_x[ci_vec_base];
            MD_FLOAT* ci_f       = &atom->cl_f[ci_vec_base];
            int numneighs        = neighbor->numneigh[ci];
            int numneighs_masked = neighbor->numneigh_masked[ci];

            MD_SIMD_FLOAT xi0_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(0)]);
            MD_SIMD_FLOAT xi1_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(1)]);
            MD_SIMD_FLOAT xi2_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(2)]);
            MD_SIMD_FLOAT xi3_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(3)]);
            MD_SIMD_FLOAT yi0_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(0)]);
            MD_SIMD_FLOAT yi1_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(1)]);
            MD_SIMD_FLOAT yi2_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(2)]);
            MD_SIMD_FLOAT yi3_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(3)]);
            MD_SIMD_FLOAT zi0_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(0)]);
            MD_SIMD_FLOAT zi1_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(1)]);
            MD_SIMD_FLOAT zi2_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(2)]);
            MD_SIMD_FLOAT zi3_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(3)]);
            MD_SIMD_FLOAT fix0    = simd_real_zero();
            MD_SIMD_FLOAT fiy0    = simd_real_zero();
            MD_SIMD_FLOAT fiz0    = simd_real_zero();
            MD_SIMD_FLOAT fix1    = simd_real_zero();
            MD_SIMD_FLOAT fiy1    = simd_real_zero();
            MD_SIMD_FLOAT fiz1    = simd_real_zero();
            MD_SIMD_FLOAT fix2    = simd_real_zero();
            MD_SIMD_FLOAT fiy2    = simd_real_zero();
            MD_SIMD_FLOAT fiz2    = simd_real_zero();
            MD_SIMD_FLOAT fix3    = simd_real_zero();
            MD_SIMD_FLOAT fiy3    = simd_real_zero();
            MD_SIMD_FLOAT fiz3    = simd_real_zero();

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int ci_sca_base       = CI_SCALAR_BASE_INDEX(ci);
            int* ci_t             = &atom->cl_t[ci_sca_base];
            MD_SIMD_INT tbase0    = simd_i32_broadcast(ci_t[0] * atom->ntypes);
            MD_SIMD_INT tbase1    = simd_i32_broadcast(ci_t[1] * atom->ntypes);
            MD_SIMD_INT tbase2    = simd_i32_broadcast(ci_t[2] * atom->ntypes);
            MD_SIMD_INT tbase3    = simd_i32_broadcast(ci_t[3] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs_masked; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
                MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp    = simd_real_load(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp    = simd_real_load(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp    = simd_real_load(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0     = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0     = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0     = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx1     = simd_real_sub(xi1_tmp, xj_tmp);
                MD_SIMD_FLOAT dely1     = simd_real_sub(yi1_tmp, yj_tmp);
                MD_SIMD_FLOAT delz1     = simd_real_sub(zi1_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2     = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2     = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2     = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT delx3     = simd_real_sub(xi3_tmp, xj_tmp);
                MD_SIMD_FLOAT dely3     = simd_real_sub(yi3_tmp, yj_tmp);
                MD_SIMD_FLOAT delz3     = simd_real_sub(zi3_tmp, zj_tmp);

#if CLUSTER_M == CLUSTER_N
                unsigned int cond0      = (unsigned int)(cj == ci_cj0);
                MD_SIMD_MASK excl_mask0 = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 4 + 0]);
                MD_SIMD_MASK excl_mask1 = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 4 + 1]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 4 + 2]);
                MD_SIMD_MASK excl_mask3 = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 4 + 3]);
#else
#if CLUSTER_M < CLUSTER_N
                unsigned int cond0        = (unsigned int)((cj << 1) + 0 == ci);
                unsigned int cond1        = (unsigned int)((cj << 1) + 1 == ci);
#else
                unsigned int cond0 = (unsigned int)(cj == ci_cj0);
                unsigned int cond1 = (unsigned int)(cj == ci_cj1);
#endif
                MD_SIMD_MASK excl_mask0   = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 0]);
                MD_SIMD_MASK excl_mask1 = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 1]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 2]);
                MD_SIMD_MASK excl_mask3 = simd_mask_from_u32(
                    atom->masks_4xn_hn[cond0 * 8 + cond1 * 4 + 3]);
#endif

                MD_SIMD_FLOAT rsq0 = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq1 = simd_real_fma(delx1,
                    delx1,
                    simd_real_fma(dely1, dely1, simd_real_mul(delz1, delz1)));
                MD_SIMD_FLOAT rsq2 = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));
                MD_SIMD_FLOAT rsq3 = simd_real_fma(delx3,
                    delx3,
                    simd_real_fma(dely3, dely3, simd_real_mul(delz3, delz3)));

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp = simd_i32_load(cj_t);
                MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec1  = simd_i32_add(tbase1, tj_tmp);
                MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);
                MD_SIMD_INT tvec3  = simd_i32_add(tbase3, tj_tmp);

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq1 = simd_real_gather(tvec1,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq3 = simd_real_gather(tvec3,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT sigma6_0 = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_1 = simd_real_gather(tvec1,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2 = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_3 = simd_real_gather(tvec3,
                    atom->sigma6,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT eps0 = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps1 = simd_real_gather(tvec1,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2 = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps3 = simd_real_gather(tvec3,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq1 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq3 = cutforcesq_vec;

                MD_SIMD_FLOAT sigma6_0 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_1 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_3 = sigma6_vec;

                MD_SIMD_FLOAT eps0        = eps_vec;
                MD_SIMD_FLOAT eps1        = eps_vec;
                MD_SIMD_FLOAT eps2        = eps_vec;
                MD_SIMD_FLOAT eps3        = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_and(excl_mask0,
                    simd_mask_cond_lt(rsq0, cutforcesq0));
                MD_SIMD_MASK cutoff_mask1 = simd_mask_and(excl_mask1,
                    simd_mask_cond_lt(rsq1, cutforcesq1));
                MD_SIMD_MASK cutoff_mask2 = simd_mask_and(excl_mask2,
                    simd_mask_cond_lt(rsq2, cutforcesq2));
                MD_SIMD_MASK cutoff_mask3 = simd_mask_and(excl_mask3,
                    simd_mask_cond_lt(rsq3, cutforcesq3));

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_1 = simd_real_reciprocal(rsq1);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);
                MD_SIMD_FLOAT sr2_3 = simd_real_reciprocal(rsq3);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_1 = simd_real_mul(sr2_1,
                    simd_real_mul(sr2_1, simd_real_mul(sr2_1, sigma6_1)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));
                MD_SIMD_FLOAT sr6_3 = simd_real_mul(sr2_3,
                    simd_real_mul(sr2_3, simd_real_mul(sr2_3, sigma6_3)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force1 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_1,
                        simd_real_mul(simd_real_sub(sr6_1, c05_vec),
                            simd_real_mul(sr2_1, eps1))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));
                MD_SIMD_FLOAT force3 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_3,
                        simd_real_mul(simd_real_sub(sr6_3, c05_vec),
                            simd_real_mul(sr2_3, eps3))));

                MD_SIMD_FLOAT tx0 = simd_real_select_by_mask(simd_real_mul(delx0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT ty0 = simd_real_select_by_mask(simd_real_mul(dely0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tz0 = simd_real_select_by_mask(simd_real_mul(delz0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tx1 = simd_real_select_by_mask(simd_real_mul(delx1, force1),
                    cutoff_mask1);
                MD_SIMD_FLOAT ty1 = simd_real_select_by_mask(simd_real_mul(dely1, force1),
                    cutoff_mask1);
                MD_SIMD_FLOAT tz1 = simd_real_select_by_mask(simd_real_mul(delz1, force1),
                    cutoff_mask1);
                MD_SIMD_FLOAT tx2 = simd_real_select_by_mask(simd_real_mul(delx2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT ty2 = simd_real_select_by_mask(simd_real_mul(dely2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT tz2 = simd_real_select_by_mask(simd_real_mul(delz2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT tx3 = simd_real_select_by_mask(simd_real_mul(delx3, force3),
                    cutoff_mask3);
                MD_SIMD_FLOAT ty3 = simd_real_select_by_mask(simd_real_mul(dely3, force3),
                    cutoff_mask3);
                MD_SIMD_FLOAT tz3 = simd_real_select_by_mask(simd_real_mul(delz3, force3),
                    cutoff_mask3);

                fix0 = simd_real_add(fix0, tx0);
                fiy0 = simd_real_add(fiy0, ty0);
                fiz0 = simd_real_add(fiz0, tz0);
                fix1 = simd_real_add(fix1, tx1);
                fiy1 = simd_real_add(fiy1, ty1);
                fiz1 = simd_real_add(fiz1, tz1);
                fix2 = simd_real_add(fix2, tx2);
                fiy2 = simd_real_add(fiy2, ty2);
                fiz2 = simd_real_add(fiz2, tz2);
                fix3 = simd_real_add(fix3, tx3);
                fiy3 = simd_real_add(fiy3, ty3);
                fiz3 = simd_real_add(fiz3, tz3);

                if (cj < CJ1_FROM_CI(atom->Nclusters_local) || param->method) {
                    MD_SIMD_FLOAT tx_sum = simd_real_add(tx0,
                        simd_real_add(tx1, simd_real_add(tx2, tx3)));
                    MD_SIMD_FLOAT ty_sum = simd_real_add(ty0,
                        simd_real_add(ty1, simd_real_add(ty2, ty3)));
                    MD_SIMD_FLOAT tz_sum = simd_real_add(tz0,
                        simd_real_add(tz1, simd_real_add(tz2, tz3)));

                    simd_real_store(&cj_f[CL_X_INDEX_3D(0)],
                        simd_real_sub(simd_real_load(&cj_f[CL_X_INDEX_3D(0)]), tx_sum));
                    simd_real_store(&cj_f[CL_Y_INDEX_3D(0)],
                        simd_real_sub(simd_real_load(&cj_f[CL_Y_INDEX_3D(0)]), ty_sum));
                    simd_real_store(&cj_f[CL_Z_INDEX_3D(0)],
                        simd_real_sub(simd_real_load(&cj_f[CL_Z_INDEX_3D(0)]), tz_sum));
                }
            }

            for (int k = numneighs_masked; k < numneighs; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
                MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp = simd_real_load(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp = simd_real_load(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp = simd_real_load(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0  = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0  = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0  = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx1  = simd_real_sub(xi1_tmp, xj_tmp);
                MD_SIMD_FLOAT dely1  = simd_real_sub(yi1_tmp, yj_tmp);
                MD_SIMD_FLOAT delz1  = simd_real_sub(zi1_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2  = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2  = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2  = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT delx3  = simd_real_sub(xi3_tmp, xj_tmp);
                MD_SIMD_FLOAT dely3  = simd_real_sub(yi3_tmp, yj_tmp);
                MD_SIMD_FLOAT delz3  = simd_real_sub(zi3_tmp, zj_tmp);

                MD_SIMD_FLOAT rsq0 = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq1 = simd_real_fma(delx1,
                    delx1,
                    simd_real_fma(dely1, dely1, simd_real_mul(delz1, delz1)));
                MD_SIMD_FLOAT rsq2 = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));
                MD_SIMD_FLOAT rsq3 = simd_real_fma(delx3,
                    delx3,
                    simd_real_fma(dely3, dely3, simd_real_mul(delz3, delz3)));

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp = simd_i32_load(cj_t);
                MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec1  = simd_i32_add(tbase1, tj_tmp);
                MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);
                MD_SIMD_INT tvec3  = simd_i32_add(tbase3, tj_tmp);

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq1 = simd_real_gather(tvec1,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq3 = simd_real_gather(tvec3,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT sigma6_0 = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_1 = simd_real_gather(tvec1,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2 = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_3 = simd_real_gather(tvec3,
                    atom->sigma6,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT eps0 = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps1 = simd_real_gather(tvec1,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2 = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps3 = simd_real_gather(tvec3,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq1 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq3 = cutforcesq_vec;

                MD_SIMD_FLOAT sigma6_0 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_1 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_3 = sigma6_vec;

                MD_SIMD_FLOAT eps0      = eps_vec;
                MD_SIMD_FLOAT eps1      = eps_vec;
                MD_SIMD_FLOAT eps2      = eps_vec;
                MD_SIMD_FLOAT eps3      = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0, cutforcesq0);
                MD_SIMD_MASK cutoff_mask1 = simd_mask_cond_lt(rsq1, cutforcesq1);
                MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2, cutforcesq2);
                MD_SIMD_MASK cutoff_mask3 = simd_mask_cond_lt(rsq3, cutforcesq3);

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_1 = simd_real_reciprocal(rsq1);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);
                MD_SIMD_FLOAT sr2_3 = simd_real_reciprocal(rsq3);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_1 = simd_real_mul(sr2_1,
                    simd_real_mul(sr2_1, simd_real_mul(sr2_1, sigma6_1)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));
                MD_SIMD_FLOAT sr6_3 = simd_real_mul(sr2_3,
                    simd_real_mul(sr2_3, simd_real_mul(sr2_3, sigma6_3)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force1 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_1,
                        simd_real_mul(simd_real_sub(sr6_1, c05_vec),
                            simd_real_mul(sr2_1, eps1))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));
                MD_SIMD_FLOAT force3 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_3,
                        simd_real_mul(simd_real_sub(sr6_3, c05_vec),
                            simd_real_mul(sr2_3, eps3))));

                MD_SIMD_FLOAT tx0 = simd_real_select_by_mask(simd_real_mul(delx0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT ty0 = simd_real_select_by_mask(simd_real_mul(dely0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tz0 = simd_real_select_by_mask(simd_real_mul(delz0, force0),
                    cutoff_mask0);
                MD_SIMD_FLOAT tx1 = simd_real_select_by_mask(simd_real_mul(delx1, force1),
                    cutoff_mask1);
                MD_SIMD_FLOAT ty1 = simd_real_select_by_mask(simd_real_mul(dely1, force1),
                    cutoff_mask1);
                MD_SIMD_FLOAT tz1 = simd_real_select_by_mask(simd_real_mul(delz1, force1),
                    cutoff_mask1);
                MD_SIMD_FLOAT tx2 = simd_real_select_by_mask(simd_real_mul(delx2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT ty2 = simd_real_select_by_mask(simd_real_mul(dely2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT tz2 = simd_real_select_by_mask(simd_real_mul(delz2, force2),
                    cutoff_mask2);
                MD_SIMD_FLOAT tx3 = simd_real_select_by_mask(simd_real_mul(delx3, force3),
                    cutoff_mask3);
                MD_SIMD_FLOAT ty3 = simd_real_select_by_mask(simd_real_mul(dely3, force3),
                    cutoff_mask3);
                MD_SIMD_FLOAT tz3 = simd_real_select_by_mask(simd_real_mul(delz3, force3),
                    cutoff_mask3);

                fix0 = simd_real_add(fix0, tx0);
                fiy0 = simd_real_add(fiy0, ty0);
                fiz0 = simd_real_add(fiz0, tz0);
                fix1 = simd_real_add(fix1, tx1);
                fiy1 = simd_real_add(fiy1, ty1);
                fiz1 = simd_real_add(fiz1, tz1);
                fix2 = simd_real_add(fix2, tx2);
                fiy2 = simd_real_add(fiy2, ty2);
                fiz2 = simd_real_add(fiz2, tz2);
                fix3 = simd_real_add(fix3, tx3);
                fiy3 = simd_real_add(fiy3, ty3);
                fiz3 = simd_real_add(fiz3, tz3);

                if (cj < CJ1_FROM_CI(atom->Nclusters_local) || param->method) {
                    MD_SIMD_FLOAT tx_sum = simd_real_add(tx0,
                        simd_real_add(tx1, simd_real_add(tx2, tx3)));
                    MD_SIMD_FLOAT ty_sum = simd_real_add(ty0,
                        simd_real_add(ty1, simd_real_add(ty2, ty3)));
                    MD_SIMD_FLOAT tz_sum = simd_real_add(tz0,
                        simd_real_add(tz1, simd_real_add(tz2, tz3)));

                    simd_real_store(&cj_f[CL_X_INDEX_3D(0)],
                        simd_real_sub(simd_real_load(&cj_f[CL_X_INDEX_3D(0)]), tx_sum));
                    simd_real_store(&cj_f[CL_Y_INDEX_3D(0)],
                        simd_real_sub(simd_real_load(&cj_f[CL_Y_INDEX_3D(0)]), ty_sum));
                    simd_real_store(&cj_f[CL_Z_INDEX_3D(0)],
                        simd_real_sub(simd_real_load(&cj_f[CL_Z_INDEX_3D(0)]), tz_sum));
                }
            }

            simd_real_incr_reduced_sum(&ci_f[CL_X_INDEX_3D(0)], fix0, fix1, fix2, fix3);
            simd_real_incr_reduced_sum(&ci_f[CL_Y_INDEX_3D(0)], fiy0, fiy1, fiy2, fiy3);
            simd_real_incr_reduced_sum(&ci_f[CL_Z_INDEX_3D(0)], fiz0, fiz1, fiz2, fiz3);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters,
                (long long int)((double)numneighs * CLUSTER_M / CLUSTER_N));
        }
        if (param->method == eightShell) computeForceGhostShell(param, atom, neighbor);
        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_4xn end\n");
    return E - S;
}

double computeForceLJ4xnFullNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_4xn begin\n");
    const int nbM                = atom->Nclusters_local;
    const int nbN                = neighbor->maxneighs;
    MD_FLOAT cutforcesq          = param->cutforce * param->cutforce;
    MD_FLOAT sigma6              = param->sigma6;
    MD_FLOAT epsilon             = param->epsilon;
    MD_SIMD_FLOAT c48_vec        = simd_real_broadcast(48.0);
    MD_SIMD_FLOAT c05_vec        = simd_real_broadcast(0.5);

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_SIMD_FLOAT cutforcesq_vec = simd_real_broadcast(cutforcesq);
    MD_SIMD_FLOAT sigma6_vec     = simd_real_broadcast(sigma6);
    MD_SIMD_FLOAT eps_vec        = simd_real_broadcast(epsilon);
#endif

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ci_f  = &atom->cl_f[ci_vec_base];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ci_f[CL_X_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Y_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Z_INDEX_3D(cii)] = 0.0;
        }
    }

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_cj0           = CJ0_FROM_CI(ci);
#if CLUSTER_M > CLUSTER_N
            int ci_cj1           = CJ1_FROM_CI(ci);
#endif
            int ci_vec_base      = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x       = &atom->cl_x[ci_vec_base];
            MD_FLOAT* ci_f       = &atom->cl_f[ci_vec_base];
            int numneighs        = neighbor->numneigh[ci];
            int numneighs_masked = neighbor->numneigh_masked[ci];

            MD_SIMD_FLOAT xi0_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(0)]);
            MD_SIMD_FLOAT xi1_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(1)]);
            MD_SIMD_FLOAT xi2_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(2)]);
            MD_SIMD_FLOAT xi3_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(3)]);
            MD_SIMD_FLOAT yi0_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(0)]);
            MD_SIMD_FLOAT yi1_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(1)]);
            MD_SIMD_FLOAT yi2_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(2)]);
            MD_SIMD_FLOAT yi3_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(3)]);
            MD_SIMD_FLOAT zi0_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(0)]);
            MD_SIMD_FLOAT zi1_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(1)]);
            MD_SIMD_FLOAT zi2_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(2)]);
            MD_SIMD_FLOAT zi3_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(3)]);
            MD_SIMD_FLOAT fix0    = simd_real_zero();
            MD_SIMD_FLOAT fiy0    = simd_real_zero();
            MD_SIMD_FLOAT fiz0    = simd_real_zero();
            MD_SIMD_FLOAT fix1    = simd_real_zero();
            MD_SIMD_FLOAT fiy1    = simd_real_zero();
            MD_SIMD_FLOAT fiz1    = simd_real_zero();
            MD_SIMD_FLOAT fix2    = simd_real_zero();
            MD_SIMD_FLOAT fiy2    = simd_real_zero();
            MD_SIMD_FLOAT fiz2    = simd_real_zero();
            MD_SIMD_FLOAT fix3    = simd_real_zero();
            MD_SIMD_FLOAT fiy3    = simd_real_zero();
            MD_SIMD_FLOAT fiz3    = simd_real_zero();

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int ci_sca_base       = CI_SCALAR_BASE_INDEX(ci);
            int* ci_t             = &atom->cl_t[ci_sca_base];
            MD_SIMD_INT tbase0    = simd_i32_broadcast(ci_t[0] * atom->ntypes);
            MD_SIMD_INT tbase1    = simd_i32_broadcast(ci_t[1] * atom->ntypes);
            MD_SIMD_INT tbase2    = simd_i32_broadcast(ci_t[2] * atom->ntypes);
            MD_SIMD_INT tbase3    = simd_i32_broadcast(ci_t[3] * atom->ntypes);
#endif

            for (int k = 0; k < numneighs_masked; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp    = simd_real_load(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp    = simd_real_load(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp    = simd_real_load(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0     = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0     = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0     = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx1     = simd_real_sub(xi1_tmp, xj_tmp);
                MD_SIMD_FLOAT dely1     = simd_real_sub(yi1_tmp, yj_tmp);
                MD_SIMD_FLOAT delz1     = simd_real_sub(zi1_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2     = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2     = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2     = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT delx3     = simd_real_sub(xi3_tmp, xj_tmp);
                MD_SIMD_FLOAT dely3     = simd_real_sub(yi3_tmp, yj_tmp);
                MD_SIMD_FLOAT delz3     = simd_real_sub(zi3_tmp, zj_tmp);

#if CLUSTER_M == CLUSTER_N
                unsigned int cond0      = (unsigned int)(cj == ci_cj0);
                MD_SIMD_MASK excl_mask0 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 4 + 0]);
                MD_SIMD_MASK excl_mask1 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 4 + 1]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 4 + 2]);
                MD_SIMD_MASK excl_mask3 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 4 + 3]);
#else
#if CLUSTER_M < CLUSTER_N
                unsigned int cond0      = (unsigned int)((cj << 1) + 0 == ci);
                unsigned int cond1      = (unsigned int)((cj << 1) + 1 == ci);
#else
                unsigned int cond0 = (unsigned int)(cj == ci_cj0);
                unsigned int cond1 = (unsigned int)(cj == ci_cj1);
#endif
                MD_SIMD_MASK excl_mask0 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 0]);
                MD_SIMD_MASK excl_mask1 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 1]);
                MD_SIMD_MASK excl_mask2 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 2]);
                MD_SIMD_MASK excl_mask3 = simd_mask_from_u32(
                    atom->masks_4xn_fn[cond0 * 8 + cond1 * 4 + 3]);
#endif

                MD_SIMD_FLOAT rsq0 = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq1 = simd_real_fma(delx1,
                    delx1,
                    simd_real_fma(dely1, dely1, simd_real_mul(delz1, delz1)));
                MD_SIMD_FLOAT rsq2 = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));
                MD_SIMD_FLOAT rsq3 = simd_real_fma(delx3,
                    delx3,
                    simd_real_fma(dely3, dely3, simd_real_mul(delz3, delz3)));

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp = simd_i32_load(cj_t);
                MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec1  = simd_i32_add(tbase1, tj_tmp);
                MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);
                MD_SIMD_INT tvec3  = simd_i32_add(tbase3, tj_tmp);

                /*
                simd_print_real("tj_tmp", tj_tmp);
                simd_print_real("tvec0", tvec0);
                simd_print_real("tvec1", tvec1);
                simd_print_real("tvec2", tvec2);
                simd_print_real("tvec3", tvec3);
                */

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq1 = simd_real_gather(tvec1,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq3 = simd_real_gather(tvec3,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT sigma6_0 = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_1 = simd_real_gather(tvec1,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2 = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_3 = simd_real_gather(tvec3,
                    atom->sigma6,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT eps0 = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps1 = simd_real_gather(tvec1,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2 = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps3 = simd_real_gather(tvec3,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq1 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq3 = cutforcesq_vec;

                MD_SIMD_FLOAT sigma6_0 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_1 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_3 = sigma6_vec;

                MD_SIMD_FLOAT eps0        = eps_vec;
                MD_SIMD_FLOAT eps1        = eps_vec;
                MD_SIMD_FLOAT eps2        = eps_vec;
                MD_SIMD_FLOAT eps3        = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_and(excl_mask0,
                    simd_mask_cond_lt(rsq0, cutforcesq0));
                MD_SIMD_MASK cutoff_mask1 = simd_mask_and(excl_mask1,
                    simd_mask_cond_lt(rsq1, cutforcesq1));
                MD_SIMD_MASK cutoff_mask2 = simd_mask_and(excl_mask2,
                    simd_mask_cond_lt(rsq2, cutforcesq2));
                MD_SIMD_MASK cutoff_mask3 = simd_mask_and(excl_mask3,
                    simd_mask_cond_lt(rsq3, cutforcesq3));

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_1 = simd_real_reciprocal(rsq1);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);
                MD_SIMD_FLOAT sr2_3 = simd_real_reciprocal(rsq3);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_1 = simd_real_mul(sr2_1,
                    simd_real_mul(sr2_1, simd_real_mul(sr2_1, sigma6_1)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));
                MD_SIMD_FLOAT sr6_3 = simd_real_mul(sr2_3,
                    simd_real_mul(sr2_3, simd_real_mul(sr2_3, sigma6_3)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force1 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_1,
                        simd_real_mul(simd_real_sub(sr6_1, c05_vec),
                            simd_real_mul(sr2_1, eps1))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));
                MD_SIMD_FLOAT force3 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_3,
                        simd_real_mul(simd_real_sub(sr6_3, c05_vec),
                            simd_real_mul(sr2_3, eps3))));

                fix0 = simd_real_masked_add(fix0,
                    simd_real_mul(delx0, force0),
                    cutoff_mask0);
                fiy0 = simd_real_masked_add(fiy0,
                    simd_real_mul(dely0, force0),
                    cutoff_mask0);
                fiz0 = simd_real_masked_add(fiz0,
                    simd_real_mul(delz0, force0),
                    cutoff_mask0);
                fix1 = simd_real_masked_add(fix1,
                    simd_real_mul(delx1, force1),
                    cutoff_mask1);
                fiy1 = simd_real_masked_add(fiy1,
                    simd_real_mul(dely1, force1),
                    cutoff_mask1);
                fiz1 = simd_real_masked_add(fiz1,
                    simd_real_mul(delz1, force1),
                    cutoff_mask1);
                fix2 = simd_real_masked_add(fix2,
                    simd_real_mul(delx2, force2),
                    cutoff_mask2);
                fiy2 = simd_real_masked_add(fiy2,
                    simd_real_mul(dely2, force2),
                    cutoff_mask2);
                fiz2 = simd_real_masked_add(fiz2,
                    simd_real_mul(delz2, force2),
                    cutoff_mask2);
                fix3 = simd_real_masked_add(fix3,
                    simd_real_mul(delx3, force3),
                    cutoff_mask3);
                fiy3 = simd_real_masked_add(fiy3,
                    simd_real_mul(dely3, force3),
                    cutoff_mask3);
                fiz3 = simd_real_masked_add(fiz3,
                    simd_real_mul(delz3, force3),
                    cutoff_mask3);
            }

            for (int k = numneighs_masked; k < numneighs; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, nbN);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                MD_SIMD_FLOAT xj_tmp = simd_real_load(&cj_x[CL_X_INDEX_3D(0)]);
                MD_SIMD_FLOAT yj_tmp = simd_real_load(&cj_x[CL_Y_INDEX_3D(0)]);
                MD_SIMD_FLOAT zj_tmp = simd_real_load(&cj_x[CL_Z_INDEX_3D(0)]);
                MD_SIMD_FLOAT delx0  = simd_real_sub(xi0_tmp, xj_tmp);
                MD_SIMD_FLOAT dely0  = simd_real_sub(yi0_tmp, yj_tmp);
                MD_SIMD_FLOAT delz0  = simd_real_sub(zi0_tmp, zj_tmp);
                MD_SIMD_FLOAT delx1  = simd_real_sub(xi1_tmp, xj_tmp);
                MD_SIMD_FLOAT dely1  = simd_real_sub(yi1_tmp, yj_tmp);
                MD_SIMD_FLOAT delz1  = simd_real_sub(zi1_tmp, zj_tmp);
                MD_SIMD_FLOAT delx2  = simd_real_sub(xi2_tmp, xj_tmp);
                MD_SIMD_FLOAT dely2  = simd_real_sub(yi2_tmp, yj_tmp);
                MD_SIMD_FLOAT delz2  = simd_real_sub(zi2_tmp, zj_tmp);
                MD_SIMD_FLOAT delx3  = simd_real_sub(xi3_tmp, xj_tmp);
                MD_SIMD_FLOAT dely3  = simd_real_sub(yi3_tmp, yj_tmp);
                MD_SIMD_FLOAT delz3  = simd_real_sub(zi3_tmp, zj_tmp);

                MD_SIMD_FLOAT rsq0 = simd_real_fma(delx0,
                    delx0,
                    simd_real_fma(dely0, dely0, simd_real_mul(delz0, delz0)));
                MD_SIMD_FLOAT rsq1 = simd_real_fma(delx1,
                    delx1,
                    simd_real_fma(dely1, dely1, simd_real_mul(delz1, delz1)));
                MD_SIMD_FLOAT rsq2 = simd_real_fma(delx2,
                    delx2,
                    simd_real_fma(dely2, dely2, simd_real_mul(delz2, delz2)));
                MD_SIMD_FLOAT rsq3 = simd_real_fma(delx3,
                    delx3,
                    simd_real_fma(dely3, dely3, simd_real_mul(delz3, delz3)));

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                MD_SIMD_INT tj_tmp = simd_i32_load(cj_t);
                MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                MD_SIMD_INT tvec1  = simd_i32_add(tbase1, tj_tmp);
                MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);
                MD_SIMD_INT tvec3  = simd_i32_add(tbase3, tj_tmp);

                MD_SIMD_FLOAT cutforcesq0 = simd_real_gather(tvec0,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq1 = simd_real_gather(tvec1,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq2 = simd_real_gather(tvec2,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT cutforcesq3 = simd_real_gather(tvec3,
                    atom->cutforcesq,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT sigma6_0 = simd_real_gather(tvec0,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_1 = simd_real_gather(tvec1,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_2 = simd_real_gather(tvec2,
                    atom->sigma6,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT sigma6_3 = simd_real_gather(tvec3,
                    atom->sigma6,
                    sizeof(MD_FLOAT));

                MD_SIMD_FLOAT eps0 = simd_real_gather(tvec0,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps1 = simd_real_gather(tvec1,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps2 = simd_real_gather(tvec2,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
                MD_SIMD_FLOAT eps3 = simd_real_gather(tvec3,
                    atom->epsilon,
                    sizeof(MD_FLOAT));
#else
                MD_SIMD_FLOAT cutforcesq0 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq1 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq2 = cutforcesq_vec;
                MD_SIMD_FLOAT cutforcesq3 = cutforcesq_vec;

                MD_SIMD_FLOAT sigma6_0 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_1 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_2 = sigma6_vec;
                MD_SIMD_FLOAT sigma6_3 = sigma6_vec;

                MD_SIMD_FLOAT eps0 = eps_vec;
                MD_SIMD_FLOAT eps1 = eps_vec;
                MD_SIMD_FLOAT eps2 = eps_vec;
                MD_SIMD_FLOAT eps3 = eps_vec;
#endif

                MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0, cutforcesq0);
                MD_SIMD_MASK cutoff_mask1 = simd_mask_cond_lt(rsq1, cutforcesq1);
                MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2, cutforcesq2);
                MD_SIMD_MASK cutoff_mask3 = simd_mask_cond_lt(rsq3, cutforcesq3);

                MD_SIMD_FLOAT sr2_0 = simd_real_reciprocal(rsq0);
                MD_SIMD_FLOAT sr2_1 = simd_real_reciprocal(rsq1);
                MD_SIMD_FLOAT sr2_2 = simd_real_reciprocal(rsq2);
                MD_SIMD_FLOAT sr2_3 = simd_real_reciprocal(rsq3);

                MD_SIMD_FLOAT sr6_0 = simd_real_mul(sr2_0,
                    simd_real_mul(sr2_0, simd_real_mul(sr2_0, sigma6_0)));
                MD_SIMD_FLOAT sr6_1 = simd_real_mul(sr2_1,
                    simd_real_mul(sr2_1, simd_real_mul(sr2_1, sigma6_1)));
                MD_SIMD_FLOAT sr6_2 = simd_real_mul(sr2_2,
                    simd_real_mul(sr2_2, simd_real_mul(sr2_2, sigma6_2)));
                MD_SIMD_FLOAT sr6_3 = simd_real_mul(sr2_3,
                    simd_real_mul(sr2_3, simd_real_mul(sr2_3, sigma6_3)));

                MD_SIMD_FLOAT force0 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_0,
                        simd_real_mul(simd_real_sub(sr6_0, c05_vec),
                            simd_real_mul(sr2_0, eps0))));
                MD_SIMD_FLOAT force1 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_1,
                        simd_real_mul(simd_real_sub(sr6_1, c05_vec),
                            simd_real_mul(sr2_1, eps1))));
                MD_SIMD_FLOAT force2 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_2,
                        simd_real_mul(simd_real_sub(sr6_2, c05_vec),
                            simd_real_mul(sr2_2, eps2))));
                MD_SIMD_FLOAT force3 = simd_real_mul(c48_vec,
                    simd_real_mul(sr6_3,
                        simd_real_mul(simd_real_sub(sr6_3, c05_vec),
                            simd_real_mul(sr2_3, eps3))));

                fix0 = simd_real_masked_add(fix0,
                    simd_real_mul(delx0, force0),
                    cutoff_mask0);
                fiy0 = simd_real_masked_add(fiy0,
                    simd_real_mul(dely0, force0),
                    cutoff_mask0);
                fiz0 = simd_real_masked_add(fiz0,
                    simd_real_mul(delz0, force0),
                    cutoff_mask0);
                fix1 = simd_real_masked_add(fix1,
                    simd_real_mul(delx1, force1),
                    cutoff_mask1);
                fiy1 = simd_real_masked_add(fiy1,
                    simd_real_mul(dely1, force1),
                    cutoff_mask1);
                fiz1 = simd_real_masked_add(fiz1,
                    simd_real_mul(delz1, force1),
                    cutoff_mask1);
                fix2 = simd_real_masked_add(fix2,
                    simd_real_mul(delx2, force2),
                    cutoff_mask2);
                fiy2 = simd_real_masked_add(fiy2,
                    simd_real_mul(dely2, force2),
                    cutoff_mask2);
                fiz2 = simd_real_masked_add(fiz2,
                    simd_real_mul(delz2, force2),
                    cutoff_mask2);
                fix3 = simd_real_masked_add(fix3,
                    simd_real_mul(delx3, force3),
                    cutoff_mask3);
                fiy3 = simd_real_masked_add(fiy3,
                    simd_real_mul(dely3, force3),
                    cutoff_mask3);
                fiz3 = simd_real_masked_add(fiz3,
                    simd_real_mul(delz3, force3),
                    cutoff_mask3);
            }

            simd_real_incr_reduced_sum(&ci_f[CL_X_INDEX_3D(0)], fix0, fix1, fix2, fix3);
            simd_real_incr_reduced_sum(&ci_f[CL_Y_INDEX_3D(0)], fiy0, fiy1, fiy2, fiy3);
            simd_real_incr_reduced_sum(&ci_f[CL_Z_INDEX_3D(0)], fiz0, fiz1, fiz2, fiz3);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters, (long long int)((double)numneighs));
            // addStat(stats->force_iters, (long long int)((double)numneighs * CLUSTER_M /
            // CLUSTER_N));
        }

        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_4xn end\n");
    return E - S;
}
#endif

// Routine for eight shell method + MPI
void computeForceGhostShell(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("computeForceLJ begin\n");

    int Nshell = neighbor->Nshell;
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;

    for (int ci = 0; ci < Nshell; ci++) {
        int numneighs   = neighbor->numNeighShell[ci];
        int cs          = neighbor->listshell[ci];
        int cs_vec_base = CJ_VECTOR3_BASE_INDEX(cs);
        MD_FLOAT* cs_x  = &atom->cl_x[cs_vec_base];
        MD_FLOAT* cs_f  = &atom->cl_f[cs_vec_base];

        for (int k = 0; k < numneighs; k++) {
            const int cj    = neighs(neighbor->neighshell, ci, k, Nshell, neighbor->maxneighs);
            int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
            MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
            MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];

            for (int css = 0; css < CLUSTER_N; css++) {
                MD_FLOAT x = cs_x[CL_X_INDEX_3D(css)];
                MD_FLOAT y = cs_x[CL_Y_INDEX_3D(css)];
                MD_FLOAT z = cs_x[CL_Z_INDEX_3D(css)];

                MD_FLOAT fix = 0;
                MD_FLOAT fiy = 0;
                MD_FLOAT fiz = 0;

                for (int cjj = 0; cjj < CLUSTER_N; cjj++) {

                    MD_FLOAT delx = x - cj_x[CL_X_INDEX_3D(cjj)];
                    MD_FLOAT dely = y - cj_x[CL_Y_INDEX_3D(cjj)];
                    MD_FLOAT delz = z - cj_x[CL_Z_INDEX_3D(cjj)];

                    MD_FLOAT rsq = delx * delx + dely * dely + delz * delz;
                    if (rsq < cutforcesq) {
                        MD_FLOAT sr2   = 1.0 / rsq;
                        MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
                        MD_FLOAT force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;

                        cj_f[CL_X_INDEX_3D(cjj)] -= delx * force;
                        cj_f[CL_Y_INDEX_3D(cjj)] -= dely * force;
                        cj_f[CL_Z_INDEX_3D(cjj)] -= delz * force;

                        fix += delx * force;
                        fiy += dely * force;
                        fiz += delz * force;
                    }
                }

                cs_f[CL_X_INDEX_3D(css)] += fix;
                cs_f[CL_Y_INDEX_3D(css)] += fiy;
                cs_f[CL_Z_INDEX_3D(css)] += fiz;
            }
        }
        // addStat(stats->calculated_forces, 1);
        // addStat(stats->num_neighs, numneighs);
        // addStat(stats->force_iters, (long long int)((double)numneighs));
    }
}
