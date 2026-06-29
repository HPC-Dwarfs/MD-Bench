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
#include <ljtable.h>
#include <math.h>
#include <neighbor.h>
#include <parameter.h>
#include <simd.h>
#include <stats.h>
#include <timing.h>
#include <util.h>

void computeForceGhostShell(Parameter*, Atom*, Neighbor*);

#ifdef USE_REFERENCE_KERNEL
double computeForceLJTableRef(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJTable begin\n");
    const int nbM = atom->Nclusters_local;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif
    const int nm1         = ljtable.n - 1;
    const MD_FLOAT inv_h  = ljtable.inv_h;
    const MD_FLOAT* coeff = ljtable.coeff;

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
            int numneighs   = neighbor->numneigh_inner[ci];

#if LJ_COMB_RULE == LJ_COMB_GEOM
            int ci_sca_base       = CI_SCALAR_BASE_INDEX(ci);
            MD_FLOAT* ci_sqrt_eps = &atom->cl_sqrt_epsilon[ci_sca_base];
            MD_FLOAT* ci_sigma3   = &atom->cl_sigma3[ci_sca_base];
#elif LJ_COMB_RULE == LJ_COMB_NONE
            int ci_sca_base = CI_SCALAR_BASE_INDEX(ci);
            int* ci_t       = &atom->cl_t[ci_sca_base];
#endif

            for (int k = 0; k < numneighs; k++) {
                const int cj    = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
                int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                int any         = 0;
                MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
                MD_FLOAT* cj_f  = &atom->cl_f[cj_vec_base];

#if LJ_COMB_RULE == LJ_COMB_GEOM
                int cj_sca_base       = CJ_SCALAR_BASE_INDEX(cj);
                MD_FLOAT* cj_sqrt_eps = &atom->cl_sqrt_epsilon[cj_sca_base];
                MD_FLOAT* cj_sigma3   = &atom->cl_sigma3[cj_sca_base];
#elif LJ_COMB_RULE == LJ_COMB_NONE
                int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

                for (int cii = 0; cii < CLUSTER_M; cii++) {
#if LJ_COMB_RULE == LJ_COMB_NONE
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

#if LJ_COMB_RULE == LJ_COMB_GEOM
                            MD_FLOAT sigma6  = ci_sigma3[cii] * cj_sigma3[cjj];
                            MD_FLOAT epsilon = ci_sqrt_eps[cii] * cj_sqrt_eps[cjj];
#elif LJ_COMB_RULE == LJ_COMB_NONE
                            int type_j          = cj_t[cjj];
                            int type_index      = type_i * atom->ntypes + type_j;
                            MD_FLOAT cutforcesq = atom->cutforcesq[type_index];
                            MD_FLOAT sigma6     = atom->sigma6[type_index];
                            MD_FLOAT epsilon    = atom->epsilon[type_index];
#endif

                            if (rsq < cutforcesq) {
                                MD_FLOAT u = LJ_TABLE_COORD(rsq) * inv_h;
                                int m      = (int)u;
                                if (m > nm1) m = nm1;
                                MD_FLOAT frac      = u - (MD_FLOAT)m;
                                const MD_FLOAT* cc = &coeff[m * LJ_TABLE_STRIDE];
                                MD_FLOAT hrep      = cc[0] +
                                                frac * (cc[1] +
                                                           frac * (cc[2] + frac * cc[3]));
                                MD_FLOAT gdisp = cc[4] +
                                                 frac *
                                                     (cc[5] +
                                                         frac * (cc[6] + frac * cc[7]));
                                MD_FLOAT force = epsilon * sigma6 *
                                                 (sigma6 * hrep + gdisp);

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
    DEBUG_MESSAGE("computeForceLJTable end\n");
    return E - S;
}
#endif
