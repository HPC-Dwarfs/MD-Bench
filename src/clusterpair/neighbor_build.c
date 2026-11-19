/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// Cluster-pair neighbor-list build, a template on LJ_COMB_RULE: the Makefile compiles it
// once per rule and buildNeighborCPU() in neighbor.c calls the variant of the rule
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

#include <allocate.h>
#include <atom.h>
#include <device.h>
#include <force.h>
#include <kernel_variant.h>
#include <neighbor.h>
#include <neighbor_internal.h>
#include <parameter.h>
#include <simd.h>
#include <util.h>

/* Returns a diagonal or off-diagonal interaction mask for plain C lists */
static unsigned int get_imask(int rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=2 */
static unsigned int get_imask_simd_j2(int rdiag, int ci, int cj)
{
    return (rdiag && ci * 2 == cj
                ? NBNXN_INTERACTION_MASK_DIAG_J2_0
                : (rdiag && ci * 2 + 1 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_1
                                             : NBNXN_INTERACTION_MASK_ALL));
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=4 */
static unsigned int get_imask_simd_j4(int rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=8 */
static unsigned int get_imask_simd_j8(int rdiag, int ci, int cj)
{
    return (rdiag && ci == cj * 2
                ? NBNXN_INTERACTION_MASK_DIAG_J8_0
                : (rdiag && ci == cj * 2 + 1 ? NBNXN_INTERACTION_MASK_DIAG_J8_1
                                             : NBNXN_INTERACTION_MASK_ALL));
}

// The mask only depends on the ratio of the cluster sizes
#if CLUSTER_M == CLUSTER_N
#define get_imask_simd get_imask_simd_j4
#elif CLUSTER_M == CLUSTER_N * 2
#define get_imask_simd get_imask_simd_j2
#elif CLUSTER_N == CLUSTER_M * 2
#define get_imask_simd get_imask_simd_j8
#else
#error "Invalid cluster configuration"
#endif

/* Returns whether any atom of i-cluster ci is excluded from any atom of j-cluster cj,
 * i.e. both belong to the same molecule. Such cluster pairs are kept in the masked
 * part of the neighbor list, where the kernels evaluate the exclusions (only for systems
 * with molecules, atom->molsize > 1). */
static int clusterPairHasExclusion(Atom* atom, int ci, int cj)
{
    const int* mol_i = &atom->cl_mol[CI_SCALAR_BASE_INDEX(ci)];
    const int* mol_j = &atom->cl_mol[CJ_SCALAR_BASE_INDEX(cj)];

    for (int cii = 0; cii < CLUSTER_M; cii++) {
        for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
            if (mol_i[cii] >= 0 && mol_i[cii] == mol_j[cjj]) {
                return 1;
            }
        }
    }

    return 0;
}

void RULE_VARIANT_NAME(buildNeighborCPU)(Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("buildNeighbor start\n");

    /* extend atom arrays if necessary */
    if (atom->Nclusters_local > nmax) {
        nmax = atom->Nclusters_local;
        if (neighbor->numneigh) free(neighbor->numneigh);
        if (neighbor->numneigh_masked) free(neighbor->numneigh_masked);
        if (neighbor->numneigh_inner) free(neighbor->numneigh_inner);
        if (neighbor->numneigh_inner_masked) free(neighbor->numneigh_inner_masked);
        if (neighbor->neighbors) free(neighbor->neighbors);
        if (neighbor->neighbors_imask) free(neighbor->neighbors_imask);
        neighbor->numneigh              = (int*)allocate(ALIGNMENT, nmax * sizeof(int));
        neighbor->numneigh_masked       = (int*)allocate(ALIGNMENT, nmax * sizeof(int));
        neighbor->numneigh_inner        = (int*)allocate(ALIGNMENT, nmax * sizeof(int));
        neighbor->numneigh_inner_masked = (int*)allocate(ALIGNMENT, nmax * sizeof(int));
        neighbor->neighbors             = (int*)allocate(ALIGNMENT,
            nmax * neighbor->maxneighs * sizeof(int));
        neighbor->neighbors_imask       = (unsigned int*)allocate(ALIGNMENT,
            nmax * neighbor->maxneighs * sizeof(unsigned int));

#ifdef NBLIST_CSR
        if (neighbor->neigh_start) free(neighbor->neigh_start);
        neighbor->neigh_start = (int*)allocate(ALIGNMENT, (nmax + 1) * sizeof(int));
    } else {
        free(neighbor->neighbors);
        free(neighbor->neighbors_imask);
        neighbor->neighbors       = (int*)allocate(ALIGNMENT,
            nmax * neighbor->maxneighs * sizeof(int));
        neighbor->neighbors_imask = (unsigned int*)allocate(ALIGNMENT,
            nmax * neighbor->maxneighs * sizeof(unsigned int));
#endif
    }

    MD_FLOAT bbx    = 0.5 * (binsizex + binsizex);
    MD_FLOAT bby    = 0.5 * (binsizey + binsizey);
    MD_FLOAT rbb_sq = MAX(0.0, cutneigh - 0.5 * sqrt(bbx * bbx + bby * bby));
    rbb_sq          = rbb_sq * rbb_sq;
    int resize      = 1;

#ifdef CUDA_TARGET
    /* Ghost cj must be resolved to its real mirror before ordering, else a
     * periodic-boundary pair gets counted from both sides. */
    const int ncj_for_order = get_ncj_from_nci(atom->Nclusters_local);
#endif

    /* loop over each atom, storing neighbors */
    while (resize) {
        const int nbM     = atom->Nclusters_local;
        const int nbN     = neighbor->maxneighs;
        int new_maxneighs = neighbor->maxneighs;
        int resize_local  = 0;

#pragma omp parallel for schedule(runtime) reduction(max : new_maxneighs)                \
    reduction(| : resize_local)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_cj0 = CJ0_FROM_CI(ci);
            int n = 0, nmasked = 0;
            int ibin        = atom->cluster_bin[ci];
            int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x  = &atom->cl_x[ci_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int ci_sca_base = CI_SCALAR_BASE_INDEX(ci);
            int* ci_t       = &atom->cl_t[ci_sca_base];
#endif

            MD_FLOAT ibb_xmin = atom->iclusters[ci].bbminx;
            MD_FLOAT ibb_xmax = atom->iclusters[ci].bbmaxx;
            MD_FLOAT ibb_ymin = atom->iclusters[ci].bbminy;
            MD_FLOAT ibb_ymax = atom->iclusters[ci].bbmaxy;
            MD_FLOAT ibb_zmin = atom->iclusters[ci].bbminz;
            MD_FLOAT ibb_zmax = atom->iclusters[ci].bbmaxz;

#if defined(CLUSTERPAIR_KERNEL_2XNN)
            MD_SIMD_FLOAT xi0_tmp = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(0)]);
            MD_SIMD_FLOAT xi2_tmp = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(2)]);
            MD_SIMD_FLOAT yi0_tmp = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(0)]);
            MD_SIMD_FLOAT yi2_tmp = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(2)]);
            MD_SIMD_FLOAT zi0_tmp = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(0)]);
            MD_SIMD_FLOAT zi2_tmp = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(2)]);

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            MD_SIMD_INT tbase0 = simd_i32_load_h_dual_scaled(&ci_t[0], atom->ntypes);
            MD_SIMD_INT tbase2 = simd_i32_load_h_dual_scaled(&ci_t[2], atom->ntypes);
#else
            MD_SIMD_FLOAT cutneighsq_vec = simd_real_broadcast(cutneighsq);
#endif

#elif defined(CLUSTERPAIR_KERNEL_4XN)
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

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            MD_SIMD_INT tbase0 = simd_i32_broadcast(ci_t[0] * atom->ntypes);
            MD_SIMD_INT tbase1 = simd_i32_broadcast(ci_t[1] * atom->ntypes);
            MD_SIMD_INT tbase2 = simd_i32_broadcast(ci_t[2] * atom->ntypes);
            MD_SIMD_INT tbase3 = simd_i32_broadcast(ci_t[3] * atom->ntypes);
#else
            MD_SIMD_FLOAT cutneighsq_vec = simd_real_broadcast(cutneighsq);
#endif

#elif defined(CLUSTERPAIR_KERNEL_2XN)
            MD_SIMD_FLOAT xi0_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(0)]);
            MD_SIMD_FLOAT xi1_tmp = simd_real_broadcast(ci_x[CL_X_INDEX_3D(1)]);
            MD_SIMD_FLOAT yi0_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(0)]);
            MD_SIMD_FLOAT yi1_tmp = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(1)]);
            MD_SIMD_FLOAT zi0_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(0)]);
            MD_SIMD_FLOAT zi1_tmp = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(1)]);

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            MD_SIMD_INT tbase0 = simd_i32_broadcast(ci_t[0] * atom->ntypes);
            MD_SIMD_INT tbase1 = simd_i32_broadcast(ci_t[1] * atom->ntypes);
#else
            MD_SIMD_FLOAT cutneighsq_vec = simd_real_broadcast(cutneighsq);

#endif

#endif
            for (int k = 0; k < nstencil; k++) {
                int jbin     = ibin + stencil[k];
                int* loc_bin = &bin_clusters[jbin * clusters_per_bin];
                int cj, m = -1;
                MD_FLOAT jbb_xmin, jbb_xmax, jbb_ymin, jbb_ymax, jbb_zmin, jbb_zmax;
                const int c = bin_nclusters[jbin];

                if (c > 0) {
                    // d_bb_sq starts at 0: when the first candidate is skipped by the
                    // half-list order check, the front skip stops there
                    MD_FLOAT dl, dh, dm, dm0, d_bb_sq = 0.0;

                    do {
                        m++;
                        cj = loc_bin[m];
#ifdef CUDA_TARGET
                        if (neighbor->half_neigh &&
                            ci_cj0 > (cj < ncj_for_order
                                             ? cj
                                             : atom->border_map[cj - ncj_for_order])) {
                            continue;
                        }
#else
                        if (neighbor->half_neigh && ci_cj0 > cj) {
                            continue;
                        }
#endif
                        jbb_zmin = atom->jclusters[cj].bbminz;
                        jbb_zmax = atom->jclusters[cj].bbmaxz;
                        dl       = ibb_zmin - jbb_zmax;
                        dh       = jbb_zmin - ibb_zmax;
                        dm       = MAX(dl, dh);
                        dm0      = MAX(dm, 0.0);
                        d_bb_sq  = dm0 * dm0;
                    } while (m + 1 < c && d_bb_sq > cutneighsq);

                    jbb_xmin = atom->jclusters[cj].bbminx;
                    jbb_xmax = atom->jclusters[cj].bbmaxx;
                    jbb_ymin = atom->jclusters[cj].bbminy;
                    jbb_ymax = atom->jclusters[cj].bbmaxy;

                    while (m < c) {
                        // Clusters in a bin are sorted by bbminz, so once one starts
                        // beyond the cut-off above the i-cluster, all following do too
                        // (jbb_zmin is not loaded for clusters skipped by the half-list
                        // order check above, so read it from the cluster)
                        dh = atom->jclusters[cj].bbminz - ibb_zmax;
                        if (dh > 0.0 && dh * dh > cutneighsq) {
                            break;
                        }

#ifdef CUDA_TARGET
                        if (!neighbor->half_neigh ||
                            ci_cj0 <= (cj < ncj_for_order
                                              ? cj
                                              : atom->border_map[cj - ncj_for_order])) {
#else
                        if (!neighbor->half_neigh || ci_cj0 <= cj) {
#endif
                            dl      = ibb_zmin - jbb_zmax;
                            dh      = jbb_zmin - ibb_zmax;
                            dm      = MAX(dl, dh);
                            dm0     = MAX(dm, 0.0);
                            d_bb_sq = dm0 * dm0;

                            /*if(d_bb_sq > cutneighsq) {
                                break;
                            }*/

                            dl  = ibb_ymin - jbb_ymax;
                            dh  = jbb_ymin - ibb_ymax;
                            dm  = MAX(dl, dh);
                            dm0 = MAX(dm, 0.0);
                            d_bb_sq += dm0 * dm0;

                            dl  = ibb_xmin - jbb_xmax;
                            dh  = jbb_xmin - ibb_xmax;
                            dm  = MAX(dl, dh);
                            dm0 = MAX(dm, 0.0);
                            d_bb_sq += dm0 * dm0;

                            if (d_bb_sq < cutneighsq) {
                                int is_neighbor = (d_bb_sq < rbb_sq) ? 1 : 0;
                                if (!is_neighbor) {
                                    int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                                    MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                                    int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);
                                    int* cj_t       = &atom->cl_t[cj_sca_base];
#endif

#if defined(CLUSTERPAIR_KERNEL_2XNN)

                                    MD_SIMD_FLOAT xj_tmp = simd_real_load_h_duplicate(
                                        &cj_x[CL_X_INDEX_3D(0)]);
                                    MD_SIMD_FLOAT yj_tmp = simd_real_load_h_duplicate(
                                        &cj_x[CL_Y_INDEX_3D(0)]);
                                    MD_SIMD_FLOAT zj_tmp = simd_real_load_h_duplicate(
                                        &cj_x[CL_Z_INDEX_3D(0)]);

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                                    MD_SIMD_INT tj_tmp = simd_i32_load_h_duplicate(cj_t);
                                    MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                                    MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);

                                    MD_SIMD_FLOAT cutneighsq0 = simd_real_gather(tvec0,
                                        atom->cutneighsq,
                                        sizeof(MD_FLOAT));
                                    MD_SIMD_FLOAT cutneighsq2 = simd_real_gather(tvec2,
                                        atom->cutneighsq,
                                        sizeof(MD_FLOAT));
#else
                                    MD_SIMD_FLOAT cutneighsq0 = cutneighsq_vec;
                                    MD_SIMD_FLOAT cutneighsq2 = cutneighsq_vec;
#endif

                                    MD_SIMD_FLOAT delx0 = simd_real_sub(xi0_tmp, xj_tmp);
                                    MD_SIMD_FLOAT dely0 = simd_real_sub(yi0_tmp, yj_tmp);
                                    MD_SIMD_FLOAT delz0 = simd_real_sub(zi0_tmp, zj_tmp);
                                    MD_SIMD_FLOAT delx2 = simd_real_sub(xi2_tmp, xj_tmp);
                                    MD_SIMD_FLOAT dely2 = simd_real_sub(yi2_tmp, yj_tmp);
                                    MD_SIMD_FLOAT delz2 = simd_real_sub(zi2_tmp, zj_tmp);
                                    MD_SIMD_FLOAT rsq0  = simd_real_fma(delx0,
                                        delx0,
                                        simd_real_fma(dely0,
                                            dely0,
                                            simd_real_mul(delz0, delz0)));
                                    MD_SIMD_FLOAT rsq2  = simd_real_fma(delx2,
                                        delx2,
                                        simd_real_fma(dely2,
                                            dely2,
                                            simd_real_mul(delz2, delz2)));

                                    MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0,
                                        cutneighsq0);
                                    MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2,
                                        cutneighsq2);

                                    if (simd_test_any(cutoff_mask0) ||
                                        simd_test_any(cutoff_mask2)) {
                                        is_neighbor = 1;
                                    }

#elif defined(CLUSTERPAIR_KERNEL_4XN)

                                    MD_SIMD_FLOAT xj_tmp = simd_real_load(
                                        &cj_x[CL_X_INDEX_3D(0)]);
                                    MD_SIMD_FLOAT yj_tmp = simd_real_load(
                                        &cj_x[CL_Y_INDEX_3D(0)]);
                                    MD_SIMD_FLOAT zj_tmp = simd_real_load(
                                        &cj_x[CL_Z_INDEX_3D(0)]);
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                                    MD_SIMD_INT tj_tmp = simd_i32_load(cj_t);
                                    MD_SIMD_INT tvec0  = simd_i32_add(tbase0, tj_tmp);
                                    MD_SIMD_INT tvec1  = simd_i32_add(tbase1, tj_tmp);
                                    MD_SIMD_INT tvec2  = simd_i32_add(tbase2, tj_tmp);
                                    MD_SIMD_INT tvec3  = simd_i32_add(tbase3, tj_tmp);

                                    MD_SIMD_FLOAT cutneighsq0 = simd_real_gather(tvec0,
                                        atom->cutneighsq,
                                        sizeof(MD_FLOAT));
                                    MD_SIMD_FLOAT cutneighsq1 = simd_real_gather(tvec1,
                                        atom->cutneighsq,
                                        sizeof(MD_FLOAT));
                                    MD_SIMD_FLOAT cutneighsq2 = simd_real_gather(tvec2,
                                        atom->cutneighsq,
                                        sizeof(MD_FLOAT));
                                    MD_SIMD_FLOAT cutneighsq3 = simd_real_gather(tvec3,
                                        atom->cutneighsq,
                                        sizeof(MD_FLOAT));
#else
                                    MD_SIMD_FLOAT cutneighsq0 = cutneighsq_vec;
                                    MD_SIMD_FLOAT cutneighsq1 = cutneighsq_vec;
                                    MD_SIMD_FLOAT cutneighsq2 = cutneighsq_vec;
                                    MD_SIMD_FLOAT cutneighsq3 = cutneighsq_vec;
#endif

                                    MD_SIMD_FLOAT delx0 = simd_real_sub(xi0_tmp, xj_tmp);
                                    MD_SIMD_FLOAT dely0 = simd_real_sub(yi0_tmp, yj_tmp);
                                    MD_SIMD_FLOAT delz0 = simd_real_sub(zi0_tmp, zj_tmp);
                                    MD_SIMD_FLOAT delx1 = simd_real_sub(xi1_tmp, xj_tmp);
                                    MD_SIMD_FLOAT dely1 = simd_real_sub(yi1_tmp, yj_tmp);
                                    MD_SIMD_FLOAT delz1 = simd_real_sub(zi1_tmp, zj_tmp);
                                    MD_SIMD_FLOAT delx2 = simd_real_sub(xi2_tmp, xj_tmp);
                                    MD_SIMD_FLOAT dely2 = simd_real_sub(yi2_tmp, yj_tmp);
                                    MD_SIMD_FLOAT delz2 = simd_real_sub(zi2_tmp, zj_tmp);
                                    MD_SIMD_FLOAT delx3 = simd_real_sub(xi3_tmp, xj_tmp);
                                    MD_SIMD_FLOAT dely3 = simd_real_sub(yi3_tmp, yj_tmp);
                                    MD_SIMD_FLOAT delz3 = simd_real_sub(zi3_tmp, zj_tmp);

                                    MD_SIMD_FLOAT rsq0 = simd_real_fma(delx0,
                                        delx0,
                                        simd_real_fma(dely0,
                                            dely0,
                                            simd_real_mul(delz0, delz0)));
                                    MD_SIMD_FLOAT rsq1 = simd_real_fma(delx1,
                                        delx1,
                                        simd_real_fma(dely1,
                                            dely1,
                                            simd_real_mul(delz1, delz1)));
                                    MD_SIMD_FLOAT rsq2 = simd_real_fma(delx2,
                                        delx2,
                                        simd_real_fma(dely2,
                                            dely2,
                                            simd_real_mul(delz2, delz2)));
                                    MD_SIMD_FLOAT rsq3 = simd_real_fma(delx3,
                                        delx3,
                                        simd_real_fma(dely3,
                                            dely3,
                                            simd_real_mul(delz3, delz3)));

                                    MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0,
                                        cutneighsq0);
                                    MD_SIMD_MASK cutoff_mask1 = simd_mask_cond_lt(rsq1,
                                        cutneighsq1);
                                    MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2,
                                        cutneighsq2);
                                    MD_SIMD_MASK cutoff_mask3 = simd_mask_cond_lt(rsq3,
                                        cutneighsq3);

                                    if (simd_test_any(cutoff_mask0) ||
                                        simd_test_any(cutoff_mask1) ||
                                        simd_test_any(cutoff_mask2) ||
                                        simd_test_any(cutoff_mask3)) {
                                        is_neighbor = 1;
                                    }

#elif defined(CLUSTERPAIR_KERNEL_2XN)

                                MD_SIMD_FLOAT xj_tmp = simd_real_load(&cj_x[CL_X_OFFSET]);
                                MD_SIMD_FLOAT yj_tmp = simd_real_load(&cj_x[CL_Y_OFFSET]);
                                MD_SIMD_FLOAT zj_tmp = simd_real_load(&cj_x[CL_Z_OFFSET]);

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                                MD_SIMD_INT tj_tmp        = simd_i32_load(cj_t);
                                MD_SIMD_INT tvec0         = simd_i32_add(tbase0, tj_tmp);
                                MD_SIMD_INT tvec1         = simd_i32_add(tbase1, tj_tmp);
                                MD_SIMD_FLOAT cutneighsq0 = simd_real_gather(tvec0,
                                    atom->cutneighsq,
                                    sizeof(MD_FLOAT));
                                MD_SIMD_FLOAT cutneighsq1 = simd_real_gather(tvec1,
                                    atom->cutneighsq,
                                    sizeof(MD_FLOAT));
#else
                                MD_SIMD_FLOAT cutneighsq0 = cutneighsq_vec;
                                MD_SIMD_FLOAT cutneighsq1 = cutneighsq_vec;
#endif

                                MD_SIMD_FLOAT delx0 = simd_real_sub(xi0_tmp, xj_tmp);
                                MD_SIMD_FLOAT dely0 = simd_real_sub(yi0_tmp, yj_tmp);
                                MD_SIMD_FLOAT delz0 = simd_real_sub(zi0_tmp, zj_tmp);
                                MD_SIMD_FLOAT delx1 = simd_real_sub(xi1_tmp, xj_tmp);
                                MD_SIMD_FLOAT dely1 = simd_real_sub(yi1_tmp, yj_tmp);
                                MD_SIMD_FLOAT delz1 = simd_real_sub(zi1_tmp, zj_tmp);

                                MD_SIMD_FLOAT rsq0 = simd_real_fma(delx0,
                                    delx0,
                                    simd_real_fma(dely0,
                                        dely0,
                                        simd_real_mul(delz0, delz0)));
                                MD_SIMD_FLOAT rsq1 = simd_real_fma(delx1,
                                    delx1,
                                    simd_real_fma(dely1,
                                        dely1,
                                        simd_real_mul(delz1, delz1)));

                                MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0,
                                    cutneighsq0);
                                MD_SIMD_MASK cutoff_mask1 = simd_mask_cond_lt(rsq1,
                                    cutneighsq1);

                                if (simd_test_any(cutoff_mask0) ||
                                    simd_test_any(cutoff_mask1)) {
                                    is_neighbor = 1;
                                }
#else
                                is_neighbor = 0;
                                for (int cii = 0; cii < CLUSTER_M; cii++) {
                                    for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
                                        MD_FLOAT delx = ci_x[CL_X_INDEX_3D(cii)] -
                                                        cj_x[CL_X_INDEX_3D(cjj)];
                                        MD_FLOAT dely = ci_x[CL_Y_INDEX_3D(cii)] -
                                                        cj_x[CL_Y_INDEX_3D(cjj)];
                                        MD_FLOAT delz = ci_x[CL_Z_INDEX_3D(cii)] -
                                                        cj_x[CL_Z_INDEX_3D(cjj)];

                                        if (delx * delx + dely * dely + delz * delz <
                                            cutneighsq) {
                                            is_neighbor = 1;
                                        }
                                    }
                                }
#endif
                                }

                                if (is_neighbor) {
                                    // We use true (1) for rdiag because we only care if
                                    // there are masks at all, and when this is set to
                                    // false (0) the self-exclusions are not accounted
                                    // for, which  makes the optimized version to not
                                    // work!
                                    unsigned int imask = get_imask_simd(1, ci, cj);

                                    if (n < neighbor->maxneighs) {
                                        if (imask == NBNXN_INTERACTION_MASK_ALL &&
                                            (atom->molsize <= 1 ||
                                                !clusterPairHasExclusion(atom, ci, cj))) {
                                            neighs_padded(neighbor->neighbors,
                                                ci,
                                                n,
                                                nbM,
                                                nbN) = cj;
                                            neighs_padded(neighbor->neighbors_imask,
                                                ci,
                                                n,
                                                nbM,
                                                nbN) = imask;
                                        } else {
                                            neighs_padded(neighbor->neighbors,
                                                ci,
                                                n,
                                                nbM,
                                                nbN) = neighs_padded(neighbor->neighbors,
                                                ci,
                                                nmasked,
                                                nbM,
                                                nbN);
                                            neighs_padded(neighbor->neighbors_imask,
                                                ci,
                                                n,
                                                nbM,
                                                nbN) =
                                                neighs_padded(neighbor->neighbors_imask,
                                                    ci,
                                                    nmasked,
                                                    nbM,
                                                    nbN);
                                            neighs_padded(neighbor->neighbors,
                                                ci,
                                                nmasked,
                                                nbM,
                                                nbN) = cj;
                                            neighs_padded(neighbor->neighbors_imask,
                                                ci,
                                                nmasked,
                                                nbM,
                                                nbN) = imask;
                                            nmasked++;
                                        }
                                    }

                                    n++;
                                }
                            }
                        }

                        m++;
                        if (m < c) {
                            cj       = loc_bin[m];
                            jbb_xmin = atom->jclusters[cj].bbminx;
                            jbb_xmax = atom->jclusters[cj].bbmaxx;
                            jbb_ymin = atom->jclusters[cj].bbminy;
                            jbb_ymax = atom->jclusters[cj].bbmaxy;
                            jbb_zmin = atom->jclusters[cj].bbminz;
                            jbb_zmax = atom->jclusters[cj].bbmaxz;
                        }
                    }
                }
            }

            // Fill neighbor list with dummy values to fit vector width
            /*
            if (CLUSTER_N < VECTOR_WIDTH) {
                while (n % (VECTOR_WIDTH / CLUSTER_N)) {
                    // Last cluster is always a dummy cluster
                    neighs_padded(neighbor->neighbors, ci, n, nbM, nbN) = atom->dummy_cj;
                    neighs_padded(neighbor->neighbors_imask, ci, n, nbM, nbN) = 0;
                    n++;
                }
            }
            */
            neighbor->numneigh[ci]        = n;
            neighbor->numneigh_masked[ci] = nmasked;
            if (n >= neighbor->maxneighs) {
                resize_local = 1;

                if (n > new_maxneighs) {
                    new_maxneighs = n;
                }
            }
        }

        resize = resize_local;

        if (resize) {
            neighbor->maxneighs = new_maxneighs * 1.2;
            fprintf(stdout, "RESIZE %d, PROC %d\n", neighbor->maxneighs, me);
            fflush(stdout);
            free(neighbor->neighbors);
            free(neighbor->neighbors_imask);
            free(is_inner_buf);
            neighbor->neighbors       = (int*)allocate(ALIGNMENT,
                nmax * neighbor->maxneighs * sizeof(int));
            neighbor->neighbors_imask = (unsigned int*)allocate(ALIGNMENT,
                nmax * neighbor->maxneighs * sizeof(unsigned int));
            is_inner_buf              = (int*)allocate(ALIGNMENT,
                (size_t)max_threads() * neighbor->maxneighs * sizeof(int));
#ifdef CUDA_TARGET
            growNeighborCUDA(atom, neighbor);
#endif
        }
    }

#ifdef NBLIST_CSR
    /* Compact padded build buffer into CSR layout */
    {
        const int padded_nbN     = neighbor->maxneighs;
        const int nbM            = atom->Nclusters_local;
        neighbor->neigh_start[0] = 0;
        for (int ci = 0; ci < nbM; ci++) {
            neighbor->neigh_start[ci + 1] = neighbor->neigh_start[ci] +
                                            neighbor->numneigh[ci];
        }
        int total               = neighbor->neigh_start[nbM];
        int* csr_neighbors      = (int*)allocate(ALIGNMENT, MAX(1, total) * sizeof(int));
        unsigned int* csr_imask = (unsigned int*)allocate(ALIGNMENT,
            MAX(1, total) * sizeof(unsigned int));
        for (int ci = 0; ci < nbM; ci++) {
            int nn = neighbor->numneigh[ci];
            for (int k = 0; k < nn; k++) {
                int dst            = neighbor->neigh_start[ci] + k;
                csr_neighbors[dst] = neighs_padded(neighbor->neighbors,
                    ci,
                    k,
                    nbM,
                    padded_nbN);
                csr_imask[dst]     = neighs_padded(neighbor->neighbors_imask,
                    ci,
                    k,
                    nbM,
                    padded_nbN);
            }
        }
        free(neighbor->neighbors);
        free(neighbor->neighbors_imask);
        neighbor->neighbors       = csr_neighbors;
        neighbor->neighbors_imask = (unsigned int*)csr_imask;
        int new_max               = 0;
        for (int ci = 0; ci < nbM; ci++) {
            if (neighbor->numneigh[ci] > new_max) new_max = neighbor->numneigh[ci];
        }
        if (new_max != neighbor->maxneighs) {
            neighbor->maxneighs = new_max;
            free(is_inner_buf);
            is_inner_buf = (int*)allocate(ALIGNMENT,
                (size_t)max_threads() * MAX(1, neighbor->maxneighs) * sizeof(int));
        }
    }
#elif !defined(NBLIST_AOS)
    /* Compact padded build buffer into SOA layout */
    {
        const int padded_nbN    = neighbor->maxneighs;
        const int nbM           = atom->Nclusters_local;
        int* soa_neighbors      = (int*)allocate(ALIGNMENT,
            (size_t)nmax * padded_nbN * sizeof(int));
        unsigned int* soa_imask = (unsigned int*)allocate(ALIGNMENT,
            (size_t)nmax * padded_nbN * sizeof(unsigned int));

        for (int ci = 0; ci < nbM; ci++) {
            int nn = neighbor->numneigh[ci];
            for (int k = 0; k < nn; k++) {
                soa_neighbors[k * nbM + ci] = neighs_padded(neighbor->neighbors,
                    ci,
                    k,
                    nbM,
                    padded_nbN);
                soa_imask[k * nbM + ci]     = neighs_padded(neighbor->neighbors_imask,
                    ci,
                    k,
                    nbM,
                    padded_nbN);
            }
        }

        free(neighbor->neighbors);
        free(neighbor->neighbors_imask);
        neighbor->neighbors       = soa_neighbors;
        neighbor->neighbors_imask = soa_imask;
    }
#endif

    if (method == eightShell) {
        neighborGhost(atom, neighbor);
    }

    // for (int ci = 0; ci < atom->Nclusters_local; ci++) {
    //     printf("ci=%d, numneighs=%d\n", ci, neighbor->numneigh[ci]);
    //     for(int k = 0; k < neighbor->numneigh[ci]; k++) {
    //     const int nbM     = atom->Nclusters_local;
    //     const int nbN     = neighbor->maxneighs;
    //     int cj = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
    //     printf("cj=%d\n", cj);
    //     }
    // }

    /*
    DEBUG_MESSAGE("\ncutneighsq = %f, rbb_sq = %f\n", cutneighsq, rbb_sq);
    for(int ci = 0; ci < 6; ci++) {
        int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT *ci_x = &atom->cl_x[ci_vec_base];

        DEBUG_MESSAGE("Cluster %d, bbx = {%f, %f}, bby = {%f, %f}, bbz = {%f, %f}\n",
            ci,
            atom->iclusters[ci].bbminx,
            atom->iclusters[ci].bbmaxx,
            atom->iclusters[ci].bbminy,
            atom->iclusters[ci].bbmaxy,
            atom->iclusters[ci].bbminz,
            atom->iclusters[ci].bbmaxz);

        for(int cii = 0; cii < CLUSTER_M; cii++) {
            DEBUG_MESSAGE("%f, %f, %f\n", ci_x[CL_X_INDEX_3D(cii)], ci_x[CL_Y_OFFSET +
    cii], ci_x[CL_Z_OFFSET + cii]);
        }

        DEBUG_MESSAGE("Neighbors:\n");
        for(int k = 0; k < neighbor->numneigh[ci]; k++) {
            int cj = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
            int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
            MD_FLOAT *cj_x = &atom->cl_x[cj_vec_base];

            DEBUG_MESSAGE("    Cluster %d, bbx = {%f, %f}, bby = {%f, %f}, bbz = {%f,
    %f}\n", cj, atom->jclusters[cj].bbminx, atom->jclusters[cj].bbmaxx,
                atom->jclusters[cj].bbminy,
                atom->jclusters[cj].bbmaxy,
                atom->jclusters[cj].bbminz,
                atom->jclusters[cj].bbmaxz);

            for(int cjj = 0; cjj < CLUSTER_N; cjj++) {
                DEBUG_MESSAGE("    %f, %f, %f\n", cj_x[CL_X_INDEX_3D(cjj)],
    cj_x[CL_Y_OFFSET + cjj], cj_x[CL_Z_OFFSET + cjj]);
            }
        }
    }
    */

    pruneNeighborCPU(NULL, atom, neighbor);
    DEBUG_MESSAGE("buildNeighbor end\n");
}
