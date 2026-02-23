/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <atom.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <simd.h>
#include <util.h>

#define SMALL  1.0e-6
#define FACTOR 0.999
#define eps    1.0e-9

BuildNeighborFunction buildNeighbor;
PruneNeighborFunction pruneNeighbor;
BuildClustersFunction buildClusters;

static MD_FLOAT xprd, yprd, zprd;
static MD_FLOAT bininvx, bininvy;
static int mbinxlo, mbinylo;
static int nbinx, nbiny;
static int mbinx, mbiny; // n bins in x, y
static int* bincount;
static int* bins;
static int* bin_nclusters;
static int* bin_clusters;
static int mbins;            // total number of bins
static int atoms_per_bin;    // max atoms per bin
static int clusters_per_bin; // max clusters per bin
static MD_FLOAT cutneigh;
static MD_FLOAT cutneighsq; // neighbor cutoff squared
static int nmax;
static int nstencil; // # of bins in stencil
static int* stencil; // stencil list of bin offsets
static MD_FLOAT binsizex, binsizey;

static int coord2bin(MD_FLOAT, MD_FLOAT);
static MD_FLOAT bindist(int, int);
// MPI Implementation
int me;          // rank
int method;      // method
int shellMethod; // If shell method exist
// static int ghostZone(Atom*, int);
static int halfZoneCluster(Atom*, int);
static int ghostClusterinRange(Atom*, int, int, MD_FLOAT);
static void neighborGhost(Atom*, Neighbor*);

/* exported subroutines */
void initNeighbor(Neighbor* neighbor, Parameter* param) {
    MD_FLOAT neighscale = 5.0 / 6.0;
    xprd                = param->nx * param->lattice;
    yprd                = param->ny * param->lattice;
    zprd                = param->nz * param->lattice;
    cutneigh            = param->cutneigh;
    nmax                = 0;
    atoms_per_bin       = 8;
    clusters_per_bin    = (atoms_per_bin / CLUSTER_M) + 10;
    stencil             = NULL;
    bins                = NULL;
    bincount            = NULL;
    bin_clusters        = NULL;
    bin_nclusters       = NULL;
    // neighbor->half_neigh      = param->half_neigh;
    if(param->super_clustering) {
        neighbor->maxneighs       = 800;
    } else {
        neighbor->maxneighs       = 200;
    }
    neighbor->numneigh        = NULL;
    neighbor->numneigh_masked = NULL;
    neighbor->neighbors       = NULL;
    neighbor->neighbors_imask = NULL;
    // MPI Implementation
    method      = param->method;
    shellMethod = method ? 1 : 0;
    if (method) param->half_neigh = 1;
    neighbor->half_neigh = param->half_neigh;
    me                   = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
    // Eight Shell
    neighbor->Nshell        = 0;
    neighbor->numNeighShell = NULL;
    neighbor->neighshell    = NULL;
    neighbor->listshell     = NULL;

    if (param->super_clustering) {
        buildNeighbor = buildNeighborSuperclusters;
        pruneNeighbor = pruneNeighborSuperclusters;
        buildClusters = buildSuperclusters;
    } else {
        buildNeighbor = buildNeighborCPU;
        pruneNeighbor = pruneNeighborCPU;
        buildClusters = buildClustersCPU;
    }
}

void setupNeighbor(Parameter* param, Atom* atom) {
    MD_FLOAT coord;
    int mbinxhi, mbinyhi;
    int nextx, nexty, nextz;

    if (param->input_file != NULL) {
        xprd = param->xprd;
        yprd = param->yprd;
        zprd = param->zprd;
    }

    // TODO: update lo and hi for standard case and use them here instead
    MD_FLOAT xlo = 0.0;
    MD_FLOAT xhi = xprd;
    MD_FLOAT ylo = 0.0;
    MD_FLOAT yhi = yprd;
    MD_FLOAT zlo = 0.0;
    MD_FLOAT zhi = zprd;

    MD_FLOAT atom_density = ((MD_FLOAT)(atom->Natoms)) /
                            ((xhi - xlo) * (yhi - ylo) * (zhi - zlo));
    MD_FLOAT atoms_in_cell = MAX(CLUSTER_M, CLUSTER_N);
    MD_FLOAT targetsizex   = cbrt(atoms_in_cell / atom_density);
    MD_FLOAT targetsizey   = cbrt(atoms_in_cell / atom_density);

    if (param->super_clustering) {
        targetsizex *= (MD_FLOAT)SCLUSTER_SIZE_X;
        targetsizey *= (MD_FLOAT)SCLUSTER_SIZE_Y;
    }

    nbinx      = MAX(1, (int)ceil((xhi - xlo) / targetsizex));
    nbiny      = MAX(1, (int)ceil((yhi - ylo) / targetsizey));
    binsizex   = (xhi - xlo) / nbinx;
    binsizey   = (yhi - ylo) / nbiny;
    bininvx    = 1.0 / binsizex;
    bininvy    = 1.0 / binsizey;
    cutneighsq = cutneigh * cutneigh;

    coord   = xlo - cutneigh - SMALL * xprd;
    mbinxlo = (int)(coord * bininvx);
    if (coord < 0.0) {
        mbinxlo = mbinxlo - 1;
    }
    coord   = xhi + cutneigh + SMALL * xprd;
    mbinxhi = (int)(coord * bininvx);

    coord   = ylo - cutneigh - SMALL * yprd;
    mbinylo = (int)(coord * bininvy);
    if (coord < 0.0) {
        mbinylo = mbinylo - 1;
    }
    coord   = yhi + cutneigh + SMALL * yprd;
    mbinyhi = (int)(coord * bininvy);

    mbinxlo = mbinxlo - 1;
    mbinxhi = mbinxhi + 1;
    mbinx   = mbinxhi - mbinxlo + 1;

    mbinylo = mbinylo - 1;
    mbinyhi = mbinyhi + 1;
    mbiny   = mbinyhi - mbinylo + 1;

    nextx = (int)(cutneigh * bininvx);
    nexty = (int)(cutneigh * bininvy);
    if (nextx * binsizex < FACTOR * cutneigh) nextx++;
    if (nexty * binsizey < FACTOR * cutneigh) nexty++;

    if (stencil) {
        free(stencil);
    }
    stencil  = (int*)malloc((2 * nexty + 1) * (2 * nextx + 1) * sizeof(int));
    nstencil = 0;

    for (int j = -nexty; j <= nexty; j++) {
        for (int i = -nextx; i <= nextx; i++) {
            if (bindist(i, j) < cutneighsq) {
                stencil[nstencil++] = j * mbinx + i;
            }
        }
    }

    if (bincount) {
        free(bincount);
    }
    if (bins) {
        free(bins);
    }
    if (bin_nclusters) {
        free(bin_nclusters);
    }
    if (bin_clusters) {
        free(bin_clusters);
    }
    mbins         = mbinx * mbiny;
    bincount      = (int*)malloc(mbins * sizeof(int));
    bins          = (int*)malloc(mbins * atoms_per_bin * sizeof(int));
    bin_nclusters = (int*)malloc(mbins * sizeof(int));
    bin_clusters  = (int*)malloc(mbins * clusters_per_bin * sizeof(int));

    /*
    DEBUG_MESSAGE("lo, hi = (%e, %e, %e), (%e, %e, %e)\n", xlo, ylo, zlo, xhi, yhi, zhi);
    DEBUG_MESSAGE("binsize = %e, %e\n", binsizex, binsizey);
    DEBUG_MESSAGE("mbin lo, hi = (%d, %d), (%d, %d)\n", mbinxlo, mbinylo, mbinxhi,
    mbinyhi); DEBUG_MESSAGE("mbins = %d (%d x %d)\n", mbins, mbinx, mbiny);
    DEBUG_MESSAGE("nextx = %d, nexty = %d\n", nextx, nexty);
    */
}

MD_FLOAT getBoundingBoxDistanceSq(Atom* atom, int ci, int cj) {
    MD_FLOAT dl  = atom->iclusters[ci].bbminx - atom->jclusters[cj].bbmaxx;
    MD_FLOAT dh  = atom->jclusters[cj].bbminx - atom->iclusters[ci].bbmaxx;
    MD_FLOAT dm  = MAX(dl, dh);
    MD_FLOAT dm0 = MAX(dm, 0.0);
    MD_FLOAT d2  = dm0 * dm0;

    dl  = atom->iclusters[ci].bbminy - atom->jclusters[cj].bbmaxy;
    dh  = atom->jclusters[cj].bbminy - atom->iclusters[ci].bbmaxy;
    dm  = MAX(dl, dh);
    dm0 = MAX(dm, 0.0);
    d2 += dm0 * dm0;

    dl  = atom->iclusters[ci].bbminz - atom->jclusters[cj].bbmaxz;
    dh  = atom->jclusters[cj].bbminz - atom->iclusters[ci].bbmaxz;
    dm  = MAX(dl, dh);
    dm0 = MAX(dm, 0.0);
    d2 += dm0 * dm0;
    return d2;
}

/* Returns a diagonal or off-diagonal interaction mask for plain C lists */
static unsigned int get_imask(int rdiag, int ci, int cj) {
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=2 */
static unsigned int get_imask_simd_j2(int rdiag, int ci, int cj) {
    return (rdiag && ci * 2 == cj
                ? NBNXN_INTERACTION_MASK_DIAG_J2_0
                : (rdiag && ci * 2 + 1 == cj ? NBNXN_INTERACTION_MASK_DIAG_J2_1
                                             : NBNXN_INTERACTION_MASK_ALL));
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=4 */
static unsigned int get_imask_simd_j4(int rdiag, int ci, int cj) {
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

/* Returns a diagonal or off-diagonal interaction mask for cj-size=8 */
static unsigned int get_imask_simd_j8(int rdiag, int ci, int cj) {
    return (rdiag && ci == cj * 2
                ? NBNXN_INTERACTION_MASK_DIAG_J8_0
                : (rdiag && ci == cj * 2 + 1 ? NBNXN_INTERACTION_MASK_DIAG_J8_1
                                             : NBNXN_INTERACTION_MASK_ALL));
}

#if VECTOR_WIDTH == 2
#define get_imask_simd_4xn get_imask_simd_j2
#elif VECTOR_WIDTH == 4
#define get_imask_simd_4xn get_imask_simd_j4
#elif VECTOR_WIDTH == 8
#define get_imask_simd_4xn  get_imask_simd_j8
#define get_imask_simd_2xnn get_imask_simd_j4
#elif VECTOR_WIDTH == 16
#define get_imask_simd_2xnn get_imask_simd_j8
#else
#error "Invalid cluster configuration"
#endif

void buildNeighborCPU(Atom* atom, Neighbor* neighbor) {
    DEBUG_MESSAGE("buildNeighbor start\n");

    /* extend atom arrays if necessary */
    if (atom->Nclusters_local > nmax) {
        nmax = atom->Nclusters_local;
        if (neighbor->numneigh) free(neighbor->numneigh);
        if (neighbor->numneigh_masked) free(neighbor->numneigh_masked);
        if (neighbor->neighbors) free(neighbor->neighbors);
        if (neighbor->neighbors_imask) free(neighbor->neighbors_imask);
        neighbor->numneigh        = (int*)malloc(nmax * sizeof(int));
        neighbor->numneigh_masked = (int*)malloc(nmax * sizeof(int));
        neighbor->neighbors = (int*)malloc(nmax * neighbor->maxneighs * sizeof(int));
        neighbor->neighbors_imask = (unsigned int*)malloc(
            nmax * neighbor->maxneighs * sizeof(unsigned int));
    }

    MD_FLOAT bbx    = 0.5 * (binsizex + binsizex);
    MD_FLOAT bby    = 0.5 * (binsizey + binsizey);
    MD_FLOAT rbb_sq = MAX(0.0, cutneigh - 0.5 * sqrt(bbx * bbx + bby * bby));
    rbb_sq          = rbb_sq * rbb_sq;
    int resize      = 1;

    /* loop over each atom, storing neighbors */
    while (resize) {
        const int nbM     = atom->Nclusters_local;
        const int nbN     = neighbor->maxneighs;
        int new_maxneighs = neighbor->maxneighs;
        resize            = 0;
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_cj0    = CJ0_FROM_CI(ci);
            int n = 0, nmasked = 0;
            int ibin        = atom->cluster_bin[ci];
            int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x  = &atom->cl_x[ci_vec_base];

#ifndef ONE_ATOM_TYPE
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

#ifndef ONE_ATOM_TYPE
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

#ifndef ONE_ATOM_TYPE
            MD_SIMD_INT tbase0    = simd_i32_broadcast(ci_t[0] * atom->ntypes);
            MD_SIMD_INT tbase1    = simd_i32_broadcast(ci_t[1] * atom->ntypes);
            MD_SIMD_INT tbase2    = simd_i32_broadcast(ci_t[2] * atom->ntypes);
            MD_SIMD_INT tbase3    = simd_i32_broadcast(ci_t[3] * atom->ntypes);
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
                    MD_FLOAT dl, dh, dm, dm0, d_bb_sq;

                    do {
                        m++;
                        cj = loc_bin[m];
                        if (neighbor->half_neigh && ci_cj0 > cj) {
                            continue;
                        }
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
                        if (!neighbor->half_neigh || ci_cj0 <= cj) {
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

#ifndef ONE_ATOM_TYPE
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

#ifndef ONE_ATOM_TYPE
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

                                    MD_SIMD_FLOAT xj_tmp = simd_real_load(&cj_x[CL_X_INDEX_3D(0)]);
                                    MD_SIMD_FLOAT yj_tmp = simd_real_load(&cj_x[CL_Y_INDEX_3D(0)]);
                                    MD_SIMD_FLOAT zj_tmp = simd_real_load(&cj_x[CL_Z_INDEX_3D(0)]);
#ifndef ONE_ATOM_TYPE
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
                                    unsigned int imask;
#ifdef CLUSTERPAIR_KERNEL_2XNN
                                    imask = get_imask_simd_2xnn(1, ci, cj);
#else
                                    imask = get_imask_simd_4xn(1, ci, cj);
#endif

                                    if (n < neighbor->maxneighs) {
                                        if (imask == NBNXN_INTERACTION_MASK_ALL) {
                                            neighs(neighbor->neighbors, ci, n, nbM, nbN) = cj;
                                            neighs(neighbor->neighbors_imask, ci, n, nbM, nbN) = imask;
                                        } else {
                                            neighs(neighbor->neighbors, ci, n, nbM, nbN) =
                                                                    neighs(neighbor->neighbors, ci, nmasked, nbM, nbN);
                                            neighs(neighbor->neighbors_imask, ci, n, nbM, nbN) =
                                                                    neighs(neighbor->neighbors_imask, ci, nmasked, nbM, nbN);
                                            neighs(neighbor->neighbors, ci, nmasked, nbM, nbN) = cj;
                                            neighs(neighbor->neighbors_imask, ci, nmasked, nbM, nbN) = imask;
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
                    neighs(neighbor->neighbors, ci, n, nbM, nbN) = atom->dummy_cj;
                    neighs(neighbor->neighbors_imask, ci, n, nbM, nbN) = 0;
                    n++;
                }
            }
            */
            neighbor->numneigh[ci]        = n;
            neighbor->numneigh_masked[ci] = nmasked;
            if (n >= neighbor->maxneighs) {
                resize = 1;

                if (n >= new_maxneighs) {
                    new_maxneighs = n;
                }
            }
        }

        if (resize) {
            neighbor->maxneighs = new_maxneighs * 1.2;
            fprintf(stdout, "RESIZE %d, PROC %d\n", neighbor->maxneighs, me);
            fflush(stdout);
            free(neighbor->neighbors);
            free(neighbor->neighbors_imask);
            neighbor->neighbors = (int*)malloc(nmax * neighbor->maxneighs * sizeof(int));
            neighbor->neighbors_imask = (unsigned int*)malloc(
                nmax * neighbor->maxneighs * sizeof(unsigned int));
#ifdef CUDA_TARGET
            growNeighborCUDA(atom, neighbor);
#endif
        }
    }

    if (method == eightShell) {
        neighborGhost(atom, neighbor);
    }

    // for (int ci = 0; ci < atom->Nclusters_local; ci++) {
    //     printf("ci=%d, numneighs=%d\n", ci, neighbor->numneigh[ci]);
    //     for(int k = 0; k < neighbor->numneigh[ci]; k++) {
    //     const int nbM     = atom->Nclusters_local;
    //     const int nbN     = neighbor->maxneighs;
    //     int cj = neighs(neighbor->neighbors, ci, k, nbM, nbN);
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
            int cj = neighs(neighbor->neighbors, ci, k, nbM, nbN);
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

    DEBUG_MESSAGE("buildNeighbor end\n");
}

// TODO For future parallelization on GPU
void buildNeighborSuperclusters(Atom* atom, Neighbor* neighbor) {
    DEBUG_MESSAGE("buildNeighborSuperclusters start\n");

    /* extend atom arrays if necessary */
    if (atom->Nclusters_local > nmax) {
        nmax = atom->Nclusters_local;
        if (neighbor->numneigh) free(neighbor->numneigh);
        if (neighbor->neighbors) free(neighbor->neighbors);
        neighbor->numneigh  = (int*)malloc(nmax * sizeof(int));
        neighbor->neighbors = (int*)malloc(nmax * neighbor->maxneighs * sizeof(int*));
    }

    MD_FLOAT bbx    = 0.5 * (binsizex + binsizex);
    MD_FLOAT bby    = 0.5 * (binsizey + binsizey);
    MD_FLOAT rbb_sq = MAX(0.0, cutneigh - 0.5 * sqrt(bbx * bbx + bby * bby));
    rbb_sq          = rbb_sq * rbb_sq;
    int resize      = 1;

    /* loop over each atom, storing neighbors */
    while (resize) {
        const int nbM     = atom->Nclusters_local;
        const int nbN     = neighbor->maxneighs;
        int new_maxneighs = neighbor->maxneighs;
        resize            = 0;

        for (int sci = 0; sci < atom->Nclusters_local; sci++) {
            const int sci_vec_base = SCI_VECTOR_BASE_INDEX(sci);
            int n             = 0;
            int ibin          = atom->cluster_bin[sci];
            MD_FLOAT ibb_xmin = atom->siclusters[sci].bbminx;
            MD_FLOAT ibb_xmax = atom->siclusters[sci].bbmaxx;
            MD_FLOAT ibb_ymin = atom->siclusters[sci].bbminy;
            MD_FLOAT ibb_ymax = atom->siclusters[sci].bbmaxy;
            MD_FLOAT ibb_zmin = atom->siclusters[sci].bbminz;
            MD_FLOAT ibb_zmax = atom->siclusters[sci].bbmaxz;

            for (int k = 0; k < nstencil; k++) {
                int jbin     = ibin + stencil[k];
                int* loc_bin = &bin_clusters[jbin * clusters_per_bin];
                int cj, m = -1;
                MD_FLOAT jbb_xmin, jbb_xmax, jbb_ymin, jbb_ymax, jbb_zmin, jbb_zmax;
                const int c = bin_nclusters[jbin];

                if (c > 0) {
                    MD_FLOAT dl, dh, dm, dm0, d_bb_sq;

                    do {
                        m++;
                        cj = loc_bin[m];
                        if (neighbor->half_neigh && sci > SCI_FROM_CJ(cj)) {
                            continue;
                        }

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
                        if (!neighbor->half_neigh || sci <= SCI_FROM_CJ(cj)) {
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
                                    for (int sci_ci = 0; sci_ci < atom->siclusters[sci].nclusters; sci_ci++) {
                                        const int ci = sci * SCLUSTER_SIZE + sci_ci;
                                        int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
#ifdef SOA_SUP
                                        MD_FLOAT* ci_x  = &atom->cl_x[sci_vec_base + sci_ci * CLUSTER_M];                             
#else                                     
                                        MD_FLOAT* ci_x  = &atom->cl_x[sci_vec_base + sci_ci * CLUSTER_M * ATOM_DIM];
#endif
                                        MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];

                                        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
                                            for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
                                                MD_FLOAT delx = ci_x[CL_X_INDEX(cii)] - cj_x[CL_X_INDEX(cjj)];
                                                MD_FLOAT dely = ci_x[CL_Y_INDEX(cii)] - cj_x[CL_Y_INDEX(cjj)];
                                                MD_FLOAT delz = ci_x[CL_Z_INDEX(cii)] - cj_x[CL_Z_INDEX(cjj)];

                                                if (delx * delx + dely * dely + delz * delz < cutneighsq) {
                                                    is_neighbor = 1;
                                                    goto distance_check_out;
                                                }
                                            }
                                        }
                                    }
                                }

distance_check_out:
                                if (is_neighbor) {
                                    neighs(neighbor->neighbors, sci, n, nbM, nbN) = cj;
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
            if (CLUSTER_N < VECTOR_WIDTH) {
                while (n % (VECTOR_WIDTH / CLUSTER_N)) {
                    // Last cluster is always a dummy cluster
                    neighs(neighbor->neighbors, sci, n, nbM, nbN) = atom->dummy_cj;
                    n++;
                }
            }

            neighbor->numneigh[sci] = n;
            if (n >= neighbor->maxneighs) {
                resize = 1;

                if (n >= new_maxneighs) {
                    new_maxneighs = n;
                }
            }
        }

        if (resize) {
            neighbor->maxneighs = new_maxneighs * 1.2;
            fprintf(stdout, "RESIZE %d\n", neighbor->maxneighs);
            free(neighbor->neighbors);
            neighbor->neighbors = (int*)malloc(
                atom->Nmax * neighbor->maxneighs * sizeof(int));
        }
    }
    // if (atom->Nclusters_local > 0) debug_check_supercluster_neighbors(atom, neighbor, 0);


    // int sci = 0;
    //  for (int sci = 0; sci < atom->Nclusters_local; sci++) {
    //     int total_neighs = neighbor->numneigh[sci];
    //     int ghost_neighs = 0;
    //     int local_neighs = 0;
    //     for(int k = 0; k < neighbor->numneigh[sci]; k++) {
    //         const int nbM     = atom->Nclusters_local;
    //         const int nbN     = neighbor->maxneighs;
    //         int cj = neighs(neighbor->neighbors, sci, k, nbM, nbN);
    //         if(cj >= (atom->Nclusters_local * SCLUSTER_SIZE)){
    //             ghost_neighs++;
    //         }else{
    //             local_neighs++;
    //         }
    //     }
   
    // printf("sci=%d, total neighs:%d, local:%d, ghost:%d\n", sci, total_neighs, local_neighs, ghost_neighs);
    // }
//  printf("\n");

//     bool cj_check[MAX_CJ_VALUE] = {0}; 
    

//     bool cj_out_of_range = false; 
//   int sci = 0;
//     printf("sci=%d, numneighs=%d\n", sci, neighbor->numneigh[sci]);


    
//     for(int k = 0; k < neighbor->numneigh[sci]; k++) {
//         const int nbM     = atom->Nclusters_local;
//         const int nbN     = neighbor->maxneighs;
//         int cj = neighs(neighbor->neighbors, sci, k, nbM, nbN);

//         // printf("cj=%d\n", cj);

//         if(cj >= (atom->Nclusters_local * SCLUSTER_SIZE)){
//             printf("👻 cj is ghost cluster: %d\n", cj);
//         }

        
//         if (cj >= 0 && cj < MAX_CJ_VALUE) {

//             if (cj_check[cj]) {
//                 fprintf(stderr, "🚨 cj repeated\n", sci, cj);

//             } else {
//                 cj_check[cj] = true;
//             }
//         } else {
//             // cj 值不合法或超出预设范围
//             if (!cj_out_of_range) {
//                 fprintf(stderr, "⚠️ increase MAX_CJ_VALUE size。\n", cj, MAX_CJ_VALUE);
//                 cj_out_of_range = true;
//             }
//         }
//     }

  
    //  for (int sci = 0; sci < atom->Nclusters_local; sci++) {
    // const int sci_vec_base = SCI_VECTOR_BASE_INDEX(sci);
    // for (int sci_ci = 0; sci_ci < atom->siclusters[sci].nclusters; sci_ci++){
    //     const int ci = sci * SCLUSTER_SIZE + sci_ci;
    //     MD_FLOAT* ci_x  = &atom->cl_x[sci_vec_base + sci_ci * CLUSTER_M * ATOM_DIM];
    //     for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
    //         int id = ci *  CLUSTER_M + cii;
    //         printf("atom %d: x=%f y=%f z=%f\n", id, ci_x[CL_X_INDEX_3D(cii)], ci_x[CL_Y_INDEX_3D(cii)], ci_x[CL_Z_INDEX_3D(cii)]);
    //     }
        
        // printf("\n");
        // printf("Cluster %d, bbx = {%f, %f}, bby = {%f, %f}, bbz = {%f, %f}\n",
        //     ci,
        //     atom->iclusters[ci].bbminx,
        //     atom->iclusters[ci].bbmaxx,
        //     atom->iclusters[ci].bbminy,
        //     atom->iclusters[ci].bbmaxy,
        //     atom->iclusters[ci].bbminz,
        //     atom->iclusters[ci].bbmaxz);
        // }
        // }
    // }
    // exit(0);

    // printf("\n");
    // printf("sci=%d, numneighs=%d\n", sci, neighbor->numneigh[sci]);
    // for(int k = 0; k < neighbor->numneigh[sci]; k++) {
    //     const int nbM     = atom->Nclusters_local;
    //     const int nbN     = neighbor->maxneighs;
    //     int cj = neighs(neighbor->neighbors, sci, k, nbM, nbN);
    //     // int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
    //     // MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
    //     if(cj==1440){
    //         for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
    //             int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
    //             MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];
    //             printf("%f \n",cj_x[CL_Z_INDEX_3D(cjj)]);
    //         }
    //     }

    //     // printf("cj=%d\n", cj);
    // }
    // // }
    // printf("buildNeighborSuperclusters end\n");
    // exit(0);

    /*
    DEBUG_MESSAGE("\ncutneighsq = %f, rbb_sq = %f\n", cutneighsq, rbb_sq);
    for(int ci = 0; ci < 6; ci++) {
        int ci_vec_base = CI_VECTOR_BASE_INDEX(ci);
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
            DEBUG_MESSAGE("%f, %f, %f\n", ci_x[CL_X_INDEX(cii)], ci_x[CL_Y_INDEX(cii)], ci_x[CL_Z_INDEX(cii)]);
        }

        DEBUG_MESSAGE("Neighbors:\n");
        for(int k = 0; k < neighbor->numneigh[ci]; k++) {
            int cj = neighs(neighbor->neighbors, ci, k, nbM, nbN);
            int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
            MD_FLOAT *cj_x = &atom->cl_x[cj_vec_base];

            DEBUG_MESSAGE("    Cluster %d, bbx = {%f, %f}, bby = {%f, %f}, bbz = {%f,%f}\n", cj, atom->jclusters[cj].bbminx, atom->jclusters[cj].bbmaxx,
                atom->jclusters[cj].bbminy,
                atom->jclusters[cj].bbmaxy,
                atom->jclusters[cj].bbminz,
                atom->jclusters[cj].bbmaxz);

            for(int cjj = 0; cjj < CLUSTER_N; cjj++) {
                DEBUG_MESSAGE("    %f, %f, %f\n", cj_x[CL_X_INDEX(cjj)],cj_x[CL_Y_INDEX(cjj)], cj_x[CL_Z_INDEX(cjj)]);
            }
        }
    }
    */

    DEBUG_MESSAGE("buildNeighborSuperclusters end\n");
}

void pruneNeighborCPU(Parameter* param, Atom* atom, Neighbor* neighbor) {
    DEBUG_MESSAGE("pruneNeighbor start\n");
    // MD_FLOAT cutsq = param->cutforce * param->cutforce;
    MD_FLOAT cutsq = cutneighsq;
    const int nbM = atom->Nclusters_local;
    const int nbN = neighbor->maxneighs;

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int numneighs              = neighbor->numneigh[ci];
        int numneighs_masked       = neighbor->numneigh_masked[ci];
        int k                      = 0;
        int ci_vec_base            = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ci_x             = &atom->cl_x[ci_vec_base];

#if defined(CLUSTERPAIR_KERNEL_2XNN)
        MD_SIMD_FLOAT cutneighsq_vec = simd_real_broadcast(cutsq);
        MD_SIMD_FLOAT xi0_tmp        = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(0)]);
        MD_SIMD_FLOAT xi2_tmp        = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(2)]);
        MD_SIMD_FLOAT yi0_tmp        = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(0)]);
        MD_SIMD_FLOAT yi2_tmp        = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(2)]);
        MD_SIMD_FLOAT zi0_tmp        = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(0)]);
        MD_SIMD_FLOAT zi2_tmp        = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(2)]);
#elif defined(CLUSTERPAIR_KERNEL_4XN)
        MD_SIMD_FLOAT cutneighsq_vec = simd_real_broadcast(cutsq);
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
#endif

        // Remove dummy clusters if necessary
        if (CLUSTER_N < VECTOR_WIDTH) {
            while (neighs(neighbor->neighbors, ci, numneighs - 1, nbM, nbN) == atom->dummy_cj) {
                numneighs--;
            }
        }

        while (k < numneighs) {
            int cj                 = neighs(neighbor->neighbors, ci, k, nbM, nbN);
            int cj_vec_base        = CJ_VECTOR3_BASE_INDEX(cj);
            MD_FLOAT* cj_x         = &atom->cl_x[cj_vec_base];
            int atom_dist_in_range = 0;

#if defined(CLUSTERPAIR_KERNEL_2XNN)

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

            MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0, cutneighsq_vec);
            MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2, cutneighsq_vec);

            if (simd_test_any(cutoff_mask0) || simd_test_any(cutoff_mask2)) {
                atom_dist_in_range = 1;
            }

#elif defined(CLUSTERPAIR_KERNEL_4XN)

            MD_SIMD_FLOAT xj_tmp = simd_real_load(&cj_x[CL_X_INDEX_3D(0)]);
            MD_SIMD_FLOAT yj_tmp = simd_real_load(&cj_x[CL_Y_INDEX_3D(0)]);
            MD_SIMD_FLOAT zj_tmp = simd_real_load(&cj_x[CL_Z_INDEX_3D(0)]);
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

            MD_SIMD_MASK cutoff_mask0 = simd_mask_cond_lt(rsq0, cutneighsq_vec);
            MD_SIMD_MASK cutoff_mask1 = simd_mask_cond_lt(rsq1, cutneighsq_vec);
            MD_SIMD_MASK cutoff_mask2 = simd_mask_cond_lt(rsq2, cutneighsq_vec);
            MD_SIMD_MASK cutoff_mask3 = simd_mask_cond_lt(rsq3, cutneighsq_vec);

            if (simd_test_any(cutoff_mask0) || simd_test_any(cutoff_mask1) ||
                simd_test_any(cutoff_mask2) || simd_test_any(cutoff_mask3)) {
                atom_dist_in_range = 1;
            }
#else
            for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
                for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
                    MD_FLOAT delx = ci_x[CL_X_INDEX_3D(cii)] - cj_x[CL_X_INDEX_3D(cjj)];
                    MD_FLOAT dely = ci_x[CL_Y_INDEX_3D(cii)] - cj_x[CL_Y_INDEX_3D(cjj)];
                    MD_FLOAT delz = ci_x[CL_Z_INDEX_3D(cii)] - cj_x[CL_Z_INDEX_3D(cjj)];
                    if (delx * delx + dely * dely + delz * delz < cutsq) {
                        atom_dist_in_range = 1;
                        break;
                    }
                }
            }
#endif

            if (atom_dist_in_range) {
                k++;
            } else {
                numneighs--;
                if (k < numneighs_masked) {
                    numneighs_masked--;
                }
                neighs(neighbor->neighbors, ci, k, nbM, nbN) = neighs(neighbor->neighbors, ci, numneighs, nbM, nbN);
            }
        }

        // Readd dummy clusters if necessary
        if (CLUSTER_N < VECTOR_WIDTH) {
            while (numneighs % (VECTOR_WIDTH / CLUSTER_N)) {
                neighs(neighbor->neighbors, ci, numneighs, nbM, nbN) =
                    atom->dummy_cj; // Last cluster is always a dummy cluster
                neighs(neighbor->neighbors_imask, ci, numneighs, nbM, nbN) = 0;
                numneighs++;
            }
        }

        neighbor->numneigh[ci]        = numneighs;
        neighbor->numneigh_masked[ci] = numneighs_masked;
    }

    DEBUG_MESSAGE("pruneNeighbor end\n");
}

void pruneNeighborSuperclusters(Parameter* param, Atom* atom, Neighbor* neighbor) {
    DEBUG_MESSAGE("pruneNeighbor start\n");
    // MD_FLOAT cutsq = param->cutforce * param->cutforce;
    MD_FLOAT cutsq = cutneighsq;
    const int nbM = atom->Nclusters_local;
    const int nbN = neighbor->maxneighs;

    for (int sci = 0; sci < atom->Nclusters_local; sci++) {
        for (int scii = 0; scii < atom->siclusters[sci].nclusters; scii++) {
            int numneighs = neighbor->numneigh[sci];
            int k         = 0;

            // Remove dummy clusters if necessary
            if (CLUSTER_N < VECTOR_WIDTH) {
                while (neighs(neighbor->neighbors, sci, numneighs - 1, nbM, nbN) == atom->dummy_cj) {
                    numneighs--;
                }
            }

            while (k < numneighs) {
                int is_neighbor = 0;
                int cj          = neighs(neighbor->neighbors, sci, k, nbM, nbN);

                for (int sci_ci = 0; sci_ci < atom->siclusters[sci].nclusters; sci_ci++) {
                    const int ci    = sci * SCLUSTER_SIZE + sci_ci;
                    int ci_vec_base = SCI_VECTOR_BASE_INDEX(sci) + sci_ci * CLUSTER_M;
                    int cj_vec_base = CJ_VECTOR3_BASE_INDEX(cj);
                    MD_FLOAT* ci_x  = &atom->cl_x[ci_vec_base];
                    MD_FLOAT* cj_x  = &atom->cl_x[cj_vec_base];

                    for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
                        for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
                            MD_FLOAT delx = ci_x[CL_X_INDEX(cii)] -
                                            cj_x[CL_X_INDEX(cjj)];
                            MD_FLOAT dely = ci_x[CL_Y_INDEX(cii)] -
                                            cj_x[CL_Y_INDEX(cjj)];
                            MD_FLOAT delz = ci_x[CL_Z_INDEX(cii)] -
                                            cj_x[CL_Z_INDEX(cjj)];
                            if (delx * delx + dely * dely + delz * delz < cutsq) {
                                is_neighbor = 1;
                                break;
                            }
                        }
                    }
                }

                if (is_neighbor) {
                    k++;
                } else {
                    numneighs--;
                    neighs(neighbor->neighbors, sci, k, nbM, nbN) = neighs(neighbor->neighbors, sci, numneighs, nbM, nbN);
                }
            }

            // Readd dummy clusters if necessary
            if (CLUSTER_N < VECTOR_WIDTH) {
                while (numneighs % (VECTOR_WIDTH / CLUSTER_N)) {
                    neighs(neighbor->neighbors, sci, numneighs, nbM, nbN) =
                        atom->dummy_cj; // Last cluster is always a dummy cluster
                    neighs(neighbor->neighbors_imask, sci, numneighs, nbM, nbN) = 0;
                    numneighs++;
                }
            }

            neighbor->numneigh[sci] = numneighs;
        }
    }

    DEBUG_MESSAGE("pruneNeighbor end\n");
}

/* internal subroutines */
MD_FLOAT bindist(int i, int j) {
    MD_FLOAT delx, dely, delz;

    if (i > 0) {
        delx = (i - 1) * binsizex;
    } else if (i == 0) {
        delx = 0.0;
    } else {
        delx = (i + 1) * binsizex;
    }

    if (j > 0) {
        dely = (j - 1) * binsizey;
    } else if (j == 0) {
        dely = 0.0;
    } else {
        dely = (j + 1) * binsizey;
    }

    return (delx * delx + dely * dely);
}

int coord2bin(MD_FLOAT xin, MD_FLOAT yin) {
    int ix, iy;

    if (xin >= xprd) {
        ix = (int)((xin - xprd) * bininvx) + nbinx - mbinxlo;
    } else if (xin >= 0.0) {
        ix = (int)(xin * bininvx) - mbinxlo;
    } else {
        ix = (int)(xin * bininvx) - mbinxlo - 1;
    }

    if (yin >= yprd) {
        iy = (int)((yin - yprd) * bininvy) + nbiny - mbinylo;
    } else if (yin >= 0.0) {
        iy = (int)(yin * bininvy) - mbinylo;
    } else {
        iy = (int)(yin * bininvy) - mbinylo - 1;
    }

    return (iy * mbinx + ix + 1);
}

void coord2bin2D(MD_FLOAT xin, MD_FLOAT yin, int* ix, int* iy) {
    if (xin >= xprd) {
        *ix = (int)((xin - xprd) * bininvx) + nbinx - mbinxlo;
    } else if (xin >= 0.0) {
        *ix = (int)(xin * bininvx) - mbinxlo;
    } else {
        *ix = (int)(xin * bininvx) - mbinxlo - 1;
    }

    if (yin >= yprd) {
        *iy = (int)((yin - yprd) * bininvy) + nbiny - mbinylo;
    } else if (yin >= 0.0) {
        *iy = (int)(yin * bininvy) - mbinylo;
    } else {
        *iy = (int)(yin * bininvy) - mbinylo - 1;
    }
}

void binAtoms(Atom* atom) {
    DEBUG_MESSAGE("binAtoms start\n");
    int resize = 1;

    while (resize > 0) {
        resize = 0;

        for (int i = 0; i < mbins; i++) {
            bincount[i] = 0;
        }

        for (int i = 0; i < atom->Nlocal; i++) {
            int ibin = coord2bin(atom_x(i), atom_y(i));
            if (ibin < 0 || ibin >= mbins) {
                fprintf(
                    stderr,
                    "Error: Atom %d assigned to invalid bin %d, position: %.4f, %.4f\n",
                    i, ibin, atom_x(i), atom_y(i));

                exit(EXIT_FAILURE);
            }
            if (bincount[ibin] < atoms_per_bin) {
                int ac                          = bincount[ibin]++;
                bins[ibin * atoms_per_bin + ac] = i;
            } else {
                resize = 1;
            }
        }

        if (resize) {
            free(bins);
            atoms_per_bin *= 2;
            bins = (int*)malloc(mbins * atoms_per_bin * sizeof(int));
        }
    }

    DEBUG_MESSAGE("binAtoms end\n");
}

// TODO: Use pigeonhole sorting
void sortAtomsByZCoord(Atom* atom) {
    DEBUG_MESSAGE("sortAtomsByZCoord start\n");
    for (int bin = 0; bin < mbins; bin++) {
        int c        = bincount[bin];
        int* bin_ptr = &bins[bin * atoms_per_bin];

        for (int ac_i = 0; ac_i < c; ac_i++) {
            int i          = bin_ptr[ac_i];
            int min_ac     = ac_i;
            int min_idx    = i;
            MD_FLOAT min_z = atom_z(i);

            for (int ac_j = ac_i + 1; ac_j < c; ac_j++) {
                int j       = bin_ptr[ac_j];
                MD_FLOAT zj = atom_z(j);
                if (zj < min_z) {
                    min_ac  = ac_j;
                    min_idx = j;
                    min_z   = zj;
                }
            }

            bin_ptr[ac_i]   = min_idx;
            bin_ptr[min_ac] = i;
        }
    }

    DEBUG_MESSAGE("sortAtomsByZCoord end\n");
}

// TODO: Use pigeonhole sorting
void sortAtomsByCoord(Atom* atom, int dim, int bin, int start_index, int end_index) {
    // DEBUG_MESSAGE("sortAtomsByCoord start\n");
    int* bin_ptr = &bins[bin * atoms_per_bin];

    for (int ac_i = start_index; ac_i <= end_index; ac_i++) {
        int i              = bin_ptr[ac_i];
        int min_ac         = ac_i;
        int min_idx        = i;
        MD_FLOAT min_coord = (dim == 0) ? atom_x(i) :
                             (dim == 1) ? atom_y(i) :
                                          atom_z(i);

        for (int ac_j = ac_i + 1; ac_j <= end_index; ac_j++) {
            int j           = bin_ptr[ac_j];
            MD_FLOAT coordj = (dim == 0) ? atom_x(j) :
                              (dim == 1) ? atom_y(j) :
                                           atom_z(j);

            if (coordj < min_coord) {
                min_ac    = ac_j;
                min_idx   = j;
                min_coord = coordj;
            }
        }

        bin_ptr[ac_i]   = min_idx;
        bin_ptr[min_ac] = i;
    }

    // DEBUG_MESSAGE("sortAtomsByCoord end\n");
}

void buildClustersCPU(Atom* atom) {
    DEBUG_MESSAGE("buildClusters start\n");
    atom->Nclusters_local = 0;

    /* bin local atoms */
    binAtoms(atom);
    sortAtomsByZCoord(atom);

    for (int bin = 0; bin < mbins; bin++) {
        int c         = bincount[bin];
        int ac        = 0;
        int nclusters = ((c + CLUSTER_M - 1) / CLUSTER_M);
        if (CLUSTER_N > CLUSTER_M && nclusters % 2) {
            nclusters++;
        }
        for (int cl = 0; cl < nclusters; cl++) {
            const int ci = atom->Nclusters_local;
            if (ci >= atom->Nclusters_max) {
                growClusters(atom, 0);
            }

            int ci_sca_base = CI_SCALAR_BASE_INDEX(ci);
            int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
            MD_FLOAT* ci_x  = &atom->cl_x[ci_vec_base];
            MD_FLOAT* ci_v  = &atom->cl_v[ci_vec_base];
            int* ci_t       = &atom->cl_t[ci_sca_base];
            MD_FLOAT bbminx = INFINITY, bbmaxx = -INFINITY;
            MD_FLOAT bbminy = INFINITY, bbmaxy = -INFINITY;
            MD_FLOAT bbminz = INFINITY, bbmaxz = -INFINITY;

            atom->iclusters[ci].natoms = 0;
            for (int cii = 0; cii < CLUSTER_M; cii++) {
                if (ac < c) {
                    int i         = bins[bin * atoms_per_bin + ac];
                    MD_FLOAT xtmp = atom_x(i);
                    MD_FLOAT ytmp = atom_y(i);
                    MD_FLOAT ztmp = atom_z(i);

                    ci_x[CL_X_INDEX_3D(cii)] = xtmp;
                    ci_x[CL_Y_INDEX_3D(cii)] = ytmp;
                    ci_x[CL_Z_INDEX_3D(cii)] = ztmp;
                    ci_v[CL_X_INDEX_3D(cii)] = atom->vx[i];
                    ci_v[CL_Y_INDEX_3D(cii)] = atom->vy[i];
                    ci_v[CL_Z_INDEX_3D(cii)] = atom->vz[i];

                    // TODO: To create the bounding boxes faster, we can use SIMD
                    // operations
                    if (bbminx > xtmp) {
                        bbminx = xtmp;
                    }
                    if (bbmaxx < xtmp) {
                        bbmaxx = xtmp;
                    }
                    if (bbminy > ytmp) {
                        bbminy = ytmp;
                    }
                    if (bbmaxy < ytmp) {
                        bbmaxy = ytmp;
                    }
                    if (bbminz > ztmp) {
                        bbminz = ztmp;
                    }
                    if (bbmaxz < ztmp) {
                        bbmaxz = ztmp;
                    }

                    ci_t[cii] = atom->type[i];
                    atom->iclusters[ci].natoms++;
                } else {
                    ci_x[CL_X_INDEX_3D(cii)] = INFINITY;
                    ci_x[CL_Y_INDEX_3D(cii)] = INFINITY;
                    ci_x[CL_Z_INDEX_3D(cii)] = INFINITY;
                    ci_t[cii]               = 0;
                }

                ac++;
            }

            atom->cluster_bin[ci]      = bin;
            atom->iclusters[ci].bbminx = bbminx;
            atom->iclusters[ci].bbmaxx = bbmaxx;
            atom->iclusters[ci].bbminy = bbminy;
            atom->iclusters[ci].bbmaxy = bbmaxy;
            atom->iclusters[ci].bbminz = bbminz;
            atom->iclusters[ci].bbmaxz = bbmaxz;
            atom->Nclusters_local++;
        }
    }

    DEBUG_MESSAGE("buildClusters end\n");
}

void buildSuperclusters(Atom* atom) {
    DEBUG_MESSAGE("buildSuperclusters start\n");
    atom->Nclusters_local  = 0;

    /* bin local atoms */
    binAtoms(atom);

    for (int bin = 0; bin < mbins; bin++) {
        int c          = bincount[bin];
        int ac         = 0;
        int natoms_sc  = SCLUSTER_SIZE * CLUSTER_M;
        int nsclusters = ((c + natoms_sc - 1) / natoms_sc);

        // Sort atoms in the Z dimension
        sortAtomsByCoord(atom, 2, bin, 0, c - 1);

        for (int scl = 0; scl < nsclusters; scl++) {
            const int sci = atom->Nclusters_local;
            if (sci >= atom->Nclusters_max) {
                growClusters(atom, 1);
            }

            int scl_bin_offset = scl * natoms_sc;
            MD_FLOAT sc_bbminx = INFINITY, sc_bbmaxx = -INFINITY;
            MD_FLOAT sc_bbminy = INFINITY, sc_bbmaxy = -INFINITY;
            MD_FLOAT sc_bbminz = INFINITY, sc_bbmaxz = -INFINITY;
            atom->siclusters[sci].nclusters = 0;

            for (int scl_z = 0; scl_z < SCLUSTER_SIZE_Z; scl_z++) {
                const int scl_z_size = SCLUSTER_SIZE_X * SCLUSTER_SIZE_Y * CLUSTER_M;
                const int atom_CL_Z_OFFSET = scl_bin_offset + scl_z * scl_z_size;
                const int atom_scl_z_end_idx = MIN(atom_CL_Z_OFFSET + scl_z_size - 1, c - 1);

                // Sort atoms in the Y dimension
                sortAtomsByCoord(atom, 1, bin, atom_CL_Z_OFFSET, atom_scl_z_end_idx);

                for (int scl_y = 0; scl_y < SCLUSTER_SIZE_Y; scl_y++) {
                    const int scl_y_size = SCLUSTER_SIZE_X * CLUSTER_M;
                    const int atom_CL_Y_OFFSET = scl_bin_offset +
                                                  scl_z * scl_z_size +
                                                  scl_y * scl_y_size;
                    const int atom_scl_y_end_idx = MIN(atom_CL_Y_OFFSET + scl_y_size - 1, c - 1);

                    // Sort atoms in the X dimension
                    sortAtomsByCoord(atom, 0, bin, atom_CL_Y_OFFSET, atom_scl_y_end_idx);

                    for (int scl_x = 0; scl_x < SCLUSTER_SIZE_X; scl_x++) {
                        const int sci_ci = atom->siclusters[sci].nclusters;
                        const int ci = sci * SCLUSTER_SIZE + sci_ci;
                        int sci_vec_base = SCI_VECTOR3_BASE_INDEX(sci);
                        int sci_sca_base = SCI_SCALAR_BASE_INDEX(sci);
                        MD_FLOAT* sci_x  = &atom->cl_x[SCI_VECTOR_BASE_INDEX(sci)];
                        MD_FLOAT* sci_v  = &atom->cl_v[sci_vec_base];
                        int* sci_t       = &atom->cl_t[sci_sca_base];

                        MD_FLOAT bbminx = INFINITY, bbmaxx = -INFINITY;
                        MD_FLOAT bbminy = INFINITY, bbmaxy = -INFINITY;
                        MD_FLOAT bbminz = INFINITY, bbmaxz = -INFINITY;
                        atom->iclusters[ci].natoms = 0;

                        for (int cii = 0; cii < CLUSTER_M; cii++) {
                            if (ac < c) {
                                int i         = bins[bin * atoms_per_bin + ac];
                                MD_FLOAT xtmp = atom_x(i);
                                MD_FLOAT ytmp = atom_y(i);
                                MD_FLOAT ztmp = atom_z(i);

                                sci_x[CL_X_INDEX(sci_ci * CLUSTER_M + cii)] = xtmp;
                                sci_x[CL_Y_INDEX(sci_ci * CLUSTER_M + cii)] = ytmp;
                                sci_x[CL_Z_INDEX(sci_ci * CLUSTER_M + cii)] = ztmp;
                                sci_v[CL_X_INDEX_3D(sci_ci * CLUSTER_M + cii)] = atom->vx[i];
                                sci_v[CL_Y_INDEX_3D(sci_ci * CLUSTER_M + cii)] = atom->vy[i];
                                sci_v[CL_Z_INDEX_3D(sci_ci * CLUSTER_M + cii)] = atom->vz[i];
                                sci_t[sci_ci * CLUSTER_M + cii] = atom->type[i];

                                // TODO: To create the bounding boxes faster, we can use
                                // SIMD operations
                                if (bbminx > xtmp) {
                                    bbminx = xtmp;
                                }
                                if (bbmaxx < xtmp) {
                                    bbmaxx = xtmp;
                                }
                                if (bbminy > ytmp) {
                                    bbminy = ytmp;
                                }
                                if (bbmaxy < ytmp) {
                                    bbmaxy = ytmp;
                                }
                                if (bbminz > ztmp) {
                                    bbminz = ztmp;
                                }
                                if (bbmaxz < ztmp) {
                                    bbmaxz = ztmp;
                                }

                                atom->iclusters[ci].natoms++;
                            } else {
                                sci_x[CL_X_INDEX(sci_ci * CLUSTER_M + cii)] = INFINITY;
                                sci_x[CL_Y_INDEX(sci_ci * CLUSTER_M + cii)] = INFINITY;
                                sci_x[CL_Z_INDEX(sci_ci * CLUSTER_M + cii)] = INFINITY;
                                sci_t[sci_ci * CLUSTER_M + cii]                = 0;
                            }

                            ac++;
                        }

                        atom->iclusters[ci].bbminx = bbminx;
                        atom->iclusters[ci].bbmaxx = bbmaxx;
                        atom->iclusters[ci].bbminy = bbminy;
                        atom->iclusters[ci].bbmaxy = bbmaxy;
                        atom->iclusters[ci].bbminz = bbminz;
                        atom->iclusters[ci].bbmaxz = bbmaxz;

                        // TODO: To create the bounding boxes faster, we can use SIMD
                        // operations
                        if (sc_bbminx > bbminx) {
                            sc_bbminx = bbminx;
                        }
                        if (sc_bbmaxx < bbmaxx) {
                            sc_bbmaxx = bbmaxx;
                        }
                        if (sc_bbminy > bbminy) {
                            sc_bbminy = bbminy;
                        }
                        if (sc_bbmaxy < bbmaxy) {
                            sc_bbmaxy = bbmaxy;
                        }
                        if (sc_bbminz > bbminz) {
                            sc_bbminz = bbminz;
                        }
                        if (sc_bbmaxz < bbmaxz) {
                            sc_bbmaxz = bbmaxz;
                        }

                        atom->siclusters[sci].nclusters++;
                    }
                }
            }

            atom->cluster_bin[sci]       = bin;
            atom->siclusters[sci].bbminx = sc_bbminx;
            atom->siclusters[sci].bbmaxx = sc_bbmaxx;
            atom->siclusters[sci].bbminy = sc_bbminy;
            atom->siclusters[sci].bbmaxy = sc_bbmaxy;
            atom->siclusters[sci].bbminz = sc_bbminz;
            atom->siclusters[sci].bbmaxz = sc_bbmaxz;
            atom->Nclusters_local++;
        }
    }

    /*
    for(int sci = 0; sci < atom->Nclusters_local; sci++) {
        printf("Supercluster %d: nclusters=%d, bbminx=%.6f, bbmaxx=%.6f, bbminy=%.6f, bbmaxy=%.6f, bbminz=%.6f, bbmaxz=%.6f\n",
            sci,
            atom->siclusters[sci].nclusters,
            atom->siclusters[sci].bbminx,
            atom->siclusters[sci].bbmaxx,
            atom->siclusters[sci].bbminy,
            atom->siclusters[sci].bbmaxy,
            atom->siclusters[sci].bbminz,
            atom->siclusters[sci].bbmaxz);
        int sci_vec_base = SCI_VECTOR_BASE_INDEX(sci);
        MD_FLOAT* sci_x = &atom->cl_x[sci_vec_base];
        for (int sci_ci = 0; sci_ci < atom->siclusters[sci].nclusters; sci_ci++) {
            const int ci = sci * SCLUSTER_SIZE + sci_ci;
            printf("  Cluster %d (in supercluster %d):\n", ci, sci);
            for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
                printf("    Atom %d: x=%.6f y=%.6f z=%.6f\n",
                    cii,
                    sci_x[CL_X_INDEX(sci_ci * CLUSTER_M + cii)],
                    sci_x[CL_Y_INDEX(sci_ci * CLUSTER_M + cii)],
                    sci_x[CL_Z_INDEX(sci_ci * CLUSTER_M + cii)]);
            }
        }
    }
    */

    DEBUG_MESSAGE("buildSuperclusters end\n");
}

void defineJClusters(Parameter* param, Atom* atom) {
    DEBUG_MESSAGE("defineJClusters start\n");

    const int jfac            = MAX(1, CLUSTER_N / CLUSTER_M);
    atom->ncj                 = atom->Nclusters_local * SCLUSTER_SIZE / jfac;

    for (int ci = 0; ci < atom->Nclusters_local * SCLUSTER_SIZE; ci++) {
        int cj0 = CJ0_FROM_CI(ci);

        if (CLUSTER_M == CLUSTER_N) {
            atom->jclusters[cj0].bbminx = atom->iclusters[ci].bbminx;
            atom->jclusters[cj0].bbmaxx = atom->iclusters[ci].bbmaxx;
            atom->jclusters[cj0].bbminy = atom->iclusters[ci].bbminy;
            atom->jclusters[cj0].bbmaxy = atom->iclusters[ci].bbmaxy;
            atom->jclusters[cj0].bbminz = atom->iclusters[ci].bbminz;
            atom->jclusters[cj0].bbmaxz = atom->iclusters[ci].bbmaxz;
            atom->jclusters[cj0].natoms = atom->iclusters[ci].natoms;

        } else if (CLUSTER_M > CLUSTER_N) {
            int cj1         = CJ1_FROM_CI(ci);
            int ci_vec_base = CI_VECTOR_BASE_INDEX(ci);
            MD_FLOAT* ci_x  = &atom->cl_x[ci_vec_base];
            MD_FLOAT bbminx = INFINITY, bbmaxx = -INFINITY;
            MD_FLOAT bbminy = INFINITY, bbmaxy = -INFINITY;
            MD_FLOAT bbminz = INFINITY, bbmaxz = -INFINITY;

            for (int cii = 0; cii < MIN(atom->iclusters[ci].natoms, CLUSTER_N); cii++) {
                MD_FLOAT xtmp = ci_x[CL_X_INDEX(cii)];
                MD_FLOAT ytmp = ci_x[CL_Y_INDEX(cii)];
                MD_FLOAT ztmp = ci_x[CL_Z_INDEX(cii)];

                // TODO: To create the bounding boxes faster, we can use SIMD operations
                if (bbminx > xtmp) {
                    bbminx = xtmp;
                }
                if (bbmaxx < xtmp) {
                    bbmaxx = xtmp;
                }
                if (bbminy > ytmp) {
                    bbminy = ytmp;
                }
                if (bbmaxy < ytmp) {
                    bbmaxy = ytmp;
                }
                if (bbminz > ztmp) {
                    bbminz = ztmp;
                }
                if (bbmaxz < ztmp) {
                    bbmaxz = ztmp;
                }
            }

            atom->jclusters[cj0].bbminx = bbminx;
            atom->jclusters[cj0].bbmaxx = bbmaxx;
            atom->jclusters[cj0].bbminy = bbminy;
            atom->jclusters[cj0].bbmaxy = bbmaxy;
            atom->jclusters[cj0].bbminz = bbminz;
            atom->jclusters[cj0].bbmaxz = bbmaxz;
            atom->jclusters[cj0].natoms = MIN(atom->iclusters[ci].natoms, CLUSTER_N);

            bbminx = INFINITY, bbmaxx = -INFINITY;
            bbminy = INFINITY, bbmaxy = -INFINITY;
            bbminz = INFINITY, bbmaxz = -INFINITY;

            for (int cii = CLUSTER_N; cii < atom->iclusters[ci].natoms; cii++) {
                MD_FLOAT xtmp = ci_x[CL_X_INDEX(cii)];
                MD_FLOAT ytmp = ci_x[CL_Y_INDEX(cii)];
                MD_FLOAT ztmp = ci_x[CL_Z_INDEX(cii)];

                // TODO: To create the bounding boxes faster, we can use SIMD operations
                if (bbminx > xtmp) {
                    bbminx = xtmp;
                }
                if (bbmaxx < xtmp) {
                    bbmaxx = xtmp;
                }
                if (bbminy > ytmp) {
                    bbminy = ytmp;
                }
                if (bbmaxy < ytmp) {
                    bbmaxy = ytmp;
                }
                if (bbminz > ztmp) {
                    bbminz = ztmp;
                }
                if (bbmaxz < ztmp) {
                    bbmaxz = ztmp;
                }
            }

            atom->jclusters[cj1].bbminx = bbminx;
            atom->jclusters[cj1].bbmaxx = bbmaxx;
            atom->jclusters[cj1].bbminy = bbminy;
            atom->jclusters[cj1].bbmaxy = bbmaxy;
            atom->jclusters[cj1].bbminz = bbminz;
            atom->jclusters[cj1].bbmaxz = bbmaxz;
            atom->jclusters[cj1].natoms = MAX(0, atom->iclusters[ci].natoms - CLUSTER_N);
        } else {
            if (ci % 2 == 0) {
                const int ci1               = ci + 1;
                atom->jclusters[cj0].bbminx = MIN(atom->iclusters[ci].bbminx,
                    atom->iclusters[ci1].bbminx);
                atom->jclusters[cj0].bbmaxx = MAX(atom->iclusters[ci].bbmaxx,
                    atom->iclusters[ci1].bbmaxx);
                atom->jclusters[cj0].bbminy = MIN(atom->iclusters[ci].bbminy,
                    atom->iclusters[ci1].bbminy);
                atom->jclusters[cj0].bbmaxy = MAX(atom->iclusters[ci].bbmaxy,
                    atom->iclusters[ci1].bbmaxy);
                atom->jclusters[cj0].bbminz = MIN(atom->iclusters[ci].bbminz,
                    atom->iclusters[ci1].bbminz);
                atom->jclusters[cj0].bbmaxz = MAX(atom->iclusters[ci].bbmaxz,
                    atom->iclusters[ci1].bbmaxz);
                atom->jclusters[cj0].natoms = atom->iclusters[ci].natoms +
                                              atom->iclusters[ci1].natoms;
            }
        }
    }

    DEBUG_MESSAGE("defineJClusters end\n");
}

void binJClusters(Parameter* param, Atom* atom) {
    DEBUG_MESSAGE("binJClusters start\n");

    /*
    DEBUG_MESSAGE("Nghost = %d\n", atom->Nclusters_ghost);
    for(int ci = atom->Nclusters_local; ci < atom->Nclusters_local + 4; ci++) {
        MD_FLOAT *cptr = cluster_pos_ptr(ci);
        DEBUG_MESSAGE("Cluster %d:\n", ci);
        DEBUG_MESSAGE("bin=%d, Natoms=%d, bbox={%f,%f},{%f,%f},{%f,%f}\n",
            atom->cluster_bin[ci],
            atom->clusters[ci].natoms,
            atom->clusters[ci].bbminx,
            atom->clusters[ci].bbmaxx,
            atom->clusters[ci].bbminy,
            atom->clusters[ci].bbmaxy,
            atom->clusters[ci].bbminz,
            atom->clusters[ci].bbmaxz);

        for(int cii = 0; cii < CLUSTER_M; cii++) {
            DEBUG_MESSAGE("%f, %f, %f\n", cluster_x(cptr, cii), cluster_y(cptr, cii),
    cluster_z(cptr, cii));
        }
    }
    */

    const int nlocal = atom->Nclusters_local;
    const int ncj    = get_ncj_from_nci(nlocal);
    int resize       = 1;
    while (resize > 0) {
        resize = 0;

        for (int bin = 0; bin < mbins; bin++) {
            bin_nclusters[bin] = 0;
        }

        for (int ci = 0; ci < nlocal && !resize; ci++) {
            if(param->super_clustering) {
                for(int sci_cj = 0; sci_cj < SCLUSTER_SIZE; sci_cj++) {
                    int bin = atom->cluster_bin[ci];
                    int c   = bin_nclusters[bin];
                    if (c + 1 < clusters_per_bin) {
                        bin_clusters[bin * clusters_per_bin + c] = ci * SCLUSTER_SIZE + sci_cj;
                        bin_nclusters[bin]++;
                    } else {
                        resize = 1;
                    }
                }
            } else {
                // Assure we add this j-cluster only once in the bin
                if (CLUSTER_M >= CLUSTER_N || ci % 2 == 0) {
                    int bin = atom->cluster_bin[ci];
                    int c   = bin_nclusters[bin];
                    if (c + 1 < clusters_per_bin) {
                        bin_clusters[bin * clusters_per_bin + c] = CJ0_FROM_CI(ci);
                        bin_nclusters[bin]++;

                        if (CLUSTER_M > CLUSTER_N) {
                            int cj1 = CJ1_FROM_CI(ci);
                            if (atom->jclusters[cj1].natoms > 0) {
                                bin_clusters[bin * clusters_per_bin + c + 1] = cj1;
                                bin_nclusters[bin]++;
                            }
                        }
                    } else {
                        resize = 1;
                    }
                }
            }
        }

        for (int cg = 0; cg < atom->Nclusters_ghost && !resize; cg++) {
            const int cj = ncj + cg;
            int ix = -1, iy = -1;
            MD_FLOAT xtmp, ytmp;

            if (shellMethod == halfShell && !halfZoneCluster(atom, cj)) {
                continue;
            }

            if (atom->jclusters[cj].natoms > 0) {
                int cj_vec_base  = CJ_VECTOR_BASE_INDEX(cj);
                MD_FLOAT* cj_x   = &atom->cl_x[cj_vec_base];
                MD_FLOAT cj_minz = atom->jclusters[cj].bbminz;

                xtmp = cj_x[CL_X_INDEX(0)];
                ytmp = cj_x[CL_Y_INDEX(0)];

                coord2bin2D(xtmp, ytmp, &ix, &iy);
                ix = MAX(MIN(ix, mbinx - 1), 0);
                iy = MAX(MIN(iy, mbiny - 1), 0);

                for (int cjj = 1; cjj < atom->jclusters[cj].natoms; cjj++) {
                    int nix, niy;
                    xtmp = cj_x[CL_X_INDEX(cjj)];
                    ytmp = cj_x[CL_Y_INDEX(cjj)];

                    coord2bin2D(xtmp, ytmp, &nix, &niy);
                    nix = MAX(MIN(nix, mbinx - 1), 0);
                    niy = MAX(MIN(niy, mbiny - 1), 0);

                    // Always put the cluster on the bin of its innermost atom so
                    // the cluster should be closer to local clusters
                    if (atom->PBCx[cg] > 0 && ix > nix) {
                        ix = nix;
                    }
                    if (atom->PBCx[cg] < 0 && ix < nix) {
                        ix = nix;
                    }
                    if (atom->PBCy[cg] > 0 && iy > niy) {
                        iy = niy;
                    }
                    if (atom->PBCy[cg] < 0 && iy < niy) {
                        iy = niy;
                    }
                }

                int bin = iy * mbinx + ix + 1;
                int c   = bin_nclusters[bin];
                if (c < clusters_per_bin) {
                    // Insert the current ghost cluster in the bin keeping clusters
                    // sorted by z coordinate
                    int inserted = 0;
                    for (int i = 0; i < c; i++) {
                        int last_cl = bin_clusters[bin * clusters_per_bin + i];
                        if (atom->jclusters[last_cl].bbminz > cj_minz) {
                            bin_clusters[bin * clusters_per_bin + i] = cj;

                            for (int j = i + 1; j <= c; j++) {
                                int tmp = bin_clusters[bin * clusters_per_bin + j];
                                bin_clusters[bin * clusters_per_bin + j] = last_cl;
                                last_cl                                  = tmp;
                            }

                            inserted = 1;
                            break;
                        }
                    }

                    if (!inserted) {
                        bin_clusters[bin * clusters_per_bin + c] = cj;
                    }

                    bin_nclusters[bin]++;
                } else {
                    resize = 1;
                }
            }
        }

        if (resize) {
            free(bin_clusters);
            clusters_per_bin *= 2;
            bin_clusters = (int*)malloc(mbins * clusters_per_bin * sizeof(int));
        }
    }

    /*
    DEBUG_MESSAGE("bin_nclusters\n");
    for(int i = 0; i < mbins; i++) { DEBUG_MESSAGE("%d, ", bin_nclusters[i]); }
    DEBUG_MESSAGE("\n");
    */

    DEBUG_MESSAGE("binJClusters stop\n");
}

void updateSingleAtoms(Parameter* param, Atom* atom) {
    DEBUG_MESSAGE("updateSingleAtoms start\n");
    int Natom = 0;

    if(param->super_clustering) {
        for (int sci = 0; sci < atom->Nclusters_local; sci++) {
            int sci_vec_base = SCI_VECTOR3_BASE_INDEX(sci);
            int sci_sca_base = SCI_SCALAR_BASE_INDEX(sci);
            MD_FLOAT* sci_x  = &atom->cl_x[SCI_VECTOR_BASE_INDEX(sci)];
            MD_FLOAT* sci_v  = &atom->cl_v[sci_vec_base];
            int* sci_t       = &atom->cl_t[sci_sca_base];

            for(int sci_ci = 0; sci_ci < atom->siclusters[sci].nclusters; sci_ci++) {
                const int ci = sci * SCLUSTER_SIZE + sci_ci;
                for(int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
                    atom_x(Natom)     = sci_x[CL_X_INDEX(sci_ci * CLUSTER_M + cii)];
                    atom_y(Natom)     = sci_x[CL_Y_INDEX(sci_ci * CLUSTER_M + cii)];
                    atom_z(Natom)     = sci_x[CL_Z_INDEX(sci_ci * CLUSTER_M + cii)];
                    atom_vx(Natom)    = sci_v[CL_X_INDEX_3D(sci_ci * CLUSTER_M + cii)];
                    atom_vy(Natom)    = sci_v[CL_Y_INDEX_3D(sci_ci * CLUSTER_M + cii)];
                    atom_vz(Natom)    = sci_v[CL_Z_INDEX_3D(sci_ci * CLUSTER_M + cii)];
                    atom->type[Natom] = sci_t[sci_ci * CLUSTER_M + cii];
                    Natom++;
                }
            }
        }
    } else {
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
            int ci_sca_base = CI_SCALAR_BASE_INDEX(ci);
            MD_FLOAT* ci_x  = &atom->cl_x[CI_VECTOR_BASE_INDEX(ci)];
            MD_FLOAT* ci_v  = &atom->cl_v[ci_vec_base];
            int* ci_t       = &atom->cl_t[ci_sca_base];

            for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
                atom_x(Natom)     = ci_x[CL_X_INDEX(cii)];
                atom_y(Natom)     = ci_x[CL_Y_INDEX(cii)];
                atom_z(Natom)     = ci_x[CL_Z_INDEX(cii)];
                atom->vx[Natom]   = ci_v[CL_X_INDEX_3D(cii)];
                atom->vy[Natom]   = ci_v[CL_Y_INDEX_3D(cii)];
                atom->vz[Natom]   = ci_v[CL_Z_INDEX_3D(cii)];
                atom->type[Natom] = ci_t[cii];
                Natom++;
            }
        }
    }

    if (Natom != atom->Nlocal) {
        fprintf(stderr, "updateSingleAtoms(): Number of atoms changed!\n");
    }

    DEBUG_MESSAGE("updateSingleAtoms stop\n");
}

// MPI Shell Methods
static int eightZoneCluster(Atom* atom, int cj) {
    // Mapping: 0->0, 1->1, 2->2, 3->6, 4->3, 5->5, 6->4, 7->7
    int zoneMapping[] = { 0, 1, 2, 6, 3, 5, 4, 7 };
    int zone          = 0;
    MD_FLOAT* hi      = atom->mybox.hi;

    if (atom->jclusters[cj].bbminx >= hi[0]) {
        zone += 1;
    }

    if (atom->jclusters[cj].bbminy >= hi[1]) {
        zone += 2;
    }

    if (atom->jclusters[cj].bbminz >= hi[2]) {
        zone += 4;
    }

    return zoneMapping[zone];
}

static int halfZoneCluster(Atom* atom, int cj) {
    MD_FLOAT* hi = atom->mybox.hi;
    MD_FLOAT* lo = atom->mybox.lo;

    if (atom->jclusters[cj].bbmaxx < lo[0] && atom->jclusters[cj].bbmaxy < hi[1] &&
        atom->jclusters[cj].bbmaxz < hi[2]) {
        return 0;
    } else if (atom->jclusters[cj].bbmaxy < lo[1] &&
               atom->jclusters[cj].bbmaxz < hi[2]) {
        return 0;
    } else if (atom->jclusters[cj].bbmaxz < lo[2]) {
        return 0;
    } else {
        return 1;
    }
}

int BoxGhostDistance(Atom* atom, int ci, int cj) {
    MD_FLOAT dl  = atom->jclusters[ci].bbminx - atom->jclusters[cj].bbmaxx;
    MD_FLOAT dh  = atom->jclusters[cj].bbminx - atom->jclusters[ci].bbmaxx;
    MD_FLOAT dm  = MAX(dl, dh);
    MD_FLOAT dm0 = MAX(dm, 0.0);
    MD_FLOAT dx2 = dm0 * dm0;

    dl           = atom->jclusters[ci].bbminy - atom->jclusters[cj].bbmaxy;
    dh           = atom->jclusters[cj].bbminy - atom->jclusters[ci].bbmaxy;
    dm           = MAX(dl, dh);
    dm0          = MAX(dm, 0.0);
    MD_FLOAT dy2 = dm0 * dm0;

    dl           = atom->jclusters[ci].bbminz - atom->jclusters[cj].bbmaxz;
    dh           = atom->jclusters[cj].bbminz - atom->jclusters[ci].bbmaxz;
    dm           = MAX(dl, dh);
    dm0          = MAX(dm, 0.0);
    MD_FLOAT dz2 = dm0 * dm0;

    return dx2 > cutneighsq ? 0 : dy2 > cutneighsq ? 0 : dz2 > cutneighsq ? 0 : 1;
}

static int ghostClusterinRange(Atom* atom, int cs, int cg, MD_FLOAT rsq) {
    int cs_vec_base = CJ_VECTOR_BASE_INDEX(cs);
    int cj_vec_base = CJ_VECTOR_BASE_INDEX(cg);
    MD_FLOAT* cs_x  = &atom->cl_x[cs_vec_base];
    MD_FLOAT* cg_x  = &atom->cl_x[cj_vec_base];

    for (int cii = 0; cii < atom->jclusters[cs].natoms; cii++) {
        for (int cjj = 0; cjj < atom->jclusters[cg].natoms; cjj++) {
            MD_FLOAT delx = cs_x[CL_X_INDEX(cii)] - cg_x[CL_X_INDEX(cjj)];
            MD_FLOAT dely = cs_x[CL_Y_INDEX(cii)] - cg_x[CL_Y_INDEX(cjj)];
            MD_FLOAT delz = cs_x[CL_Z_INDEX(cii)] - cg_x[CL_Z_INDEX(cjj)];
            if (delx * delx + dely * dely + delz * delz < rsq) {
                return 1;
            }
        }
    }

    return 0;
}

static void neighborGhost(Atom* atom, Neighbor* neighbor){
    int Nshell         = 0;
    int Ncluster_local = atom->Nclusters_local;
    int Nclusterghost  = atom->Nclusters_ghost;
    if (neighbor->listshell) free(neighbor->listshell);
    neighbor->listshell = (int*)malloc(Nclusterghost * sizeof(int));
    int* listzone       = (int*)malloc(8 * Nclusterghost * sizeof(int));
    int countCluster[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

    // Selecting ghost atoms for interaction and putting them into regions
    for (int cg = atom->ncj; cg < atom->ncj + Nclusterghost; cg++) {
        int czone = eightZoneCluster(atom, cg);
        int* list = &listzone[Nclusterghost * czone];
        int n     = countCluster[czone];
        list[n]   = cg;
        countCluster[czone]++;
        // It is only necessary to find neighbour particles for 3 regions
        // if(czone == 1 || czone == 2 || czone == 3)
        // neighbor->listshell[Nshell++] = cg;
    }

    for (int zone = 1; zone <= 3; zone++) {
        int* list = &listzone[Nclusterghost * zone];
        for (int n = 0; n < countCluster[zone]; n++)
            neighbor->listshell[Nshell++] = list[n];
    }

    neighbor->Nshell = Nshell;
    if (neighbor->numNeighShell) free(neighbor->numNeighShell);
    if (neighbor->neighshell) free(neighbor->neighshell);
    neighbor->neighshell    = (int*)malloc(Nshell * neighbor->maxneighs * sizeof(int));
    neighbor->numNeighShell = (int*)malloc(Nshell * sizeof(int));

    int resize = 1;

    while (resize) {
        resize = 0;
        for (int ic = 0; ic < Nshell; ic++) {
            int* neighshell = &(neighbor->neighshell[ic * neighbor->maxneighs]);
            int n           = 0;
            int icluster    = neighbor->listshell[ic];
            int iczone      = eightZoneCluster(atom, icluster);

            for (int jczone = 0; jczone < 8; jczone++) {

                if (jczone <= iczone) continue;
                if (iczone == 1 && (jczone == 5 || jczone == 6 || jczone == 7)) continue;
                if (iczone == 2 && (jczone == 4 || jczone == 6 || jczone == 7)) continue;
                if (iczone == 3 && (jczone == 4 || jczone == 5 || jczone == 7)) continue;

                int Ncluster  = countCluster[jczone];
                int* loc_zone = &listzone[jczone * Nclusterghost];

                for (int k = 0; k < Ncluster; k++) {
                    int jcluster = loc_zone[k];

                    if (BoxGhostDistance(atom, icluster, jcluster)) {
                        if (ghostClusterinRange(atom, icluster, jcluster, cutneighsq))
                            neighshell[n++] = jcluster;
                    }
                }
            }
            neighbor->numNeighShell[ic] = n;

            if (n >= neighbor->maxneighs) {
                resize              = 1;
                neighbor->maxneighs = n * 1.2;
                fprintf(stdout,
                    "RESIZE EIGHT SHELL %d, PROC %d\n",
                    neighbor->maxneighs,
                    me);
                break;
            }
        }

        if (resize) {
            free(neighbor->neighshell);
            neighbor->neighshell = (int*)malloc(
                Nshell * neighbor->maxneighs * sizeof(int));
        }
    }

    free(listzone);
}
