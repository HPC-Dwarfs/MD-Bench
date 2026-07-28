/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <allocate.h>
#include <atom.h>
#include <neighbor.h>
#include <parameter.h>
#include <util.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// USE_SIMD_NEIGHBOR controls the vectorized neighbor-list build independently
// of USE_SIMD_KERNEL (which only affects the force kernel). Under AOS/SOA it
// writes surviving neighbor indices through the neighs() macro directly into
// the padded per-atom slot; under CSR it appends them to that atom's
// growable per-thread buffer instead (see buildNeighborCPU) -- the
// vectorized bin-scan itself (gather candidate positions, build a lane mask,
// scalar bit-scan over survivors) is identical either way, only the storage
// step differs.
#ifdef __SIMD_NEIGHBOR__
#include <simd.h>
#if !defined(__ISA_AVX512__) && !defined(__ISA_AVX2__) && !defined(__ISA_NEON__) \
    && !defined(__ISA_SVE__) && !defined(__ISA_SVE2__)
#error "USE_SIMD_NEIGHBOR requires SIMD=AVX2/AVX512/NEON/SVE/SVE2 (needs " \
    "simd_mask_not/simd_mask_i32_cond_eq, not implemented for SSE/AVX so far)"
#endif
#endif

#define SMALL  1.0e-6
#define FACTOR 0.999

#ifdef CUDA_TARGET
BuildNeighborFunction buildNeighbor = buildNeighborCUDA;
#else
BuildNeighborFunction buildNeighbor = buildNeighborCPU;
#endif
PruneNeighborFunction pruneNeighbor = pruneNeighborCPU;

MD_FLOAT xprd, yprd, zprd;
MD_FLOAT bininvx, bininvy, bininvz;
int mbinxlo, mbinylo, mbinzlo; // not used
int nbinx, nbiny, nbinz;
int mbinx, mbiny, mbinz; // n bins in x, y, z
int* bincount;
int* bins;
int mbins;         // total number of bins
int atoms_per_bin; // max atoms per bin
MD_FLOAT cutneigh;
MD_FLOAT cutneighsq;        // neighbor (outer) cutoff squared
MD_FLOAT cutneigh_inner_sq; // inner cutoff squared (cutforce + skin)
int dcut_enabled;           // double-cutoff pruning active for this run
int* is_inner_buf;          // reusable scratch for pruneNeighborCPU, sized by maxneighs
int nmax;
int nstencil; // # of bins in stencil
int* stencil; // stencil list of bin offsets
MD_FLOAT binsizex, binsizey, binsizez;
static int coord2bin(MD_FLOAT, MD_FLOAT, MD_FLOAT);
static MD_FLOAT bindist(int, int, int);

// MPI version
int pad_x, pad_y, pad_z; // version MPI
int me;                  // rank
int method;              // method
int half_stencil;        // If half stencil exist
int shellMethod;         // If shell method exist
static int ghostZone(Atom*, int);
static int eightZone(Atom*, int);
static int halfZone(Atom*, int);
static void neighborGhost(Atom*, Neighbor*);
static inline int skipNeigh(Atom* atom, int i, int j);

/* exported subroutines */
void initNeighbor(Neighbor* neighbor, Parameter* param)
{
    MD_FLOAT neighscale      = 5.0 / 6.0;
    xprd                     = param->nx * param->lattice;
    yprd                     = param->ny * param->lattice;
    zprd                     = param->nz * param->lattice;
    cutneigh                 = param->cutneigh;
    nbinx                    = MAX(1, neighscale * param->nx);
    nbiny                    = MAX(1, neighscale * param->ny);
    nbinz                    = MAX(1, neighscale * param->nz);
    nmax                     = 0;
    atoms_per_bin            = 8;
    stencil                  = NULL;
    bins                     = NULL;
    bincount                 = NULL;
    neighbor->maxneighs      = 100;
    neighbor->numneigh       = NULL;
    neighbor->numneigh_inner = NULL;
    neighbor->neighbors      = NULL;
    neighbor->neigh_start    = NULL;
    is_inner_buf = (int*)allocate(ALIGNMENT, neighbor->maxneighs * sizeof(int));
    method = param->method;
    if (method == halfShell || method == eightShell) {
        param->half_neigh = 1;
        shellMethod       = 1;
    }

    if (method == halfStencil) {
        param->half_neigh = 0;
        half_stencil      = 1;
    }
    neighbor->half_neigh = param->half_neigh;

    me = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
    neighbor->Nshell        = 0;
    neighbor->numNeighShell = NULL;
    neighbor->neighshell    = NULL;
    neighbor->listshell     = NULL;
}

void setupNeighbor(Parameter* param)
{
    MD_FLOAT coord;
    int mbinxhi, mbinyhi, mbinzhi;
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

    cutneighsq               = cutneigh * cutneigh;
    const MD_FLOAT cut_inner = param->cutforce + param->skin;
    cutneigh_inner_sq        = cut_inner * cut_inner;
    // Only prune when an outer skin actually widens the list beyond the inner
    // cutoff; this also keeps EAM (which sets its own cutneigh) on the full list.
    dcut_enabled = (param->outer_skin > 0.0) && (cutneigh_inner_sq < cutneighsq);

    if (param->input_file != NULL) {
        binsizex = cutneigh * 0.5;
        binsizey = cutneigh * 0.5;
        binsizez = cutneigh * 0.5;
        nbinx    = (int)((param->xhi - param->xlo) / binsizex);
        nbiny    = (int)((param->yhi - param->ylo) / binsizey);
        nbinz    = (int)((param->zhi - param->zlo) / binsizez);
        if (nbinx == 0) {
            nbinx = 1;
        }
        if (nbiny == 0) {
            nbiny = 1;
        }
        if (nbinz == 0) {
            nbinz = 1;
        }
        bininvx = nbinx / (param->xhi - param->xlo);
        bininvy = nbiny / (param->yhi - param->ylo);
        bininvz = nbinz / (param->zhi - param->zlo);
    } else {
        binsizex = xprd / nbinx;
        binsizey = yprd / nbiny;
        binsizez = zprd / nbinz;
        bininvx  = 1.0 / binsizex;
        bininvy  = 1.0 / binsizey;
        bininvz  = 1.0 / binsizez;
    }

    pad_x = (int)(cutneigh * bininvx);
    while (pad_x * binsizex < FACTOR * cutneigh)
        pad_x++;
    pad_y = (int)(cutneigh * bininvy);
    while (pad_y * binsizey < FACTOR * cutneigh)
        pad_y++;
    pad_z = (int)(cutneigh * bininvz);
    while (pad_z * binsizez < FACTOR * cutneigh)
        pad_z++;

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

    coord   = zlo - cutneigh - SMALL * zprd;
    mbinzlo = (int)(coord * bininvz);
    if (coord < 0.0) {
        mbinzlo = mbinzlo - 1;
    }
    coord   = zhi + cutneigh + SMALL * zprd;
    mbinzhi = (int)(coord * bininvz);

    mbinxlo = mbinxlo - 1;
    mbinxhi = mbinxhi + 1;
    mbinylo = mbinylo - 1;
    mbinyhi = mbinyhi + 1;

    mbinzlo = mbinzlo - 1;
    mbinzhi = mbinzhi + 1;

    /*
        mbinxlo = mbinxlo - 1;
        mbinxhi = mbinxhi + 1;
        mbinx   = mbinxhi - mbinxlo + 1;

        mbinylo = mbinylo - 1;
        mbinyhi = mbinyhi + 1;
        mbiny   = mbinyhi - mbinylo + 1;

        mbinzlo = mbinzlo - 1;
        mbinzhi = mbinzhi + 1;
        mbinz   = mbinzhi - mbinzlo + 1;
    */
    nextx = (int)(cutneigh * bininvx);
    while (nextx * binsizex < FACTOR * cutneigh) {
        nextx++;
        pad_x++;
    }

    nexty = (int)(cutneigh * bininvy);
    while (nexty * binsizey < FACTOR * cutneigh) {
        nexty++;
        pad_y++;
    }

    nextz = (int)(cutneigh * bininvz);
    while (nextz * binsizez < FACTOR * cutneigh) {
        nextz++;
        pad_z++;
    }

    mbinx = MAX(1, nbinx + 4 * pad_x);
    mbiny = MAX(1, nbiny + 4 * pad_y);
    mbinz = MAX(1, nbinz + 4 * pad_z);

    if (stencil) {
        free(stencil);
    }
    stencil    = (int*)allocate(ALIGNMENT,
        (2 * nextz + 1) * (2 * nexty + 1) * (2 * nextx + 1) * sizeof(int));
    nstencil   = 0;
    int kstart = -nextz;
    int ibin   = 0;
    for (int k = kstart; k <= nextz; k++) {
        for (int j = -nexty; j <= nexty; j++) {
            for (int i = -nextx; i <= nextx; i++) {
                if (bindist(i, j, k) < cutneighsq) {
                    int jbin = k * mbiny * mbinx + j * mbinx + i;
                    if (ibin > jbin && half_stencil) continue;
                    stencil[nstencil++] = jbin;
                    // stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
                }
            }
        }
    }

    mbins = mbinx * mbiny * mbinz;
    if (bincount) {
        free(bincount);
    }
    bincount = (int*)allocate(ALIGNMENT, mbins * sizeof(int));
    if (bins) {
        free(bins);
    }
    bins = (int*)allocate(ALIGNMENT, mbins * atoms_per_bin * sizeof(int));
}

static void buildNeighborLists(Atom* atom, Neighbor* neighbor);

void buildNeighborCPU(Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("buildNeighborCPU begin\n");
    int nall = atom->Nlocal + atom->Nghost;

    if (nall > nmax) {
        nmax = nall;
        if (neighbor->numneigh) free(neighbor->numneigh);
        if (neighbor->numneigh_inner) free(neighbor->numneigh_inner);
        neighbor->numneigh       = (int*)allocate(ALIGNMENT, nmax * sizeof(int));
        neighbor->numneigh_inner = (int*)allocate(ALIGNMENT, nmax * sizeof(int));
#ifdef NBLIST_CSR
        if (neighbor->neigh_start) free(neighbor->neigh_start);
        neighbor->neigh_start = (int*)allocate(ALIGNMENT, (nmax + 1) * sizeof(int));
#else
        if (neighbor->neighbors) free(neighbor->neighbors);
        neighbor->neighbors = (int*)allocate(ALIGNMENT,
            nmax * neighbor->maxneighs * sizeof(int));
        memset(neighbor->neighbors, 0, nmax * neighbor->maxneighs * sizeof(int));
#endif
    }

    binatoms(atom);
    buildNeighborLists(atom, neighbor);

    if (method == eightShell) {
        neighborGhost(atom, neighbor);
    }

    pruneNeighborCPU(NULL, atom, neighbor);
    DEBUG_MESSAGE("buildNeighborCPU end\n");
}

#if !defined(NBLIST_CSR) && !defined(__SIMD_NEIGHBOR__)

/* AOS/SOA, scalar bin scan. Layout only affects the neighs() macro used to
 * store each surviving index, so this one function covers both. */
static void buildNeighborLists(Atom* atom, Neighbor* neighbor)
{
    int resize = 1;

    while (resize) {
        int new_maxneighs = neighbor->maxneighs;
        int resize_local  = 0;

#pragma omp parallel for schedule(runtime) reduction(max                                 \
                                                     : new_maxneighs)                    \
    reduction(|                                                                          \
              : resize_local)
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
                        neighs(neighbor->neighbors, i, n, atom->Nlocal, neighbor) = j;
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
            memset(neighbor->neighbors, 0, atom->Nmax * neighbor->maxneighs * sizeof(int));
            is_inner_buf = (int*)allocate(ALIGNMENT, neighbor->maxneighs * sizeof(int));
        }
    }
}

#elif !defined(NBLIST_CSR) && defined(__SIMD_NEIGHBOR__)

/* AOS/SOA, vectorized bin scan: gather candidate positions, build a lane
 * mask from (i != j) & (half_neigh ? j >= i : true) & (rsq < cutoff), then
 * compact survivors via a scalar bit-scan. */
static void buildNeighborLists(Atom* atom, Neighbor* neighbor)
{
    int resize = 1;

    while (resize) {
        int new_maxneighs = neighbor->maxneighs;
        int resize_local  = 0;

#pragma omp parallel for schedule(runtime) reduction(max                                 \
                                                     : new_maxneighs)                    \
    reduction(|                                                                          \
              : resize_local)
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
                            neighs(neighbor->neighbors, i, n, atom->Nlocal, neighbor) = j;
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
                    MD_SIMD_FLOAT xj = simd_real_gather(j3, &(atom->x[0]), sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT yj = simd_real_gather(j3, &(atom->x[1]), sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT zj = simd_real_gather(j3, &(atom->x[2]), sizeof(MD_FLOAT));
#else
                    MD_SIMD_FLOAT xj = simd_real_gather(jvec, atom->x, sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT yj = simd_real_gather(jvec, atom->y, sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT zj = simd_real_gather(jvec, atom->z, sizeof(MD_FLOAT));
#endif
                    MD_SIMD_FLOAT delx = simd_real_sub(xtmp_vec, xj);
                    MD_SIMD_FLOAT dely = simd_real_sub(ytmp_vec, yj);
                    MD_SIMD_FLOAT delz = simd_real_sub(ztmp_vec, zj);
                    MD_SIMD_FLOAT rsq  = simd_real_fma(
                        delx, delx, simd_real_fma(dely, dely, simd_real_mul(delz, delz)));

                    MD_SIMD_MASK keep = simd_mask_and(
                        bound_mask, simd_mask_not(simd_mask_i32_cond_eq(jvec, ivec)));
                    if (neighbor->half_neigh) {
                        keep = simd_mask_and(
                            keep, simd_mask_not(simd_mask_i32_cond_lt(jvec, ivec)));
                    }

#if LJ_COMB_RULE == LJ_COMB_SINGLE
                    MD_SIMD_FLOAT cutoff_vec = cutoff_vec_single;
#else
                    MD_SIMD_INT type_j_vec = simd_i32_gather(jvec, atom->type, sizeof(int));
                    MD_SIMD_INT tij        = simd_i32_add(tbase_i, type_j_vec);
                    MD_SIMD_FLOAT cutoff_vec =
                        simd_real_gather(tij, atom->cutneighsq, sizeof(MD_FLOAT));
#endif
                    keep = simd_mask_and(keep, simd_mask_cond_lt(rsq, cutoff_vec));

                    unsigned int bits = simd_mask_to_u32(keep);
                    if (bits) {
                        int tmp_j[VECTOR_WIDTH] __attribute__((aligned(64)));
                        simd_i32_store(tmp_j, jvec);
                        while (bits) {
                            int lane = __builtin_ctz(bits);
                            bits &= bits - 1;
                            neighs(neighbor->neighbors, i, n, atom->Nlocal, neighbor) =
                                tmp_j[lane];
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
            memset(neighbor->neighbors, 0, atom->Nmax * neighbor->maxneighs * sizeof(int));
            is_inner_buf = (int*)allocate(ALIGNMENT, neighbor->maxneighs * sizeof(int));
        }
    }
}

#elif defined(NBLIST_CSR) && !defined(__SIMD_NEIGHBOR__)

/* CSR, scalar bin scan: each thread fills its own growable buffer, then a
 * prefix sum over per-thread totals produces global offsets and buffers are
 * merged into the final compact array. */
static void buildNeighborLists(Atom* atom, Neighbor* neighbor)
{
    int Nlocal          = atom->Nlocal;
    int max_threads     = 1;
    int actual_nthreads = 1;
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#endif
    int** t_buf    = (int**)malloc(max_threads * sizeof(int*));
    int*  t_cap    = (int*)malloc(max_threads * sizeof(int));
    int*  t_pos    = (int*)malloc(max_threads * sizeof(int));
    int*  t_owner  = (int*)allocate(ALIGNMENT, Nlocal * sizeof(int));
    int   init_cap = MAX(1, Nlocal / max_threads) * MAX(1, neighbor->maxneighs);
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

#define CSR_APPEND(jval)                                                       \
    do {                                                                       \
        if (t_pos[tid] + n >= t_cap[tid]) {                                    \
            t_cap[tid] *= 2;                                                   \
            t_buf[tid]  = (int*)realloc(t_buf[tid], t_cap[tid] * sizeof(int)); \
        }                                                                      \
        t_buf[tid][t_pos[tid] + n] = (jval);                                   \
        n++;                                                                   \
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
                    int type_j        = atom->type[j];
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
    t_offsets[0] = 0;
    for (int t = 0; t < actual_nthreads; t++)
        t_offsets[t + 1] = t_offsets[t] + t_pos[t];
    int total = t_offsets[actual_nthreads];

    if (neighbor->neighbors) free(neighbor->neighbors);
    neighbor->neighbors = (int*)allocate(
        ALIGNMENT, (MAX(1, total) + VECTOR_WIDTH) * sizeof(int));
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
        is_inner_buf = (int*)allocate(ALIGNMENT, MAX(1, neighbor->maxneighs) * sizeof(int));
    }

    free(t_owner);
    free(t_buf);
    free(t_cap);
    free(t_pos);
    free(t_offsets);
}

#else /* defined(NBLIST_CSR) && defined(__SIMD_NEIGHBOR__) */

/* CSR, vectorized bin scan: same per-thread growable-buffer/prefix-sum
 * structure as the scalar CSR build above, with the same vectorized scan
 * technique as the AOS/SOA SIMD build -- only the storage step differs
 * (CSR_APPEND instead of a neighs() macro write). */
static void buildNeighborLists(Atom* atom, Neighbor* neighbor)
{
    int Nlocal          = atom->Nlocal;
    int max_threads     = 1;
    int actual_nthreads = 1;
#ifdef _OPENMP
    max_threads = omp_get_max_threads();
#endif
    int** t_buf    = (int**)malloc(max_threads * sizeof(int*));
    int*  t_cap    = (int*)malloc(max_threads * sizeof(int));
    int*  t_pos    = (int*)malloc(max_threads * sizeof(int));
    int*  t_owner  = (int*)allocate(ALIGNMENT, Nlocal * sizeof(int));
    int   init_cap = MAX(1, Nlocal / max_threads) * MAX(1, neighbor->maxneighs);
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

#define CSR_APPEND(jval)                                                       \
    do {                                                                       \
        if (t_pos[tid] + n >= t_cap[tid]) {                                    \
            t_cap[tid] *= 2;                                                   \
            t_buf[tid]  = (int*)realloc(t_buf[tid], t_cap[tid] * sizeof(int)); \
        }                                                                      \
        t_buf[tid][t_pos[tid] + n] = (jval);                                   \
        n++;                                                                   \
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
                    MD_SIMD_INT j3 = simd_i32_add(simd_i32_add(jvec, jvec), jvec);
                    MD_SIMD_FLOAT xj =
                        simd_real_gather(j3, &(atom->x[0]), sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT yj =
                        simd_real_gather(j3, &(atom->x[1]), sizeof(MD_FLOAT));
                    MD_SIMD_FLOAT zj =
                        simd_real_gather(j3, &(atom->x[2]), sizeof(MD_FLOAT));
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

                    MD_SIMD_MASK keep = simd_mask_and(
                        bound_mask, simd_mask_not(simd_mask_i32_cond_eq(jvec, ivec)));
                    if (neighbor->half_neigh) {
                        keep = simd_mask_and(
                            keep, simd_mask_not(simd_mask_i32_cond_lt(jvec, ivec)));
                    }

#if LJ_COMB_RULE == LJ_COMB_SINGLE
                    MD_SIMD_FLOAT cutoff_vec = cutoff_vec_single;
#else
                    MD_SIMD_INT type_j_vec = simd_i32_gather(jvec, atom->type, sizeof(int));
                    MD_SIMD_INT tij        = simd_i32_add(tbase_i, type_j_vec);
                    MD_SIMD_FLOAT cutoff_vec =
                        simd_real_gather(tij, atom->cutneighsq, sizeof(MD_FLOAT));
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
    t_offsets[0] = 0;
    for (int t = 0; t < actual_nthreads; t++)
        t_offsets[t + 1] = t_offsets[t] + t_pos[t];
    int total = t_offsets[actual_nthreads];

    if (neighbor->neighbors) free(neighbor->neighbors);
    neighbor->neighbors = (int*)allocate(
        ALIGNMENT, (MAX(1, total) + VECTOR_WIDTH) * sizeof(int));
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
        is_inner_buf = (int*)allocate(ALIGNMENT, MAX(1, neighbor->maxneighs) * sizeof(int));
    }

    free(t_owner);
    free(t_buf);
    free(t_cap);
    free(t_pos);
    free(t_offsets);
}

#endif

void pruneNeighborCPU(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("pruneNeighborCPU begin\n");
    const int nlocal = atom->Nlocal;

    if (!dcut_enabled) {
        // Disabled: mirror outer counts so the force kernel reads the full list.
        for (int i = 0; i < nlocal; i++) {
            neighbor->numneigh_inner[i] = neighbor->numneigh[i];
        }
        DEBUG_MESSAGE("pruneNeighborCPU end\n");
        return;
    }

    const MD_FLOAT cutsq = cutneigh_inner_sq;

#pragma omp parallel
    {
        int* is_inner = is_inner_buf;
#ifdef _OPENMP
        // Each thread needs its own scratch for the in-place partition.
        is_inner = (int*)allocate(ALIGNMENT, neighbor->maxneighs * sizeof(int));
#endif

#pragma omp for schedule(runtime)
        for (int i = 0; i < nlocal; i++) {
            const int numneighs = neighbor->numneigh[i];
            const MD_FLOAT xtmp = atom_x(i);
            const MD_FLOAT ytmp = atom_y(i);
            const MD_FLOAT ztmp = atom_z(i);

            for (int k = 0; k < numneighs; k++) {
                int j         = neighs(neighbor->neighbors, i, k, nlocal, neighbor);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;
                is_inner[k]   = (rsq < cutsq);
            }

            // Partition the list in place: inner neighbors first.
            int lo = 0;
            for (int hi = 0; hi < numneighs; hi++) {
                if (is_inner[hi]) {
                    if (hi != lo) {
                        int t_j  = neighs(neighbor->neighbors, i, lo, nlocal, neighbor);
                        neighs(neighbor->neighbors,
                            i,
                            lo,
                            nlocal,
                            neighbor) = neighs(neighbor->neighbors, i, hi, nlocal, neighbor);
                        neighs(neighbor->neighbors, i, hi, nlocal, neighbor) = t_j;
                        is_inner[hi]                                         = is_inner[lo];
                        is_inner[lo]                                         = 1;
                    }
                    lo++;
                }
            }

            neighbor->numneigh_inner[i] = lo;
        }

#ifdef _OPENMP
        free(is_inner);
#endif
    }

    DEBUG_MESSAGE("pruneNeighborCPU end\n");
}

/* internal subroutines */
MD_FLOAT bindist(int i, int j, int k)
{
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

    if (k > 0) {
        delz = (k - 1) * binsizez;
    } else if (k == 0) {
        delz = 0.0;
    } else {
        delz = (k + 1) * binsizez;
    }

    return (delx * delx + dely * dely + delz * delz);
}

int coord2bin(MD_FLOAT xin, MD_FLOAT yin, MD_FLOAT zin)
{
    int ix, iy, iz;
    MD_FLOAT eps = 1e-9;
    MD_FLOAT xlo = 0.0;
    MD_FLOAT ylo = 0.0;
    MD_FLOAT zlo = 0.0;
    xlo          = fabs(xlo - pad_x * binsizex) + eps;
    ylo          = fabs(ylo - pad_y * binsizey) + eps;
    zlo          = fabs(zlo - pad_z * binsizez) + eps;
    ix           = (int)((xin + xlo) * bininvx);
    iy           = (int)((yin + ylo) * bininvy);
    iz           = (int)((zin + zlo) * bininvz);

    return (iz * mbiny * mbinx + iy * mbinx + ix);

    /*
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

    if (zin >= zprd) {
        iz = (int)((zin - zprd) * bininvz) + nbinz - mbinzlo;
    } else if (zin >= 0.0) {
        iz = (int)(zin * bininvz) - mbinzlo;
    } else {
        iz = (int)(zin * bininvz) - mbinzlo - 1;
    }

    return (iz * mbiny * mbinx + iy * mbinx + ix + 1);
    */
}

void binatoms(Atom* atom)
{
    DEBUG_MESSAGE("binatoms begin\n");

    int nall   = atom->Nlocal + atom->Nghost;
    int resize = 1;

    while (resize > 0) {
        resize = 0;

        for (int i = 0; i < mbins; i++) {
            bincount[i] = 0;
        }

        for (int i = 0; i < nall; i++) {
            int ibin = coord2bin(atom_x(i), atom_y(i), atom_z(i));
            if (shellMethod && !ghostZone(atom, i)) continue;
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
            bins = (int*)allocate(ALIGNMENT, mbins * atoms_per_bin * sizeof(int));
        }
    }

    DEBUG_MESSAGE("binatoms end\n");
}

void sortAtom(Atom* atom)
{
    DEBUG_MESSAGE("sortAtom begin\n");

    binatoms(atom);
    int Nmax    = atom->Nmax;
    int* binpos = bincount;

    for (int i = 1; i < mbins; i++) {
        binpos[i] += binpos[i - 1];
    }

#ifdef ATOM_POSITION_AOS
    MD_FLOAT* new_x  = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT) * 3);
    MD_FLOAT* new_vx = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT) * 3);
#else
    MD_FLOAT* new_x  = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_y  = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_z  = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vx = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vy = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vz = (MD_FLOAT*)allocate(ALIGNMENT, Nmax * sizeof(MD_FLOAT));
#endif
    MD_FLOAT* old_x  = atom->x;
    MD_FLOAT* old_y  = atom->y;
    MD_FLOAT* old_z  = atom->z;
    MD_FLOAT* old_vx = atom->vx;
    MD_FLOAT* old_vy = atom->vy;
    MD_FLOAT* old_vz = atom->vz;

    for (int mybin = 0; mybin < mbins; mybin++) {
        int start = mybin > 0 ? binpos[mybin - 1] : 0;
        int count = binpos[mybin] - start;
        for (int k = 0; k < count; k++) {
            int new_i = start + k;
            int old_i = bins[mybin * atoms_per_bin + k];
#ifdef ATOM_POSITION_AOS
            new_x[new_i * 3 + 0]  = old_x[old_i * 3 + 0];
            new_x[new_i * 3 + 1]  = old_x[old_i * 3 + 1];
            new_x[new_i * 3 + 2]  = old_x[old_i * 3 + 2];
            new_vx[new_i * 3 + 0] = old_vx[old_i * 3 + 0];
            new_vx[new_i * 3 + 1] = old_vx[old_i * 3 + 1];
            new_vx[new_i * 3 + 2] = old_vx[old_i * 3 + 2];
#else
            new_x[new_i]  = old_x[old_i];
            new_y[new_i]  = old_y[old_i];
            new_z[new_i]  = old_z[old_i];
            new_vx[new_i] = old_vx[old_i];
            new_vy[new_i] = old_vy[old_i];
            new_vz[new_i] = old_vz[old_i];
#endif
        }
    }

    free(atom->x);
    free(atom->vx);
    atom->x  = new_x;
    atom->vx = new_vx;
#ifndef ATOM_POSITION_AOS
    free(atom->y);
    free(atom->z);
    free(atom->vy);
    free(atom->vz);
    atom->y  = new_y;
    atom->z  = new_z;
    atom->vy = new_vy;
    atom->vz = new_vz;
#endif

    DEBUG_MESSAGE("sortAtom end\n");
}

/* internal subroutines
Added with MPI*/

static int ghostZone(Atom* atom, int i)
{
    if (i < atom->Nlocal) return 1;
    else if (method == halfShell)
        return halfZone(atom, i);
    else if (method == eightShell)
        return eightZone(atom, i);
    else
        return 0;
}

static int eightZone(Atom* atom, int i)
{
    // Mapping: 0->0, 1->1, 2->2, 3->6, 4->3, 5->5, 6->4, 7->7
    int zoneMapping[] = { 0, 1, 2, 6, 3, 5, 4, 7 };
    MD_FLOAT* hi      = atom->mybox.hi;
    int zone          = 0;

    if (BigOrEqual(atom_x(i), hi[0])) {
        zone += 1;
    }

    if (BigOrEqual(atom_y(i), hi[1])) {
        zone += 2;
    }

    if (BigOrEqual(atom_z(i), hi[2])) {
        zone += 4;
    }

    return zoneMapping[zone];
}

static int halfZone(Atom* atom, int i)
{
    MD_FLOAT* hi = atom->mybox.hi;
    MD_FLOAT* lo = atom->mybox.lo;

    if (atom_x(i) < lo[0] && atom_y(i) < hi[1] && atom_z(i) < hi[2]) {
        return 0;
    } else if (atom_y(i) < lo[1] && atom_z(i) < hi[2]) {
        return 0;
    } else if (atom_z(i) < lo[2]) {
        return 0;
    } else {
        return 1;
    }
}

static void neighborGhost(Atom* atom, Neighbor* neighbor)
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

static inline int skipNeigh(Atom* atom, int i, int j)
{
    if (i > j && j < atom->Nlocal) {
        return 1;
    } else if (atom_z(i) > atom_z(j) && j >= atom->Nlocal) {
        return 1;
    } else if (Equal(atom_z(i), atom_z(j)) && atom_y(i) > atom_y(j) &&
               j >= atom->Nlocal) {
        return 1;
    } else if (Equal(atom_z(i), atom_z(j)) && Equal(atom_y(i), atom_y(j)) &&
               atom_x(i) > atom_x(j) && j >= atom->Nlocal) {
        return 1;
    } else {
        return 0;
    }
}
