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
#include <device.h>
#include <kernel_variant.h>
#include <neighbor.h>
#include <neighbor_internal.h>
#include <parameter.h>
#include <util.h>

#ifdef _OPENMP
#include <omp.h>
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
static MD_FLOAT bindist(int, int, int);

// MPI version
int pad_x, pad_y, pad_z; // version MPI
int me;                  // rank
int method;              // method
int half_stencil;        // If half stencil exist
int shellMethod;         // If shell method exist

/* exported subroutines */
void initNeighbor(Neighbor* neighbor, Parameter* param)
{
    // The list build variant, setup() sets it again after resolving --lj-comb auto
    neighbor->lj_comb_rule = param->lj_comb_rule;
    neighbor->simd         = param->neighbor_simd;
#ifndef HAVE_SIMD_NEIGHBOR
    if (neighbor->simd != SIMD_MODE_AUTOVEC) {
        fprintf(stderr,
            "Error: --neighbor-simd %s needs SIMD=AVX2/AVX512/NEON/SVE!\n",
            simdmode2str(neighbor->simd));
        exit(EXIT_FAILURE);
    }
#endif
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
    neighbor->tpa_log2       = 0;
    neighbor->numneigh       = NULL;
    neighbor->numneigh_inner = NULL;
    neighbor->neighbors      = NULL;
    neighbor->neigh_start    = NULL;
    is_inner_buf = (int*)allocate(ALIGNMENT, neighbor->maxneighs * sizeof(int));
    method       = param->method;
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

    if (param->input_file != NULL || param->gmxbenchmark) {
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

    if (param->input_file != NULL || param->gmxbenchmark) {
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
    // The list build is compiled once per LJ combination rule (neighbor_build.c)
#ifdef HAVE_SIMD_NEIGHBOR
    if (neighbor->simd == SIMD_MODE_INTRINSICS) {
        CALL_RULE_VARIANT(neighbor->lj_comb_rule, buildNeighborListsSimd, atom, neighbor);
    } else {
        CALL_RULE_VARIANT(neighbor->lj_comb_rule, buildNeighborLists, atom, neighbor);
    }
#else
    CALL_RULE_VARIANT(neighbor->lj_comb_rule, buildNeighborLists, atom, neighbor);
#endif

    if (method == eightShell) {
        CALL_RULE_VARIANT(neighbor->lj_comb_rule, neighborGhost, atom, neighbor);
    }

    pruneNeighborCPU(NULL, atom, neighbor);
    DEBUG_MESSAGE("buildNeighborCPU end\n");
}

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
                        int t_j = neighs(neighbor->neighbors, i, lo, nlocal, neighbor);
                        neighs(neighbor->neighbors, i, lo, nlocal, neighbor) =
                            neighs(neighbor->neighbors, i, hi, nlocal, neighbor);
                        neighs(neighbor->neighbors, i, hi, nlocal, neighbor) = t_j;
                        is_inner[hi] = is_inner[lo];
                        is_inner[lo] = 1;
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
            if (ibin < 0 || ibin >= mbins) {
                fprintf(stderr,
                    "Error: Atom %d assigned to invalid bin %d, position: %.4f, %.4f, "
                    "%.4f\n",
                    i,
                    ibin,
                    atom_x(i),
                    atom_y(i),
                    atom_z(i));
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
    MD_FLOAT* new_x  = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT) * 3);
    MD_FLOAT* new_vx = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT) * 3);
#else
    MD_FLOAT* new_x  = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_y  = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_z  = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vx = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vy = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vz = (MD_FLOAT*)allocateHostPinned(Nmax * sizeof(MD_FLOAT));
#endif
    MD_FLOAT* old_x  = atom->x;
    MD_FLOAT* old_y  = atom->y;
    MD_FLOAT* old_z  = atom->z;
    MD_FLOAT* old_vx = atom->vx;
    MD_FLOAT* old_vy = atom->vy;
    MD_FLOAT* old_vz = atom->vz;
    // The type and id move with the atom, the rest of its data derives from them
    int* new_type = (int*)allocate(ALIGNMENT, Nmax * sizeof(int));
    int* new_id   = (int*)allocate(ALIGNMENT, Nmax * sizeof(int));

    for (int mybin = 0; mybin < mbins; mybin++) {
        int start = mybin > 0 ? binpos[mybin - 1] : 0;
        int count = binpos[mybin] - start;
        for (int k = 0; k < count; k++) {
            int new_i       = start + k;
            int old_i       = bins[mybin * atoms_per_bin + k];
            new_type[new_i] = atom->type[old_i];
            new_id[new_i]   = atom->id[old_i];
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

    freeHostPinned(atom->x);
    freeHostPinned(atom->vx);
    atom->x  = new_x;
    atom->vx = new_vx;
#ifndef ATOM_POSITION_AOS
    freeHostPinned(atom->y);
    freeHostPinned(atom->z);
    freeHostPinned(atom->vy);
    freeHostPinned(atom->vz);
    atom->y  = new_y;
    atom->z  = new_z;
    atom->vy = new_vy;
    atom->vz = new_vz;
#endif

    for (int i = 0; i < atom->Nlocal; i++) {
        setAtomTypeId(atom, i, new_type[i], new_id[i]);
    }
    free(new_type);
    free(new_id);

    DEBUG_MESSAGE("sortAtom end\n");
}

/* internal subroutines
Added with MPI*/

int ghostZone(Atom* atom, int i)
{
    if (i < atom->Nlocal) return 1;
    else if (method == halfShell)
        return halfZone(atom, i);
    else if (method == eightShell)
        return eightZone(atom, i);
    else
        return 0;
}

int eightZone(Atom* atom, int i)
{
    // Mapping: 0->0, 1->1, 2->2, 3->6, 4->3, 5->5, 6->4, 7->7
    int zoneMapping[] = { 0, 1, 2, 6, 3, 5, 4, 7 };
    MD_FLOAT* hi      = atom->mybox.hi;
    MD_FLOAT* lo      = atom->mybox.lo;
    int zone          = 0;

    // The eighth shell only imports ghosts from the positive directions. Without MPI the
    // periodic ghosts surround the box on all sides, the others are not used
    if (atom_x(i) < lo[0] || atom_y(i) < lo[1] || atom_z(i) < lo[2]) {
        return 0;
    }

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

int halfZone(Atom* atom, int i)
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
