/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//---

#include <device.h>

extern "C" {
#include <atom.h>
#include <neighbor.h>
#include <parameter.h>
#include <util.h>
}

extern MD_FLOAT xprd, yprd, zprd;
extern MD_FLOAT bininvx, bininvy, bininvz;
extern int mbinxlo, mbinylo, mbinzlo;
extern int nbinx, nbiny, nbinz;
extern int mbinx, mbiny, mbinz; // n bins in x, y, z
extern int mbins;               // total number of bins
extern int atoms_per_bin;       // max atoms per bin
extern MD_FLOAT cutneighsq;     // neighbor cutoff squared
extern int nmax;
extern int nstencil; // # of bins in stencil
extern int* stencil; // stencil list of bin offsets
static int* c_stencil       = NULL;
static int* c_resize_needed = NULL;
static int* c_new_maxneighs = NULL;
static Binning c_binning {
    .bincount = NULL, .bins = NULL, .mbins = 0, .atoms_per_bin = 0
};

// Multi GPU
extern int pad_x, pad_y, pad_z;
extern MD_FLOAT binsizex, binsizey, binsizez;

__device__ int coord2bin_device(
    MD_FLOAT xin, MD_FLOAT yin, MD_FLOAT zin, Neighbor_params np)
{
    /*
    int ix, iy, iz;

    if (xin >= np.xprd) {
        ix = (int)((xin - np.xprd) * np.bininvx) + np.nbinx - np.mbinxlo;
    } else if (xin >= 0.0) {
        ix = (int)(xin * np.bininvx) - np.mbinxlo;
    } else {
        ix = (int)(xin * np.bininvx) - np.mbinxlo - 1;
    }

    if (yin >= np.yprd) {
        iy = (int)((yin - np.yprd) * np.bininvy) + np.nbiny - np.mbinylo;
    } else if (yin >= 0.0) {
        iy = (int)(yin * np.bininvy) - np.mbinylo;
    } else {
        iy = (int)(yin * np.bininvy) - np.mbinylo - 1;
    }

    if (zin >= np.zprd) {
        iz = (int)((zin - np.zprd) * np.bininvz) + np.nbinz - np.mbinzlo;
    } else if (zin >= 0.0) {
        iz = (int)(zin * np.bininvz) - np.mbinzlo;
    } else {
        iz = (int)(zin * np.bininvz) - np.mbinzlo - 1;
    }

    return (iz * np.mbiny * np.mbinx + iy * np.mbinx + ix + 1);
    */

    int ix, iy, iz;
    MD_FLOAT eps = 1e-9;
    MD_FLOAT xlo = 0.0;
    MD_FLOAT ylo = 0.0;
    MD_FLOAT zlo = 0.0;
    xlo          = fabs(xlo - np.pad_x * np.binsizex) + eps;
    ylo          = fabs(ylo - np.pad_y * np.binsizey) + eps;
    zlo          = fabs(zlo - np.pad_z * np.binsizez) + eps;
    ix           = (int)((xin + xlo) * np.bininvx);
    iy           = (int)((yin + ylo) * np.bininvy);
    iz           = (int)((zin + zlo) * np.bininvz);
    return (iz * np.mbiny * np.mbinx + iy * np.mbinx + ix);
}

/* sorts the contents of a bin to make it comparable to the CPU version */
/* uses bubble sort since atoms per bin should be relatively small and can be done in situ
 */
__global__ void sort_bin_contents_kernel(
    int* bincount, int* bins, int mbins, int atoms_per_bin)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= mbins) {
        return;
    }

    int atoms_in_bin = bincount[i];
    int* bin_ptr     = &bins[i * atoms_per_bin];
    int sorted;
    do {
        sorted = 1;
        int tmp;
        for (int index = 0; index < atoms_in_bin - 1; index++) {
            if (bin_ptr[index] > bin_ptr[index + 1]) {
                tmp                = bin_ptr[index];
                bin_ptr[index]     = bin_ptr[index + 1];
                bin_ptr[index + 1] = tmp;
                sorted             = 0;
            }
        }
    } while (!sorted);
}

__global__ void binatoms_kernel(DeviceAtom a,
    int nall,
    int* bincount,
    int* bins,
    int atoms_per_bin,
    Neighbor_params np,
    int* resize_needed)
{
    DeviceAtom* atom = &a;
    const int i      = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nall) {
        return;
    }

    MD_FLOAT x = atom_x(i);
    MD_FLOAT y = atom_y(i);
    MD_FLOAT z = atom_z(i);
    int ibin   = coord2bin_device(x, y, z, np);
    int ac     = atomicAdd(&bincount[ibin], 1);

    if (ac < atoms_per_bin) {
        bins[ibin * atoms_per_bin + ac] = i;
    } else {
        atomicMax(resize_needed, ac);
    }
}

__global__ void compute_neighborhood(DeviceAtom a,
    DeviceNeighbor neigh,
    Neighbor_params np,
    int nlocal,
    int maxneighs,
    int halfneigh,
    int nstencil,
    int* stencil,
    int* bins,
    int atoms_per_bin,
    int* bincount,
    int* new_maxneighs,
    MD_FLOAT cutneighsq,
    int ntypes)
{

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nlocal) {
        return;
    }

    DeviceAtom* atom         = &a;
    DeviceNeighbor* neighbor = &neigh;

    int n         = 0;
    MD_FLOAT xtmp = atom_x(i);
    MD_FLOAT ytmp = atom_y(i);
    MD_FLOAT ztmp = atom_z(i);
    int ibin      = coord2bin_device(xtmp, ytmp, ztmp, np);
#if LJ_COMB_RULE != LJ_COMB_SINGLE
    int type_i = atom->type[i];
#endif
    for (int k = 0; k < nstencil; k++) {
        int jbin     = ibin + stencil[k];
        int* loc_bin = &bins[jbin * atoms_per_bin];

        for (int m = 0; m < bincount[jbin]; m++) {
            int j = loc_bin[m];
            if (j == i || (halfneigh && (j < i))) {
                continue;
            }

            MD_FLOAT delx = xtmp - atom_x(j);
            MD_FLOAT dely = ytmp - atom_y(j);
            MD_FLOAT delz = ztmp - atom_z(j);
            MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            int type_j            = atom->type[j];
            const MD_FLOAT cutoff = atom->cutneighsq[type_i * ntypes + type_j];
#else
            const MD_FLOAT cutoff = cutneighsq;
#endif

            if (rsq <= cutoff) {
                neighs(neighbor->neighbors, i, n, nlocal, maxneighs) = j;
                n++;
            }
        }
    }

    neighbor->numneigh[i] = n;
    if (n > maxneighs) {
        atomicMax(new_maxneighs, n);
    }
}

void binatoms_cuda(Atom* atom,
    Binning* c_binning,
    int* c_resize_needed,
    Neighbor_params* np,
    const int threads_per_block)
{
    DEBUG_MESSAGE("binatoms_cuda begin\n");

    int nall             = atom->Nlocal + atom->Nghost;
    int resize           = 1;
    const int num_blocks = ceil((float)nall / (float)threads_per_block);

    while (resize > 0) {
        resize = 0;
        memsetGPU(c_binning->bincount, 0, c_binning->mbins * sizeof(int));
        memsetGPU(c_resize_needed, 0, sizeof(int));
        binatoms_kernel<<<num_blocks, threads_per_block>>>(atom->d_atom,
            atom->Nlocal + atom->Nghost,
            c_binning->bincount,
            c_binning->bins,
            c_binning->atoms_per_bin,
            *np,
            c_resize_needed);
        cuda_assert("binatoms", cudaPeekAtLastError());
        cuda_assert("binatoms", cudaDeviceSynchronize());
        memcpyFromGPU(&resize, c_resize_needed, sizeof(int));
        if (resize) {
            c_binning->atoms_per_bin *= 2;
            c_binning->bins = (int*)reallocateGPU(c_binning->bins,
                c_binning->mbins * c_binning->atoms_per_bin * sizeof(int));
        }
    }

    atoms_per_bin        = c_binning->atoms_per_bin;
    const int sortBlocks = ceil((float)mbins / (float)threads_per_block);
    sort_bin_contents_kernel<<<sortBlocks, threads_per_block>>>(c_binning->bincount,
        c_binning->bins,
        c_binning->mbins,
        c_binning->atoms_per_bin);
    cuda_assert("sort_bin", cudaPeekAtLastError());
    cuda_assert("sort_bin", cudaDeviceSynchronize());

    DEBUG_MESSAGE("binatoms_cuda end\n");
}

void buildNeighborCUDA(Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("buildNeighborCUDA begin\n");
    DeviceNeighbor* d_neighbor      = &(neighbor->d_neighbor);
    const int num_threads_per_block = get_cuda_num_threads();
    cudaProfilerStart();

    int nall = atom->Nlocal + atom->Nghost;
    if (nall > nmax) {
        nmax                  = nall;
        d_neighbor->neighbors = (int*)reallocateGPU(d_neighbor->neighbors,
            nmax * neighbor->maxneighs * sizeof(int*));
        d_neighbor->numneigh  = (int*)reallocateGPU(d_neighbor->numneigh,
            nmax * sizeof(int));
    }

    // TODO move all of this initialization into its own method
    if (c_stencil == NULL) {
        c_stencil = (int*)allocateGPU(nstencil * sizeof(int));
        memcpyToGPU(c_stencil, stencil, nstencil * sizeof(int));
    }

    if (c_binning.mbins == 0) {
        c_binning.mbins         = mbins;
        c_binning.atoms_per_bin = atoms_per_bin;
        c_binning.bincount      = (int*)allocateGPU(c_binning.mbins * sizeof(int));
        c_binning.bins          = (int*)allocateGPU(
            c_binning.mbins * c_binning.atoms_per_bin * sizeof(int));
    }

    Neighbor_params np { .xprd = xprd,
        .yprd                  = yprd,
        .zprd                  = zprd,
        .bininvx               = bininvx,
        .bininvy               = bininvy,
        .bininvz               = bininvz,
        .mbinxlo               = mbinxlo,
        .mbinylo               = mbinylo,
        .mbinzlo               = mbinzlo,
        .nbinx                 = nbinx,
        .nbiny                 = nbiny,
        .nbinz                 = nbinz,
        .mbinx                 = mbinx,
        .mbiny                 = mbiny,
        .mbinz                 = mbinz,
        // MultiGPU
        .pad_x    = pad_x,
        .pad_y    = pad_y,
        .pad_z    = pad_z,
        .binsizex = binsizex,
        .binsizey = binsizey,
        .binsizez = binsizez };

    if (c_resize_needed == NULL) {
        c_resize_needed = (int*)allocateGPU(sizeof(int));
    }

    /* bin local & ghost atoms */
    binatoms_cuda(atom, &c_binning, c_resize_needed, &np, num_threads_per_block);
    if (c_new_maxneighs == NULL) {
        c_new_maxneighs = (int*)allocateGPU(sizeof(int));
    }

    int resize = 1;

    /* loop over each atom, storing neighbors */
    while (resize) {
        resize = 0;
        memsetGPU(c_new_maxneighs, 0, sizeof(int));
        const int num_blocks = ceil((float)atom->Nlocal / (float)num_threads_per_block);
        compute_neighborhood<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
            *d_neighbor,
            np,
            atom->Nlocal,
            neighbor->maxneighs,
            neighbor->half_neigh,
            nstencil,
            c_stencil,
            c_binning.bins,
            c_binning.atoms_per_bin,
            c_binning.bincount,
            c_new_maxneighs,
            cutneighsq,
            atom->ntypes);

        cuda_assert("compute_neighborhood", cudaPeekAtLastError());
        cuda_assert("compute_neighborhood", cudaDeviceSynchronize());

        int new_maxneighs;
        memcpyFromGPU(&new_maxneighs, c_new_maxneighs, sizeof(int));
        if (new_maxneighs > neighbor->maxneighs) {
            resize = 1;
        }

        if (resize) {
            printf("RESIZE %d\n", neighbor->maxneighs);
            neighbor->maxneighs = new_maxneighs * 1.2;
            printf("NEW SIZE %d\n", neighbor->maxneighs);
            neighbor->neighbors = (int*)reallocateGPU(neighbor->neighbors,
                atom->Nmax * neighbor->maxneighs * sizeof(int));
        }
    }

    cudaProfilerStop();
    DEBUG_MESSAGE("buildNeighborCUDA end\n");
}
