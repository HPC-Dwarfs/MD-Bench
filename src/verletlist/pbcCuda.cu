/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>

// C headers must come before device.h (which includes atom.h without C
// linkage) so that growAtom resolves against the C-compiled atom.o
extern "C" {
#include <allocate.h>
#include <atom.h>
#include <pbc.h>
#include <util.h>
}

#include <cub/cub.cuh>
#include <device.h>

extern "C" {

static int c_NmaxGhost = 0;
static int *c_PBCx = NULL, *c_PBCy = NULL, *c_PBCz = NULL;
// Exclusive per-atom ghost offsets (Nlocal+1 ints; slot [Nlocal] = total count)
static int* c_ghostOffsets  = NULL;
static int c_offsetsCap     = 0;
static void* c_scanTemp     = NULL;
static size_t c_scanTempCap = 0;

__global__ void computeAtomsPbcUpdate(
    DeviceAtom a, int nlocal, MD_FLOAT xprd, MD_FLOAT yprd, MD_FLOAT zprd)
{
    const int i      = blockIdx.x * blockDim.x + threadIdx.x;
    DeviceAtom* atom = &a;
    if (i >= nlocal) {
        return;
    }

    if (atom_x(i) < 0.0) {
        atom_x(i) += xprd;
    } else if (atom_x(i) >= xprd) {
        atom_x(i) -= xprd;
    }

    if (atom_y(i) < 0.0) {
        atom_y(i) += yprd;
    } else if (atom_y(i) >= yprd) {
        atom_y(i) -= yprd;
    }

    if (atom_z(i) < 0.0) {
        atom_z(i) += zprd;
    } else if (atom_z(i) >= zprd) {
        atom_z(i) -= zprd;
    }
}

__global__ void computePbcUpdate(DeviceAtom a,
    int nlocal,
    int nghost,
    int* PBCx,
    int* PBCy,
    int* PBCz,
    MD_FLOAT xprd,
    MD_FLOAT yprd,
    MD_FLOAT zprd)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nghost) {
        return;
    }

    DeviceAtom* atom   = &a;
    int* border_map    = atom->border_map;
    atom_x(nlocal + i) = atom_x(border_map[i]) + PBCx[i] * xprd;
    atom_y(nlocal + i) = atom_y(border_map[i]) + PBCy[i] * yprd;
    atom_z(nlocal + i) = atom_z(border_map[i]) + PBCz[i] * zprd;
    // Copy per-atom type and LJ parameters for the combination rules
    atom->type[nlocal + i]         = atom->type[border_map[i]];
    atom->sqrt_epsilon[nlocal + i] = atom->sqrt_epsilon[border_map[i]];
    atom->sigma3[nlocal + i]       = atom->sigma3[border_map[i]];
}

/* Enumerate the ghost images of one atom in the exact order of the host
 * setupPbc (pbc.c): 6 planes, 8 corners, 12 edges. Count and fill kernels
 * both call this helper so the two passes can never diverge. An atom
 * produces at most 7 ghosts (3 planes + 3 edges + 1 corner). */
__device__ static inline int computeGhostShifts(MD_FLOAT x,
    MD_FLOAT y,
    MD_FLOAT z,
    MD_FLOAT xprd,
    MD_FLOAT yprd,
    MD_FLOAT zprd,
    MD_FLOAT cutneigh,
    int pbc_x,
    int pbc_y,
    int pbc_z,
    int* sx,
    int* sy,
    int* sz)
{
    int n = 0;
#define ADDSHIFT(dx, dy, dz)                                                             \
    {                                                                                    \
        sx[n] = dx;                                                                      \
        sy[n] = dy;                                                                      \
        sz[n] = dz;                                                                      \
        n++;                                                                             \
    }

    /* 6 planes */
    if (pbc_x != 0) {
        if (x < cutneigh) {
            ADDSHIFT(+1, 0, 0);
        }
        if (x >= (xprd - cutneigh)) {
            ADDSHIFT(-1, 0, 0);
        }
    }

    if (pbc_y != 0) {
        if (y < cutneigh) {
            ADDSHIFT(0, +1, 0);
        }
        if (y >= (yprd - cutneigh)) {
            ADDSHIFT(0, -1, 0);
        }
    }

    if (pbc_z != 0) {
        if (z < cutneigh) {
            ADDSHIFT(0, 0, +1);
        }
        if (z >= (zprd - cutneigh)) {
            ADDSHIFT(0, 0, -1);
        }
    }

    /* 8 corners */
    if (pbc_x != 0 && pbc_y != 0 && pbc_z != 0) {
        if (x < cutneigh && y < cutneigh && z < cutneigh) {
            ADDSHIFT(+1, +1, +1);
        }
        if (x < cutneigh && y >= (yprd - cutneigh) && z < cutneigh) {
            ADDSHIFT(+1, -1, +1);
        }
        if (x < cutneigh && y < cutneigh && z >= (zprd - cutneigh)) {
            ADDSHIFT(+1, +1, -1);
        }
        if (x < cutneigh && y >= (yprd - cutneigh) && z >= (zprd - cutneigh)) {
            ADDSHIFT(+1, -1, -1);
        }
        if (x >= (xprd - cutneigh) && y < cutneigh && z < cutneigh) {
            ADDSHIFT(-1, +1, +1);
        }
        if (x >= (xprd - cutneigh) && y >= (yprd - cutneigh) && z < cutneigh) {
            ADDSHIFT(-1, -1, +1);
        }
        if (x >= (xprd - cutneigh) && y < cutneigh && z >= (zprd - cutneigh)) {
            ADDSHIFT(-1, +1, -1);
        }
        if (x >= (xprd - cutneigh) && y >= (yprd - cutneigh) &&
            z >= (zprd - cutneigh)) {
            ADDSHIFT(-1, -1, -1);
        }
    }

    /* 12 edges */
    if (pbc_x != 0 && pbc_z != 0) {
        if (x < cutneigh && z < cutneigh) {
            ADDSHIFT(+1, 0, +1);
        }
        if (x < cutneigh && z >= (zprd - cutneigh)) {
            ADDSHIFT(+1, 0, -1);
        }
        if (x >= (xprd - cutneigh) && z < cutneigh) {
            ADDSHIFT(-1, 0, +1);
        }
        if (x >= (xprd - cutneigh) && z >= (zprd - cutneigh)) {
            ADDSHIFT(-1, 0, -1);
        }
    }

    if (pbc_y != 0 && pbc_z != 0) {
        if (y < cutneigh && z < cutneigh) {
            ADDSHIFT(0, +1, +1);
        }
        if (y < cutneigh && z >= (zprd - cutneigh)) {
            ADDSHIFT(0, +1, -1);
        }
        if (y >= (yprd - cutneigh) && z < cutneigh) {
            ADDSHIFT(0, -1, +1);
        }
        if (y >= (yprd - cutneigh) && z >= (zprd - cutneigh)) {
            ADDSHIFT(0, -1, -1);
        }
    }

    if (pbc_x != 0 && pbc_y != 0) {
        if (y < cutneigh && x < cutneigh) {
            ADDSHIFT(+1, +1, 0);
        }
        if (y < cutneigh && x >= (xprd - cutneigh)) {
            ADDSHIFT(-1, +1, 0);
        }
        if (y >= (yprd - cutneigh) && x < cutneigh) {
            ADDSHIFT(+1, -1, 0);
        }
        if (y >= (yprd - cutneigh) && x >= (xprd - cutneigh)) {
            ADDSHIFT(-1, -1, 0);
        }
    }
#undef ADDSHIFT
    return n;
}

__global__ void computeGhostCount(DeviceAtom a,
    int nlocal,
    MD_FLOAT xprd,
    MD_FLOAT yprd,
    MD_FLOAT zprd,
    MD_FLOAT cutneigh,
    int pbc_x,
    int pbc_y,
    int pbc_z,
    int* counts)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i > nlocal) {
        return;
    }

    if (i == nlocal) {
        // sentinel so the exclusive scan yields the total in counts[nlocal]
        counts[nlocal] = 0;
        return;
    }

    DeviceAtom* atom = &a;
    int sx[7], sy[7], sz[7];
    counts[i] = computeGhostShifts(atom_x(i),
        atom_y(i),
        atom_z(i),
        xprd,
        yprd,
        zprd,
        cutneigh,
        pbc_x,
        pbc_y,
        pbc_z,
        sx,
        sy,
        sz);
}

__global__ void computeGhostFill(DeviceAtom a,
    int nlocal,
    MD_FLOAT xprd,
    MD_FLOAT yprd,
    MD_FLOAT zprd,
    MD_FLOAT cutneigh,
    int pbc_x,
    int pbc_y,
    int pbc_z,
    const int* offsets,
    int* PBCx,
    int* PBCy,
    int* PBCz)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nlocal) {
        return;
    }

    DeviceAtom* atom = &a;
    int sx[7], sy[7], sz[7];
    const int n = computeGhostShifts(atom_x(i),
        atom_y(i),
        atom_z(i),
        xprd,
        yprd,
        zprd,
        cutneigh,
        pbc_x,
        pbc_y,
        pbc_z,
        sx,
        sy,
        sz);

    const int off   = offsets[i];
    int* border_map = atom->border_map;
    for (int k = 0; k < n; k++) {
        border_map[off + k] = i;
        PBCx[off + k]       = sx[k];
        PBCy[off + k]       = sy[k];
        PBCz[off + k]       = sz[k];
    }
}

/* fully device-resident setupPbc: builds border_map and the PBC shift maps
 * on the GPU (count -> exclusive scan -> fill), preserving the exact ghost
 * ordering of the host setupPbc. Only the total ghost count is read back
 * (one int) — no positions leave the device and no maps are uploaded. */
void setupPbcCUDA(Atom* atom, Parameter* param)
{
    const int num_threads_per_block = get_cuda_num_threads();
    const int nlocal                = atom->Nlocal;

    if (c_offsetsCap < nlocal + 1) {
        c_offsetsCap   = nlocal + 1;
        c_ghostOffsets = (int*)reallocateGPU(c_ghostOffsets,
            c_offsetsCap * sizeof(int));
    }

    int num_blocks = ceil((float)(nlocal + 1) / (float)num_threads_per_block);
    computeGhostCount<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
        nlocal,
        param->xprd,
        param->yprd,
        param->zprd,
        param->cutneigh,
        param->pbc_x,
        param->pbc_y,
        param->pbc_z,
        c_ghostOffsets);
    cuda_assert("computeGhostCount", cudaPeekAtLastError());

    // in-place exclusive scan turns per-atom counts into ghost offsets
    size_t tempBytes = 0;
    cuda_assert("setupPbc.scanSize",
        cub::DeviceScan::ExclusiveSum(NULL,
            tempBytes,
            c_ghostOffsets,
            c_ghostOffsets,
            nlocal + 1));
    if (tempBytes > c_scanTempCap) {
        c_scanTempCap = tempBytes;
        c_scanTemp    = reallocateGPU(c_scanTemp, tempBytes);
    }
    cuda_assert("setupPbc.scan",
        cub::DeviceScan::ExclusiveSum(c_scanTemp,
            tempBytes,
            c_ghostOffsets,
            c_ghostOffsets,
            nlocal + 1));

    // the single remaining transfer per reneighboring (4 bytes, also syncs)
    int nghost;
    memcpyFromGPU(&nghost, &c_ghostOffsets[nlocal], sizeof(int));
    atom->Nghost = nghost;

    // capacity checks stay host-driven: growAtom preserves device contents
    // (reallocateGPUKeep); the maps are fully rewritten below, so plain
    // reallocateGPU suffices for them
    while (atom->Nmax < nlocal + nghost) {
        growAtom(atom);
    }

    if (c_NmaxGhost < nghost) {
        c_NmaxGhost = nghost + nghost / 8; // slack so Nghost jitter doesn't realloc
        c_PBCx      = (int*)reallocateGPU(c_PBCx, c_NmaxGhost * sizeof(int));
        c_PBCy      = (int*)reallocateGPU(c_PBCy, c_NmaxGhost * sizeof(int));
        c_PBCz      = (int*)reallocateGPU(c_PBCz, c_NmaxGhost * sizeof(int));
        atom->d_atom.border_map = (int*)reallocateGPU(atom->d_atom.border_map,
            c_NmaxGhost * sizeof(int));
    }

    num_blocks = ceil((float)nlocal / (float)num_threads_per_block);
    computeGhostFill<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
        nlocal,
        param->xprd,
        param->yprd,
        param->zprd,
        param->cutneigh,
        param->pbc_x,
        param->pbc_y,
        param->pbc_z,
        c_ghostOffsets,
        c_PBCx,
        c_PBCy,
        c_PBCz);
    cuda_assert("computeGhostFill", cudaPeekAtLastError());
    cuda_assert("computeGhostFill", cudaDeviceSynchronize());
}

/* update coordinates of ghost atoms */
/* uses mapping created in setupPbc */
void updatePbcCUDA(Atom* atom, Parameter* param, bool reneigh)
{
    const int num_threads_per_block = get_cuda_num_threads();

    MD_FLOAT xprd = param->xprd;
    MD_FLOAT yprd = param->yprd;
    MD_FLOAT zprd = param->zprd;

    const int num_blocks = ceil((float)atom->Nghost / (float)num_threads_per_block);
    computePbcUpdate<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
        atom->Nlocal,
        atom->Nghost,
        c_PBCx,
        c_PBCy,
        c_PBCz,
        xprd,
        yprd,
        zprd);
    cuda_assert("updatePbc", cudaPeekAtLastError());
    cuda_assert("updatePbc", cudaDeviceSynchronize());
}

void updateAtomsPbcCUDA(Atom* atom, Parameter* param, bool reneigh)
{
    const int num_threads_per_block = get_cuda_num_threads();
    MD_FLOAT xprd                   = param->xprd;
    MD_FLOAT yprd                   = param->yprd;
    MD_FLOAT zprd                   = param->zprd;

    const int num_blocks = ceil((float)atom->Nlocal / (float)num_threads_per_block);
    computeAtomsPbcUpdate<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
        atom->Nlocal,
        xprd,
        yprd,
        zprd);
    cuda_assert("computeAtomsPbcUpdate", cudaPeekAtLastError());
    cuda_assert("computeAtomsPbcUpdate", cudaDeviceSynchronize());
}
}
