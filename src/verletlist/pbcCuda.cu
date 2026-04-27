/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>

#include <device.h>

extern "C" {
#include <allocate.h>
#include <atom.h>
#include <pbc.h>
#include <util.h>

extern int NmaxGhost;
extern int *PBCx, *PBCy, *PBCz;
static int c_NmaxGhost = 0;
static int *c_PBCx = NULL, *c_PBCy = NULL, *c_PBCz = NULL;

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
    // Copy per-atom LJ parameters for geometric combination rule
    atom->sqrt_epsilon[nlocal + i] = atom->sqrt_epsilon[border_map[i]];
    atom->sigma3[nlocal + i]       = atom->sigma3[border_map[i]];
}

/* update coordinates of ghost atoms */
/* uses mapping created in setupPbc */
void updatePbcCUDA(Atom* atom, Parameter* param, bool reneigh)
{
    const int num_threads_per_block = get_cuda_num_threads();

    if (reneigh) {
        memcpyToGPU(atom->d_atom.x, atom->x, sizeof(MD_FLOAT) * atom->Nmax * 3);
        memcpyToGPU(atom->d_atom.type, atom->type, sizeof(int) * atom->Nmax);
        memcpyToGPU(atom->d_atom.sqrt_epsilon,
            atom->sqrt_epsilon,
            sizeof(MD_FLOAT) * atom->Nmax);
        memcpyToGPU(atom->d_atom.sigma3, atom->sigma3, sizeof(MD_FLOAT) * atom->Nmax);

        if (c_NmaxGhost < NmaxGhost) {
            c_NmaxGhost = NmaxGhost;
            c_PBCx      = (int*)reallocateGPU(c_PBCx, NmaxGhost * sizeof(int));
            c_PBCy      = (int*)reallocateGPU(c_PBCy, NmaxGhost * sizeof(int));
            c_PBCz      = (int*)reallocateGPU(c_PBCz, NmaxGhost * sizeof(int));
            atom->d_atom.border_map = (int*)reallocateGPU(atom->d_atom.border_map,
                NmaxGhost * sizeof(int));
        }

        memcpyToGPU(c_PBCx, PBCx, NmaxGhost * sizeof(int));
        memcpyToGPU(c_PBCy, PBCy, NmaxGhost * sizeof(int));
        memcpyToGPU(c_PBCz, PBCz, NmaxGhost * sizeof(int));
        memcpyToGPU(atom->d_atom.border_map, atom->border_map, NmaxGhost * sizeof(int));
        cuda_assert("updatePbc.reneigh", cudaPeekAtLastError());
        cuda_assert("updatePbc.reneigh", cudaDeviceSynchronize());
    }

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
    memcpyFromGPU(atom->x, atom->d_atom.x, sizeof(MD_FLOAT) * atom->Nlocal * 3);
}
}
