/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
//---
#include <likwid-marker.h>

#include <device.h>
#include <gpu_profiler.h>

extern "C" {
#include <atom.h>
#include <force.h>
#include <ljtable.h>
#include <neighbor.h>
#include <parameter.h>
#include <timing.h>
#include <util.h>
}

#include <kernel_select_gpu.h>

// Tabulated/spline-interpolated Lennard-Jones forces on the GPU (GROMACS-style).
// Mirrors computeForceLJCUDA in forceCuda.cu but replaces the analytic force
// arithmetic with a cubic-spline lookup of two universal force-shape functions
// (see src/common/ljtable.{h,c}). The coefficient table is uploaded to the
// device once and kept resident in ljtable.d_coeff.
// RULE: LJ combination rule (cutforcesq, sigma6 and epsilon are only used by single)

template <int RULE>
__global__ void computeForceLJTableCudaFullNeigh(DeviceAtom a,
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    int Nlocal,
    DeviceNeighbor neigh,
    int* neigh_neighbors,
    int* neigh_numneigh,
    int ntypes,
    int nm1,
    MD_FLOAT inv_h,
    const MD_FLOAT* coeff)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom         = &a;
    const int numneighs      = neigh_numneigh[i];

    MD_FLOAT xtmp = atom_x(i);
    MD_FLOAT ytmp = atom_y(i);
    MD_FLOAT ztmp = atom_z(i);

    MD_FLOAT fix = 0;
    MD_FLOAT fiy = 0;
    MD_FLOAT fiz = 0;

    MD_FLOAT sqrt_eps_i = (MD_FLOAT)0.0;
    MD_FLOAT sigma3_i   = (MD_FLOAT)0.0;
    int type_i          = 0;
    if constexpr (RULE == LJ_COMB_GEOM) {
        sqrt_eps_i = atom->sqrt_epsilon[i];
        sigma3_i   = atom->sigma3[i];
    } else if constexpr (RULE == LJ_COMB_FULL) {
        type_i = atom->type[i];
    }

    for (int k = 0; k < numneighs; k++) {
        int j         = neighs(neigh_neighbors, i, k, Nlocal, neighbor);
        MD_FLOAT delx = xtmp - atom_x(j);
        MD_FLOAT dely = ytmp - atom_y(j);
        MD_FLOAT delz = ztmp - atom_z(j);
        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

        MD_FLOAT cutforcesq_ij = cutforcesq;
        MD_FLOAT sigma6_ij     = sigma6;
        MD_FLOAT epsilon_ij    = epsilon;
        if constexpr (RULE == LJ_COMB_GEOM) {
            sigma6_ij  = sigma3_i * atom->sigma3[j];
            epsilon_ij = sqrt_eps_i * atom->sqrt_epsilon[j];
        } else if constexpr (RULE == LJ_COMB_FULL) {
            const int type_j  = atom->type[j];
            const int type_ij = type_i * ntypes + type_j;
            cutforcesq_ij     = atom->cutforcesq[type_ij];
            sigma6_ij         = atom->sigma6[type_ij];
            epsilon_ij        = atom->epsilon[type_ij];
        }

        if (rsq < cutforcesq_ij) {
            MD_FLOAT u = LJ_TABLE_COORD(rsq) * inv_h;
            int m      = (int)u;
            if (m > nm1) m = nm1;
            MD_FLOAT eps       = u - (MD_FLOAT)m;
            const MD_FLOAT* cc = &coeff[m * LJ_TABLE_STRIDE];
            MD_FLOAT hrep      = cc[0] + eps * (cc[1] + eps * (cc[2] + eps * cc[3]));
            MD_FLOAT gdisp     = cc[4] + eps * (cc[5] + eps * (cc[6] + eps * cc[7]));
            MD_FLOAT force     = epsilon_ij * sigma6_ij * (sigma6_ij * hrep + gdisp);
            fix += delx * force;
            fiy += dely * force;
            fiz += delz * force;
        }
    }

    atom_fx(i) = fix;
    atom_fy(i) = fiy;
    atom_fz(i) = fiz;
}

template <int RULE>
__global__ void computeForceLJTableCudaHalfNeigh(DeviceAtom a,
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    int Nlocal,
    DeviceNeighbor neigh,
    int* neigh_neighbors,
    int* neigh_numneigh,
    int ntypes,
    int nm1,
    MD_FLOAT inv_h,
    const MD_FLOAT* coeff)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom         = &a;
    const int numneighs      = neigh_numneigh[i];

    MD_FLOAT xtmp = atom_x(i);
    MD_FLOAT ytmp = atom_y(i);
    MD_FLOAT ztmp = atom_z(i);

    MD_FLOAT fix = 0;
    MD_FLOAT fiy = 0;
    MD_FLOAT fiz = 0;

    MD_FLOAT sqrt_eps_i = (MD_FLOAT)0.0;
    MD_FLOAT sigma3_i   = (MD_FLOAT)0.0;
    int type_i          = 0;
    if constexpr (RULE == LJ_COMB_GEOM) {
        sqrt_eps_i = atom->sqrt_epsilon[i];
        sigma3_i   = atom->sigma3[i];
    } else if constexpr (RULE == LJ_COMB_FULL) {
        type_i = atom->type[i];
    }

    for (int k = 0; k < numneighs; k++) {
        int j         = neighs(neigh_neighbors, i, k, Nlocal, neighbor);
        MD_FLOAT delx = xtmp - atom_x(j);
        MD_FLOAT dely = ytmp - atom_y(j);
        MD_FLOAT delz = ztmp - atom_z(j);
        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

        MD_FLOAT cutforcesq_ij = cutforcesq;
        MD_FLOAT sigma6_ij     = sigma6;
        MD_FLOAT epsilon_ij    = epsilon;
        if constexpr (RULE == LJ_COMB_GEOM) {
            sigma6_ij  = sigma3_i * atom->sigma3[j];
            epsilon_ij = sqrt_eps_i * atom->sqrt_epsilon[j];
        } else if constexpr (RULE == LJ_COMB_FULL) {
            const int type_j  = atom->type[j];
            const int type_ij = type_i * ntypes + type_j;
            cutforcesq_ij     = atom->cutforcesq[type_ij];
            sigma6_ij         = atom->sigma6[type_ij];
            epsilon_ij        = atom->epsilon[type_ij];
        }

        if (rsq < cutforcesq_ij) {
            MD_FLOAT u = LJ_TABLE_COORD(rsq) * inv_h;
            int m      = (int)u;
            if (m > nm1) m = nm1;
            MD_FLOAT eps       = u - (MD_FLOAT)m;
            const MD_FLOAT* cc = &coeff[m * LJ_TABLE_STRIDE];
            MD_FLOAT hrep      = cc[0] + eps * (cc[1] + eps * (cc[2] + eps * cc[3]));
            MD_FLOAT gdisp     = cc[4] + eps * (cc[5] + eps * (cc[6] + eps * cc[7]));
            MD_FLOAT force     = epsilon_ij * sigma6_ij * (sigma6_ij * hrep + gdisp);
            MD_FLOAT partial_force_x = delx * force;
            MD_FLOAT partial_force_y = dely * force;
            MD_FLOAT partial_force_z = delz * force;

            atomicAdd(&atom_fx(j), -partial_force_x);
            atomicAdd(&atom_fy(j), -partial_force_y);
            atomicAdd(&atom_fz(j), -partial_force_z);

            fix += partial_force_x;
            fiy += partial_force_y;
            fiz += partial_force_z;
        }
    }

    atomicAdd(&atom_fx(i), fix);
    atomicAdd(&atom_fy(i), fiy);
    atomicAdd(&atom_fz(i), fiz);
}

extern "C" {

// Upload the host-built coefficient table to the device on first use.
static void ensureLJTableOnGPU(void)
{
    if (ljtable.d_coeff != NULL) {
        return;
    }
    size_t bytes    = (size_t)(ljtable.n + 1) * LJ_TABLE_STRIDE * sizeof(MD_FLOAT);
    ljtable.d_coeff = (MD_FLOAT*)allocateGPU(bytes);
    memcpyToGPU(ljtable.d_coeff, ljtable.coeff, bytes);
}

// Release the device-resident table (counterpart to the lazy upload above).
void freeLJTableGPU(void)
{
    if (ljtable.d_coeff != NULL) {
        GPUfree(ljtable.d_coeff);
        ljtable.d_coeff = NULL;
    }
}

double computeForceLJTableCUDA(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJTableCUDA begin\n");

    const int num_threads_per_block = get_cuda_num_threads();
    int Nlocal                      = atom->Nlocal;
    int Nmax                        = atom->Nmax;
    MD_FLOAT cutforcesq             = param->cutforce * param->cutforce;
    MD_FLOAT sigma6                 = param->sigma6;
    MD_FLOAT epsilon                = param->epsilon;

    ensureLJTableOnGPU();
    const int nm1         = ljtable.n - 1;
    const MD_FLOAT inv_h  = ljtable.inv_h;
    const MD_FLOAT* coeff = ljtable.d_coeff;

    GPU_PROFILE_START("force_lj_table");
    const int num_blocks = ceil((float)Nlocal / (float)num_threads_per_block);
    double S             = getTimeStamp();
    LIKWID_MARKER_START("force");

    if (neighbor->half_neigh) {
#ifdef ATOM_POSITION_AOS
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax * 3);
#else
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fy, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fz, 0, sizeof(MD_FLOAT) * Nmax);
#endif
        auto kernel = SELECT_GPU_KERNEL_RULE(computeForceLJTableCudaHalfNeigh,
            param->lj_comb_rule);
        kernel<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
            cutforcesq,
            sigma6,
            epsilon,
            Nlocal,
            neighbor->d_neighbor,
            neighbor->d_neighbor.neighbors,
            neighbor->d_neighbor.numneigh,
            atom->ntypes,
            nm1,
            inv_h,
            coeff);
    } else {
        auto kernel = SELECT_GPU_KERNEL_RULE(computeForceLJTableCudaFullNeigh,
            param->lj_comb_rule);
        kernel<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
            cutforcesq,
            sigma6,
            epsilon,
            Nlocal,
            neighbor->d_neighbor,
            neighbor->d_neighbor.neighbors,
            neighbor->d_neighbor.numneigh,
            atom->ntypes,
            nm1,
            inv_h,
            coeff);
    }

    cuda_assert("computeForceLJTableCuda", cudaPeekAtLastError());
    cuda_assert("computeForceLJTableCuda", cudaDeviceSynchronize());
    GPU_PROFILE_STOP();

    LIKWID_MARKER_STOP("force");
    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJTableCUDA end\n");
    return E - S;
}
}
