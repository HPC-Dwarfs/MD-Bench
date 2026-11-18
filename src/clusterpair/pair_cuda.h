/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Loads of the per-atom electrostatics and exclusion data for the GPU cluster-pair
 * kernels, which evaluate each pair with pairForceScalar() (common/pair_scalar.h). */
#ifndef __PAIR_CUDA_H_
#define __PAIR_CUDA_H_

#include <atom.h>
#include <pair_scalar.h>
#include <parameter.h>

/* Charge of the cluster slot sca, 0 without electrostatics (cl_q is then not on the GPU) */
static __device__ __forceinline__ MD_FLOAT chargeCuda(const MD_FLOAT* cl_q, int sca)
{
#if COULOMB_RULE != COULOMB_NONE
    return cl_q[sca];
#else
    (void)cl_q;
    (void)sca;
    return (MD_FLOAT)0.0;
#endif
}

/* Molecule of the cluster slot sca, 0 without exclusions (cl_mol is then not copied) */
template <bool EXCL>
static __device__ __forceinline__ int moleculeCuda(const int* cl_mol, int sca)
{
    return EXCL ? cl_mol[sca] : 0;
}

// Block sum of the energies (NT threads, all must call it) added to energies[0..1]
template <int NT>
static __device__ __forceinline__ void blockEnergiesAdd(
    double vlj, double vc, double* energies)
{
    __shared__ double s_e[2][NT];
    const int t = threadIdx.y * blockDim.x + threadIdx.x;
    s_e[0][t]   = vlj;
    s_e[1][t]   = vc;
    __syncthreads();
    for (int stride = NT / 2; stride > 0; stride >>= 1) {
        if (t < stride) {
            s_e[0][t] += s_e[0][t + stride];
            s_e[1][t] += s_e[1][t + stride];
        }
        __syncthreads();
    }
    if (t == 0) {
        atomicAdd(&energies[0], s_e[0][0]);
        atomicAdd(&energies[1], s_e[1][0]);
    }
}

#endif
