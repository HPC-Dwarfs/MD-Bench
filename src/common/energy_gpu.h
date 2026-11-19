/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// Energies (-energy) of the GPU kernels: a device buffer of the LJ and Coulomb sums
#ifndef __ENERGY_GPU_H_
#define __ENERGY_GPU_H_

#include <device.h>
#include <energy.h>
#include <parameter.h>

// Cleared device buffer for the energy sums of a force call (NULL without -energy)
static inline double* deviceEnergies(Parameter* param)
{
    static double* energies = NULL;
    if (!param->compute_energy) {
        return NULL;
    }
    if (energies == NULL) {
        energies = (double*)allocateGPU(2 * sizeof(double));
    }
    memsetGPU(energies, 0, 2 * sizeof(double));
    return energies;
}

// Stores the device energies times w (1/2 for full lists) in param
static inline void storeDeviceEnergies(Parameter* param,
    const int nlocal,
    const int* type,
    const MD_FLOAT* q_per_type,
    double* energies,
    const double w)
{
    if (energies == NULL) {
        return;
    }
    double e[2];
    memcpyFromGPU(e, energies, 2 * sizeof(double));
    storeEnergies(param,
        w * e[0],
        w * e[1],
        localChargeSquaredSum(nlocal, type, q_per_type));
}

// Adds the thread energies to energies[0..1], one atomic per warp when the warp is full
static __device__ __forceinline__ void warpEnergiesAdd(
    double vlj, double vc, double* energies)
{
#ifdef __HIP_PLATFORM_AMD__
    const unsigned long long full = 0xffffffffffffffffULL;
    const unsigned long long mask = __activemask();
#else
    const unsigned int full = 0xffffffffU;
    const unsigned int mask = __activemask();
#endif
    if (mask == full) {
        for (int offset = warpSize / 2; offset > 0; offset >>= 1) {
            vlj += __shfl_down_sync(full, vlj, offset);
            vc += __shfl_down_sync(full, vc, offset);
        }
        if (((threadIdx.y * blockDim.x + threadIdx.x) % warpSize) == 0) {
            atomicAdd(&energies[0], vlj);
            atomicAdd(&energies[1], vc);
        }
    } else {
        atomicAdd(&energies[0], vlj);
        atomicAdd(&energies[1], vc);
    }
}

#endif
