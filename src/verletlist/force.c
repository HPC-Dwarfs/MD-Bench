/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include "parameter.h"
#include <force.h>
#include <stdio.h>
#include <stdlib.h>

ComputeForceFunction computeForce;

void initForce(Parameter* param)
{
#if COULOMB_RULE != COULOMB_NONE
    // Electrostatics are only in the plain C (force_lj.c) and GPU kernels
    if (param->force_field != FF_LJ) {
        fprintf(stderr, "Error: electrostatics are only supported with -f lj!\n");
        exit(EXIT_FAILURE);
    }
#ifdef __SIMD_KERNEL__
    fprintf(stderr,
        "Error: electrostatics are not supported by the verlet-list "
        "SIMD kernels, build with USE_SIMD_KERNEL=false!\n");
    exit(EXIT_FAILURE);
#endif
#endif
#ifdef CUDA_TARGET
    // The GPU neighbor lists keep all ghosts and the reaction forces of ghosts are not
    // folded back on the GPU
    if (param->method) {
        fprintf(stderr, "Error: the verlet-list GPU kernels do not support -method!\n");
        exit(EXIT_FAILURE);
    }
#endif
    // Energies (-energy) are computed by the LJ kernels only
    if (param->compute_energy && param->force_field != FF_LJ) {
        fprintf(stderr, "Error: -energy is only supported with -f lj!\n");
        exit(EXIT_FAILURE);
    }
#ifdef __SIMD_KERNEL__
    if (param->compute_energy) {
        fprintf(stderr,
            "Error: -energy is not supported by the verlet-list SIMD kernels, build with "
            "USE_SIMD_KERNEL=false!\n");
        exit(EXIT_FAILURE);
    }
#endif
    switch (param->force_field) {
    case FF_EAM:
        computeForce = computeForceEam;
        break;
    case FF_LJ_TABLE:
#ifdef CUDA_TARGET
        computeForce = computeForceLJTableCUDA;
#else
        if (param->half_neigh || param->method) {
            computeForce = computeForceLJTableHalfNeigh;
        } else {
            computeForce = computeForceLJTableFullNeigh;
        }
#endif
        break;
    case FF_LJ:
#ifdef CUDA_TARGET
        computeForce = computeForceLJCUDA;
#else
#ifdef __SIMD_KERNEL__
#ifdef __SIMD_COMPRESS__
        if (param->half_neigh || param->method) {
            computeForce = param->accurate_reciprocal
                               ? computeForceLJHalfNeigh_simd_compressAccurate
                               : computeForceLJHalfNeigh_simd_compress;
        } else {
            computeForce = param->accurate_reciprocal
                               ? computeForceLJFullNeigh_simd_compressAccurate
                               : computeForceLJFullNeigh_simd_compress;
        }
#elif defined(__SIMD_VLA__)
        if (param->half_neigh || param->method) {
            computeForce = param->accurate_reciprocal
                               ? computeForceLJHalfNeigh_simd_vlaAccurate
                               : computeForceLJHalfNeigh_simd_vla;
        } else {
            computeForce = param->accurate_reciprocal
                               ? computeForceLJFullNeigh_simd_vlaAccurate
                               : computeForceLJFullNeigh_simd_vla;
        }
#else
        if (param->half_neigh || param->method) {
            computeForce = param->accurate_reciprocal
                               ? computeForceLJHalfNeigh_simdAccurate
                               : computeForceLJHalfNeigh_simd;
        } else {
            computeForce = param->accurate_reciprocal
                               ? computeForceLJFullNeigh_simdAccurate
                               : computeForceLJFullNeigh_simd;
        }
#endif
#else
        // The copy with molecular exclusions for systems with molecules
        // (kernel_variant.h)
        if (param->half_neigh || param->method) {
            computeForce = param->molsize > 1 ? computeForceLJHalfNeighExcl
                                              : computeForceLJHalfNeigh;
        } else {
            computeForce = param->molsize > 1 ? computeForceLJFullNeighExcl
                                              : computeForceLJFullNeigh;
        }
#endif
#endif
        break;
    default:
        fprintf(stderr, "Error: Unknown force field!\n");
        exit(EXIT_FAILURE);
    }
}
