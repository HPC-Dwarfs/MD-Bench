/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <force.h>
#include <parameter.h>
#include <stdio.h>
#include <stdlib.h>
#include <util.h>

ComputeForceFunction computeForce;

void initForce(Parameter* param)
{
    switch (param->force_field) {
    case FF_EAM:
        computeForce = computeForceEam;
        break;
    case FF_LJ:
#if (COULOMB_RULE != COULOMB_NONE || USE_EXCLUSIONS) &&                                  \
    (defined(CLUSTERPAIR_KERNEL_REF) || defined(CUDA_TARGET))
        fprintf(stderr,
            "Error: electrostatics and exclusions are only supported by the SIMD "
            "kernels of the CPU (4xN, 2xNN, 2xN)!\n");
        exit(EXIT_FAILURE);
#endif
#if defined(CLUSTERPAIR_KERNEL_REF)
        computeForce = computeForceLJRef;
#elif defined(CLUSTERPAIR_KERNEL_4XN)
        if (param->half_neigh || param->method) {
            computeForce = computeForceLJ4xnHalfNeigh;
        } else {
            computeForce = computeForceLJ4xnFullNeigh;
        }
#elif defined(CLUSTERPAIR_KERNEL_2XNN)
        if (param->half_neigh || param->method) {
            computeForce = computeForceLJ2xnnHalfNeigh;
        } else {
            computeForce = computeForceLJ2xnnFullNeigh;
        }
#elif defined(CLUSTERPAIR_KERNEL_2XN)
        if (param->half_neigh) {
            computeForce = computeForceLJ2xnHalfNeigh;
        } else {
            computeForce = computeForceLJ2xnFullNeigh;
        }
#elif defined(CLUSTERPAIR_KERNEL_GPU)
        if (param->super_clustering) {
            computeForce = computeForceLJCudaSup;
        } else {
            computeForce = computeForceLJCuda;
        }
#endif
#if !defined(_MPI) && !defined(CUDA_TARGET) &&                                           \
    (defined(CLUSTERPAIR_KERNEL_4XN) || defined(CLUSTERPAIR_KERNEL_2XNN) ||              \
        defined(CLUSTERPAIR_KERNEL_2XN))
        // Without MPI the periodic images are ghost clusters of the same process. A half
        // neighbor list then keeps only half of them (half shell) and the kernels write
        // the reaction forces into the ghost clusters, which reverse() folds back. With
        // the plain half list every pair across a periodic boundary would be computed
        // from both sides.
        if (param->half_neigh && param->method == fullShell) {
            param->method = halfShell;
        }
#endif
        break;
    case FF_LJ_TABLE:
// #if defined(CLUSTERPAIR_KERNEL_REF)
#if defined(CLUSTERPAIR_KERNEL_REF)
        computeForce = computeForceLJTableRef;
#else
        fprintf(stderr,
            "Error: Tabulated force fields for SIMD/GPU cases not yet implemented!\n");
        exit(EXIT_FAILURE);
#endif
        /*
        #elif defined(CLUSTERPAIR_KERNEL_4XN)
                if (param->half_neigh || param->method) {
                    computeForce = computeForceLJTable4xnHalfNeigh;
                } else {
                    computeForce = computeForceLJTable4xnFullNeigh;
                }
        #elif defined(CLUSTERPAIR_KERNEL_2XNN)
                if (param->half_neigh || param->method) {
                    computeForce = computeForceLJTable2xnnHalfNeigh;
                } else {
                    computeForce = computeForceLJTable2xnnFullNeigh;
                }
        #elif defined(CLUSTERPAIR_KERNEL_GPU)
                if (param->super_clustering) {
                    computeForce = computeForceLJTableCudaSup;
                } else {
                    computeForce = computeForceLJTableCuda;
                }
        #endif
        */
        break;
    default:
        fprintf(stderr, "Error: Unknown force field!\n");
        exit(EXIT_FAILURE);
    }
}
