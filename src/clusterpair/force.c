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
#if COULOMB_RULE != COULOMB_NONE
    // The EAM and tabulated kernels have no electrostatics
    if (param->force_field != FF_LJ) {
        fprintf(stderr, "Error: electrostatics are only supported with -f lj!\n");
        exit(EXIT_FAILURE);
    }
#endif
    // Energies (-energy) are computed by the LJ kernels only
    if (param->compute_energy && param->force_field != FF_LJ) {
        fprintf(stderr, "Error: -energy is only supported with -f lj!\n");
        exit(EXIT_FAILURE);
    }
    switch (param->force_field) {
    case FF_EAM:
        computeForce = computeForceEam;
        break;
    case FF_LJ:
#if defined(CLUSTERPAIR_KERNEL_REF)
        computeForce = SELECT_RULE_VARIANT(param, computeForceLJRef);
#elif defined(CLUSTERPAIR_KERNEL_4XN)
        if (param->half_neigh || param->method) {
            computeForce = SELECT_SIMD_VARIANT(param, computeForceLJ4xnHalfNeigh);
        } else {
            computeForce = SELECT_SIMD_VARIANT(param, computeForceLJ4xnFullNeigh);
        }
#elif defined(CLUSTERPAIR_KERNEL_2XNN)
        if (param->half_neigh || param->method) {
            computeForce = SELECT_SIMD_VARIANT(param, computeForceLJ2xnnHalfNeigh);
        } else {
            computeForce = SELECT_SIMD_VARIANT(param, computeForceLJ2xnnFullNeigh);
        }
#elif defined(CLUSTERPAIR_KERNEL_2XN)
        if (param->half_neigh) {
            computeForce = SELECT_SIMD_VARIANT(param, computeForceLJ2xnHalfNeigh);
        } else {
            computeForce = SELECT_SIMD_VARIANT(param, computeForceLJ2xnFullNeigh);
        }
#elif defined(CLUSTERPAIR_KERNEL_GPU)
        if (param->super_clustering) {
            computeForce = computeForceLJCudaSup;
        } else {
            computeForce = computeForceLJCuda;
        }
#endif
        break;
    case FF_LJ_TABLE:
// #if defined(CLUSTERPAIR_KERNEL_REF)
#if defined(CLUSTERPAIR_KERNEL_REF)
        computeForce = SELECT_RULE_VARIANT(param, computeForceLJTableRef);
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

#if !defined(_MPI) && !defined(CUDA_TARGET)
    // Without MPI the periodic images are ghost clusters of the same process. A half
    // neighbor list then keeps only half of them (half shell) and the kernels write
    // the reaction forces into the ghost clusters, which reverse() folds back. With
    // the plain half list every pair across a periodic boundary would be computed
    // from both sides. This holds for every CPU kernel, so it must be decided here
    // and before initNeighbor(), which reads param->method.
    if (param->half_neigh && param->method == fullShell) {
        param->method = halfShell;
    }
#endif
}
