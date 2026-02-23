/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <force.h>
#include <parameter.h>
#include <stdlib.h>

ComputeForceFunction computeForce;

void initForce(Parameter* param)
{
    switch (param->force_field) {
    case FF_EAM:
        computeForce = computeForceEam;
        break;
    case FF_LJ:
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
#elif defined(CLUSTERPAIR_KERNEL_GPU)
        if (param->super_clustering) {
            computeForce = computeForceLJCudaSup;
        } else {
            computeForce = computeForceLJCuda;
        }
#endif
    }
}
