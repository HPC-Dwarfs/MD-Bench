/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdbool.h>
#include <stdlib.h>

#include <atom.h>
#include <force.h>
#include <integrate.h>
#include <parameter.h>
#include <util.h>

#ifdef CUDA_TARGET
#ifdef CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
IntegrationFunction initialIntegrate = cudaInitialIntegrateSup;
IntegrationFunction finalIntegrate   = cudaFinalIntegrateSup;
#else
IntegrationFunction initialIntegrate = cudaInitialIntegrate;
IntegrationFunction finalIntegrate   = cudaFinalIntegrate;
#endif
#else
IntegrationFunction initialIntegrate = initialIntegrateCPU;
IntegrationFunction finalIntegrate   = finalIntegrateCPU;
#endif

void initialIntegrateCPU(Parameter* param, Atom* atom) {
    DEBUG_MESSAGE("cpuInitialIntegrate start\n");

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int ciVecBase = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ciX = &atom->cl_x[ciVecBase];
        MD_FLOAT* ciV = &atom->cl_v[ciVecBase];
        MD_FLOAT* ciF = &atom->cl_f[ciVecBase];

        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ciV[CL_X_INDEX_3D(cii)] += param->dtforce * ciF[CL_X_INDEX_3D(cii)];
            ciV[CL_Y_INDEX_3D(cii)] += param->dtforce * ciF[CL_Y_INDEX_3D(cii)];
            ciV[CL_Z_INDEX_3D(cii)] += param->dtforce * ciF[CL_Z_INDEX_3D(cii)];
            ciX[CL_X_INDEX_3D(cii)] += param->dt * ciV[CL_X_INDEX_3D(cii)];
            ciX[CL_Y_INDEX_3D(cii)] += param->dt * ciV[CL_Y_INDEX_3D(cii)];
            ciX[CL_Z_INDEX_3D(cii)] += param->dt * ciV[CL_Z_INDEX_3D(cii)];
        }

        /*
        // Check if there is an invalid cluster
        MD_FLOAT threshold = 10000.0;
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            if(
                ciX[CL_X_INDEX_3D(cii)] < -threshold || ciX[CL_X_INDEX_3D(cii)] > threshold ||
                ciX[CL_Y_INDEX_3D(cii)] < -threshold || ciX[CL_Y_INDEX_3D(cii)] > threshold ||
                ciX[CL_Z_INDEX_3D(cii)] < -threshold || ciX[CL_Z_INDEX_3D(cii)] > threshold) {
                fprintf(stdout, "INVALID CLUSTER: %d\n", ci);
                fprintf(stdout, "%f, %f, %f\n", ciX[CL_X_INDEX_3D + cii], ciX[CL_Y_INDEX_3D(cii)], ciX[CL_Z_INDEX_3D(cii)]);
                exit(-1);
                break;
            }
        }
        */
    }

    DEBUG_MESSAGE("cpuInitialIntegrate end\n");
}

void finalIntegrateCPU(Parameter* param, Atom* atom) {
    DEBUG_MESSAGE("cpuFinalIntegrate start\n");

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        int ciVecBase = CI_VECTOR3_BASE_INDEX(ci);
        MD_FLOAT* ciV = &atom->cl_v[ciVecBase];
        MD_FLOAT* ciF = &atom->cl_f[ciVecBase];

        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ciV[CL_X_INDEX_3D(cii)] += param->dtforce * ciF[CL_X_INDEX_3D(cii)];
            ciV[CL_Y_INDEX_3D(cii)] += param->dtforce * ciF[CL_Y_INDEX_3D(cii)];
            ciV[CL_Z_INDEX_3D(cii)] += param->dtforce * ciF[CL_Z_INDEX_3D(cii)];
        }
    }

    DEBUG_MESSAGE("cpuFinalIntegrate end\n");
}
