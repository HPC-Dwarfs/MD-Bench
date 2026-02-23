/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <comm.h>
#include <limits.h>
#include <math.h>
#include <parameter.h>
#include <pbc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <timing.h>
#include <unistd.h>
#include <util.h>

static void addDummyCluster(Parameter*, Atom*);
#ifdef CUDA_TARGET
extern void copyGhostFromGPU(Atom*);
extern void copyGhostToGPU(Atom*);
extern void copyForceFromGPU(Atom*);
extern void copyForceToGPU(Atom*);
#endif 

double forward(Comm* comm, Atom* atom, Parameter* param) {
    double S, E;
    S = getTimeStamp();
#ifdef _MPI

#ifdef CUDA_TARGET
    copyGhostFromGPU(atom); 
#endif

    if (param->method == halfShell) {
        for (int iswap = 0; iswap < 5; iswap++)
            forwardComm(comm, atom, iswap);
    } else if (param->method == eightShell) {
        for (int iswap = 0; iswap < 6; iswap += 2)
            forwardComm(comm, atom, iswap);
    } else {
        for (int iswap = 0; iswap < 6; iswap++)
            forwardComm(comm, atom, iswap);
    }
#ifdef CUDA_TARGET
    copyGhostToGPU(atom);
#endif

#else 
    updatePbc(atom, param, false);
#endif
    E = getTimeStamp();
    return E - S;
}

double reverse(Comm* comm, Atom* atom, Parameter* param) {
    double S, E;
    S = getTimeStamp();
#ifdef _MPI

#ifdef CUDA_TARGET
    copyForceFromGPU(atom); 
#endif

    if (param->method == halfShell) {
        for (int iswap = 4; iswap >= 0; iswap--)
            reverseComm(comm, atom, iswap);
    } else if (param->method == eightShell) {
        for (int iswap = 4; iswap >= 0; iswap -= 2)
            reverseComm(comm, atom, iswap);
    } else if (param->method == halfStencil) {
        for (int iswap = 5; iswap >= 0; iswap--)
            reverseComm(comm, atom, iswap);
    } else {
    } // Full Shell Reverse does nothing

#ifdef CUDA_TARGET
    copyForceToGPU(atom);
#endif

#endif
    E = getTimeStamp();
    return E - S;
}

#ifdef _MPI
void ghostNeighbor(Comm* comm, Atom* atom, Parameter* param) {
#ifdef CLUSTER_PAIR
    atom->Nclusters_ghost = 0;
#endif
    atom->Nghost = 0;
    if (param->method == halfShell) {
        for (int iswap = 0; iswap < 5; iswap++)
            ghostComm(comm, param, atom, iswap);
    } else if (param->method == eightShell) {
        for (int iswap = 0; iswap < 6; iswap += 2)
            ghostComm(comm, param, atom, iswap);
    } else {
        for (int iswap = 0; iswap < 6; iswap++)
            ghostComm(comm, param, atom, iswap);
    }
#ifdef CLUSTER_PAIR    
    addDummyCluster(param, atom);    
#endif
}
#endif

#ifdef CLUSTER_PAIR

void addDummyCluster(Parameter* param, Atom* atom) {
    // atom->Nclusters_ghost++; // GHOST J CLUSTERS
    // atom->Nclusters = atom->Nclusters_local + Nghost + 1;
    atom->dummy_cj = LOCAL + GHOST;

    if ((LOCAL + GHOST) * JFAC >= atom->Nclusters_max) {
        growClusters(atom, param->super_clustering);
    }

    // Add dummy cluster at the end
    int cjVecBase = CJ_VECTOR_BASE_INDEX(atom->dummy_cj);
    MD_FLOAT* cjX = &atom->cl_x[cjVecBase];

    for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
        cjX[CL_X_INDEX(cjj)] = INFINITY;
        cjX[CL_Y_INDEX(cjj)] = INFINITY;
        cjX[CL_Z_INDEX(cjj)] = INFINITY;
    } 
}
#endif
