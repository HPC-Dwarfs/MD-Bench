/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <eam.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>

#ifndef __FORCE_H_
#define __FORCE_H_

typedef double (*ComputeForceFunction)(Parameter*, Atom*, Neighbor*, Stats*);
extern ComputeForceFunction computeForce;

enum forcetype { FF_LJ = 0, FF_EAM };

extern void initForce(Parameter*);
extern double computeForceLJRef(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJ4xnHalfNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJ4xnFullNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJ2xnnHalfNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJ2xnnFullNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceEam(Parameter*, Atom*, Neighbor*, Stats*);

#ifdef CUDA_TARGET
double computeForceLJCuda(Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats);
double computeForceLJCudaSup(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats);
#endif

// Nbnxn layouts (as of GROMACS):
// Simd4xN: M=4, N=VECTOR_WIDTH
// Simd2xNN: M=4, N=(VECTOR_WIDTH/2)
// Cuda: M=8, N=VECTOR_WIDTH

/* Comments from GROMACS:
 *
 * We need to choose if we want 2x(N+N) or 4xN kernels.
 * This can be controlled through CLUSTER_PAIR_KERNEL option:
 * - auto: Automatically choose based on SIMD acceleration and CPU info
 * - 4xN: Force 4xN kernel layout
 * - 2xNN: Force 2xNN kernel layout
 *
 * Auto selection behavior:
 * 4xN calculates more (zero) interactions, but has less pair-search
 * work and much better kernel instruction scheduling.
 *
 * Up till now we have only seen that on Intel Sandy/Ivy Bridge,
 * which doesn't have FMA, both the analytical and tabulated Ewald
 * kernels have similar pair rates for 4x8 and 2x(4+4), so we choose
 * 2x(4+4) because it results in significantly fewer pairs.
 * For RF, the raw pair rate of the 4x8 kernel is higher than 2x(4+4),
 * 10% with HT, 50% without HT. As we currently don't detect the actual
 * use of HT, use 4x8 to avoid a potential performance hit.
 * On Intel Haswell 4x8 is always faster.
 *
 *
 * The nbnxn SIMD 4xN and 2x(N+N) kernels can be added independently.
 * Currently the 2xNN SIMD kernels only make sense with:
 *  8-way SIMD: 4x4 setup, performance wise only useful on CPUs without FMA or on AMD Zen1
 * 16-way SIMD: 4x8 setup, used in single precision with 512 bit wide SIMD
 */

#ifdef CUDA_TARGET
extern double computeForceLJCUDA(Parameter*, Atom*, Neighbor*, Stats*);
#undef VECTOR_WIDTH
#define VECTOR_WIDTH 8
#define CLUSTERPAIR_KERNEL_GPU
#ifndef CLUSTERPAIR_KERNEL_GPU_SIMPLE
#define CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
#endif
#define KERNEL_NAME "GPU"
#define CLUSTER_M   8
#define CLUSTER_N   VECTOR_WIDTH
#define UNROLL_J    1
#else
#ifdef USE_REFERENCE_KERNEL
#define CLUSTERPAIR_KERNEL_REF
#define KERNEL_NAME "Reference"
#define CLUSTER_M   1
#define CLUSTER_N   VECTOR_WIDTH
#else
#define CLUSTER_M 4

// Auto selection based on VECTOR_WIDTH and architecture
#ifdef CLUSTERPAIR_KERNEL_AUTO
#if (VECTOR_WIDTH > (CLUSTER_M * 2))
#define CLUSTERPAIR_KERNEL_2XNN
#else
#define CLUSTERPAIR_KERNEL_4XN
#endif
#endif

// Define the kernel-specific macros based on which kernel is selected
#ifdef CLUSTERPAIR_KERNEL_4XN
#define KERNEL_NAME "Simd4xN"
#define CLUSTER_N   VECTOR_WIDTH
#define UNROLL_I    4
#define UNROLL_J    1
#endif

#ifdef CLUSTERPAIR_KERNEL_2XNN
#define KERNEL_NAME "Simd2xNN"
#define CLUSTER_N   (VECTOR_WIDTH / 2)
#define UNROLL_I    4
#define UNROLL_J    2
#endif

// Verify that one of the kernel variants is selected
#if !defined(CLUSTERPAIR_KERNEL_4XN) && !defined(CLUSTERPAIR_KERNEL_2XNN)
#error "No cluster pair kernel variant selected"
#endif

#endif
#endif

#if CLUSTER_N != 2 && CLUSTER_N != 4 && CLUSTER_N != 8
#error "Cluster N dimension can be only 2, 4 and 8"
#endif

// Super-clustering macros
#if defined(CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS) && CLUSTER_M != CLUSTER_N
#error "For super-clusters, M must be equal to N"
#endif

#ifdef CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
#define SCLUSTER_SIZE_X          2
#define SCLUSTER_SIZE_Y          2
#define SCLUSTER_SIZE_Z          2
#define SCLUSTER_SIZE            (SCLUSTER_SIZE_X * SCLUSTER_SIZE_Y * SCLUSTER_SIZE_Z)
#define SCI_BASE_INDEX(a, b)     ((a) * CLUSTER_N * SCLUSTER_SIZE * (b))
#define SCI_FROM_CJ(a)           ((a) / SCLUSTER_SIZE)
#ifdef POSITION_AOS4_SUP
#define ATOM_DIM 4
#else
#define ATOM_DIM 3
#endif
#else
#define SCLUSTER_SIZE_X          1
#define SCLUSTER_SIZE_Y          1
#define SCLUSTER_SIZE_Z          1
#define SCLUSTER_SIZE            1
#define SCI_BASE_INDEX(a, b)     ((a) * CLUSTER_M * (b))
#define SCI_FROM_CJ(a)           (a)
#define ATOM_DIM 3
#endif


#define SCI_SCALAR_BASE_INDEX(a)  (SCI_BASE_INDEX(a, 1))
#define SCI_VECTOR3_BASE_INDEX(a) (SCI_BASE_INDEX(a, 3))
#define SCI_VECTOR_BASE_INDEX(a)  (SCI_BASE_INDEX(a, ATOM_DIM))

#if defined(CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS)
#define CJ0_FROM_CI(a)      (a)
#define CJ1_FROM_CI(a)      (a)
#define CI_BASE_INDEX(a, b) ((a)*CLUSTER_N * (b))
#ifdef SOA_SUP
#define CJ_BASE_INDEX(a, b) ((((a) / SCLUSTER_SIZE) * SCLUSTER_SIZE * CLUSTER_N * (b)) + \
                             (((a) % SCLUSTER_SIZE) * CLUSTER_N))
#else
#define CJ_BASE_INDEX(a, b) ((a)*CLUSTER_N * (b))
#endif
#else
#if CLUSTER_M == CLUSTER_N
#define CJ0_FROM_CI(a)      (a)
#define CJ1_FROM_CI(a)      (a)
#define CI_BASE_INDEX(a, b) ((a)*CLUSTER_N * (b))
#define CJ_BASE_INDEX(a, b) ((a)*CLUSTER_N * (b))
#elif CLUSTER_M == CLUSTER_N * 2 // M > N
#define CJ0_FROM_CI(a)      ((a) << 1)
#define CJ1_FROM_CI(a)      (((a) << 1) | 0x1)
#define CI_BASE_INDEX(a, b) ((a)*CLUSTER_M * (b))
#define CJ_BASE_INDEX(a, b) (((a) >> 1) * CLUSTER_M * (b) + ((a)&0x1) * (CLUSTER_M >> 1))
#elif CLUSTER_M == CLUSTER_N / 2 // M < N
#define CJ0_FROM_CI(a)      ((a) >> 1)
#define CJ1_FROM_CI(a)      ((a) >> 1)
#define CI_BASE_INDEX(a, b) (((a) >> 1) * CLUSTER_N * (b) + ((a)&0x1) * (CLUSTER_N >> 1))
#define CJ_BASE_INDEX(a, b) ((a)*CLUSTER_N * (b))
#else
#error "Invalid cluster configuration!"
#endif
#endif

#define CI_SCALAR_BASE_INDEX(a)  (CI_BASE_INDEX(a, 1))
#define CJ_SCALAR_BASE_INDEX(a)  (CJ_BASE_INDEX(a, 1))
#define CI_VECTOR3_BASE_INDEX(a)  (CI_BASE_INDEX(a, 3))
#define CJ_VECTOR3_BASE_INDEX(a)  (CJ_BASE_INDEX(a, 3))
#define CI_VECTOR_BASE_INDEX(a)  (CI_BASE_INDEX(a, ATOM_DIM))
#define CJ_VECTOR_BASE_INDEX(a)  (CJ_BASE_INDEX(a, ATOM_DIM))


#ifndef CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
// CPU
#ifndef CUDA_TARGET
#if CLUSTER_M >= CLUSTER_N
#define CL_X_OFFSET (0 * CLUSTER_M)
#define CL_Y_OFFSET (1 * CLUSTER_M)
#define CL_Z_OFFSET (2 * CLUSTER_M)
#define CL_X_INDEX_3D(i) ((i) + CL_X_OFFSET)
#define CL_Y_INDEX_3D(i) ((i) + CL_Y_OFFSET)
#define CL_Z_INDEX_3D(i) ((i) + CL_Z_OFFSET)
#else
#define CL_X_OFFSET (0 * CLUSTER_N)
#define CL_Y_OFFSET (1 * CLUSTER_N)
#define CL_Z_OFFSET (2 * CLUSTER_N)
#define CL_X_INDEX_3D(i) ((i) + CL_X_OFFSET)
#define CL_Y_INDEX_3D(i) ((i) + CL_Y_OFFSET)
#define CL_Z_INDEX_3D(i) ((i) + CL_Z_OFFSET)
#endif
#define CL_X_INDEX(i) CL_X_INDEX_3D(i)
#define CL_Y_INDEX(i) CL_Y_INDEX_3D(i)
#define CL_Z_INDEX(i) CL_Z_INDEX_3D(i)
#else
//gpusimple
#if CLUSTER_M == CLUSTER_N
#ifdef SOA_SUP
#define CL_X_OFFSET (0 * CLUSTER_N)
#define CL_Y_OFFSET (1 * CLUSTER_N)
#define CL_Z_OFFSET (2 * CLUSTER_N)
#define CL_X_INDEX_3D(i) ((i) + CL_X_OFFSET)
#define CL_Y_INDEX_3D(i) ((i) + CL_Y_OFFSET)
#define CL_Z_INDEX_3D(i) ((i) + CL_Z_OFFSET)
#define CL_X_INDEX(i) CL_X_INDEX_3D(i)
#define CL_Y_INDEX(i) CL_Y_INDEX_3D(i)
#define CL_Z_INDEX(i) CL_Z_INDEX_3D(i)
#else
#define CL_X_INDEX_3D(i)  ((i) * 3 + 0)
#define CL_Y_INDEX_3D(i)  ((i) * 3 + 1)
#define CL_Z_INDEX_3D(i)  ((i) * 3 + 2)
#define CL_X_INDEX(i)  CL_X_INDEX_3D(i)
#define CL_Y_INDEX(i)  CL_Y_INDEX_3D(i)
#define CL_Z_INDEX(i)  CL_Z_INDEX_3D(i)
#endif
#else
#error "Invalid cluster configuration for CUDA!"
#endif
#endif
#else
// supercluster
#ifdef SOA_SUP
#define CL_X_OFFSET (0 * CLUSTER_N * SCLUSTER_SIZE)
#define CL_Y_OFFSET (1 * CLUSTER_N * SCLUSTER_SIZE)
#define CL_Z_OFFSET (2 * CLUSTER_N * SCLUSTER_SIZE)
#define CL_X_INDEX_3D(i) ((i) + CL_X_OFFSET)
#define CL_Y_INDEX_3D(i) ((i) + CL_Y_OFFSET)
#define CL_Z_INDEX_3D(i) ((i) + CL_Z_OFFSET)
#define CL_X_INDEX(i) CL_X_INDEX_3D(i)
#define CL_Y_INDEX(i) CL_Y_INDEX_3D(i)
#define CL_Z_INDEX(i) CL_Z_INDEX_3D(i)
#else
#define CL_X_INDEX_3D(i)  ((i) * 3 + 0)
#define CL_Y_INDEX_3D(i)  ((i) * 3 + 1)
#define CL_Z_INDEX_3D(i)  ((i) * 3 + 2)
#define CL_X_INDEX(i)  ((i) * ATOM_DIM + 0)
#define CL_Y_INDEX(i)  ((i) * ATOM_DIM + 1)
#define CL_Z_INDEX(i)  ((i) * ATOM_DIM + 2)
#endif
#endif

#endif // __FORCE_H_
