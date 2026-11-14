/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Device data of the clusterpair GPU force kernels, indexed by the atom slot
 * CJ_SCALAR_BASE_INDEX(cj) + cjj:
 *   cuda_cl_xq   {x, y, z, q}, packed every step
 *   cuda_cl_par  AtomParam, packed on list upload
 *   cuda_nbfp    {cutforcesq, sigma6, epsilon, 0} per type pair (full rule) */
#ifndef __FORCE_LJ_CUDA_COMMON_H_
#define __FORCE_LJ_CUDA_COMMON_H_

#include <pair_cuda.h>

#ifdef __HIP_PLATFORM_AMD__
#define FULL_WARP_MASK 0xffffffffffffffffULL
#else
#define FULL_WARP_MASK 0xffffffffU
#endif

// j-slices per i-cluster (grid.y): 0 = chosen per force call (auto), > 0 = fixed
#ifndef GPU_J_SLICES
#define GPU_J_SLICES 0
#endif

// Per-atom parameters (one 16-byte load), unused fields are 0
typedef struct __align__(4 * sizeof(MD_FLOAT)) {
    MD_FLOAT sqrt_epsilon; // geometric rule
    MD_FLOAT sigma3;       // geometric rule
    int mol;               // exclusions
    int type;              // full rule
} AtomParam;

// Device buffers, defined in force_lj_cuda_common.cu
extern "C" {
extern MD_FLOAT* cuda_cl_x;
extern MD_FLOAT* cuda_cl_v;
extern MD_FLOAT* cuda_cl_f;
extern int* cuda_natoms;
extern int* cuda_neighbors;
extern unsigned int* cuda_neighbors_imask;
extern int* cuda_numneigh;
extern int* cuda_numneigh_inner;
extern int* cuda_border_map;
extern int* cuda_jclusters_natoms;
extern int *cuda_PBCx, *cuda_PBCy, *cuda_PBCz;
// Unpacked per-atom data (NULL when unused)
extern int* cuda_cl_t;
extern MD_FLOAT* cuda_cl_sqrt_epsilon;
extern MD_FLOAT* cuda_cl_sigma3;
extern MD_FLOAT* cuda_cl_q;
extern int* cuda_cl_mol;
// Packed data of the force kernels
extern MD_FLOAT4* cuda_cl_xq;
extern AtomParam* cuda_cl_par;
extern MD_FLOAT4* cuda_nbfp;
}

// Arguments of the force kernels (passed by value)
typedef struct {
    const MD_FLOAT4* xq;
    const AtomParam* par;
    MD_FLOAT* f;
    int nclusters_local;
} NbAtoms;

typedef struct {
    const int* numneigh;
    const int* cj;       // gpusimple
    const int2* cj_mask; // super-cluster: {cj, mask}
    int maxneighs;
} NbList;

typedef struct {
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq, sigma6, epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq;
#else
    const MD_FLOAT4* nbfp;
    int ntypes;
#endif
} LjParams;

typedef struct {
    NbAtoms atoms;
    NbList nbl;
    LjParams lj;
    CoulombScalar cs;
} NbForceArgs;

// The super-cluster path sets nbl.cj_mask
NbForceArgs makeNbForceArgs(Parameter* param, Atom* atom, Neighbor* neighbor);

// LJ parameters of the pair i-j
typedef struct {
    MD_FLOAT cutforcesq, sigma6, epsilon;
} LjPair;

static __device__ __forceinline__ LjPair ljPair(
    const LjParams& lj, const AtomParam& pi, const AtomParam& pj)
{
    LjPair p;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    p.cutforcesq = lj.cutforcesq;
    p.sigma6     = lj.sigma6;
    p.epsilon    = lj.epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    p.cutforcesq = lj.cutforcesq;
    p.sigma6     = pi.sigma3 * pj.sigma3;
    p.epsilon    = pi.sqrt_epsilon * pj.sqrt_epsilon;
#else
    const MD_FLOAT4 t = lj.nbfp[pi.type * lj.ntypes + pj.type];
    p.cutforcesq      = t.x;
    p.sigma6          = t.y;
    p.epsilon         = t.z;
#endif
    return p;
}

// force_lj_cuda_common.cu
extern "C" void launchPackAtomsCuda(Atom* atom, int pack_par);
extern "C" int chooseJSlicesK(
    const Atom* atom, const Neighbor* neighbor, double k_slice, const void* kernel);

// force_lj_cuda.cu (gpusimple) and force_lj_cuda_sup.cu (super-cluster)
extern "C" void updatePbcCUDASimple(Atom* atom, Parameter* param);
extern "C" void updatePbcCUDASup(Atom* atom, Parameter* param);
extern "C" void reverseGhostForcesCUDASimple(Atom* atom, Parameter* param);
extern "C" void reverseGhostForcesCUDASup(Atom* atom, Parameter* param);
extern "C" void pruneNeighborCUDASimple(Parameter* param, Atom* atom, Neighbor* neighbor);
extern "C" void pruneNeighborCUDASup(Parameter* param, Atom* atom, Neighbor* neighbor);
extern "C" void launchPackNeighborsSup(Atom* atom, Neighbor* neighbor);
extern "C" void freeSupData(void);

#endif
