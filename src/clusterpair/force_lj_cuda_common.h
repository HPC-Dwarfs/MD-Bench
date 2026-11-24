/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

#ifndef __FORCE_LJ_CUDA_COMMON_H_
#define __FORCE_LJ_CUDA_COMMON_H_

#include <pair_cuda.h>

#ifdef __HIP_PLATFORM_AMD__
#define FULL_WARP_MASK   0xffffffffffffffffULL
#define DPP_I32(v, ctrl) __builtin_amdgcn_update_dpp(0, v, ctrl, 0xf, 0xf, false)
#if PRECISION == 1
#define ROW_DPP(x, ctrl) __int_as_float(DPP_I32(__float_as_int(x), ctrl))
#else
#define ROW_DPP(x, ctrl)                                                                 \
    __hiloint2double(DPP_I32(__double2hiint(x), ctrl), DPP_I32(__double2loint(x), ctrl))
#endif
// AMD: no mask check of the *_sync shuffles, DPP row shift for _ROW
#define SHFL_DOWN(x, n)     __shfl_down(x, n)
#define SHFL_UP(x, n)       __shfl_up(x, n)
#define SHFL_DOWN_ROW(x, n) ROW_DPP(x, 0x100 + (n)) // row_shl:n
#define SHFL_UP_ROW(x, n)   ROW_DPP(x, 0x110 + (n)) // row_shr:n
#else
#define FULL_WARP_MASK      0xffffffffU
#define SHFL_DOWN(x, n)     __shfl_down_sync(FULL_WARP_MASK, x, n)
#define SHFL_UP(x, n)       __shfl_up_sync(FULL_WARP_MASK, x, n)
#define SHFL_DOWN_ROW(x, n) SHFL_DOWN(x, n)
#define SHFL_UP_ROW(x, n)   SHFL_UP(x, n)
#endif

#ifndef GPU_J_SLICES
#error "GPU_J_SLICES is not defined"
#endif

// Per-atom parameters
typedef struct __align__(4 * sizeof(MD_FLOAT)) {
    MD_FLOAT sqrt_epsilon; // geometric rule
    MD_FLOAT sigma3;       // geometric rule
    int mol;               // exclusions
    int type;              // full rule
} AtomParam;


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

extern int* cuda_cl_t;
extern MD_FLOAT* cuda_cl_sqrt_epsilon;
extern MD_FLOAT* cuda_cl_sigma3;
extern MD_FLOAT* cuda_cl_q;
extern int* cuda_cl_mol;
extern MD_FLOAT* cuda_cutforcesq;
extern MD_FLOAT* cuda_sigma6;
extern MD_FLOAT* cuda_epsilon;
// Packed data 
extern MD_FLOAT4* cuda_cl_xq;
extern AtomParam* cuda_cl_par;
extern MD_FLOAT4* cuda_nbfp;
}

// Arguments of the optimized force kernels
typedef struct {
    // Atoms
    const MD_FLOAT4* xq;
    const AtomParam* par;
    MD_FLOAT* f;
    int nclusters_local;
    // Neighbor list
    const int* numneigh;
    const int* cj;       // gpusimple
    const int2* cj_mask; // super-cluster: {cj, mask}
    int maxneighs;
    // LJ parameters
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq, sigma6, epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq;
#else
    const MD_FLOAT4* nbfp;
    int ntypes;
#endif
    CoulombScalar cs;
} NbForceArgs;

NbForceArgs makeNbForceArgs(Parameter* param, Atom* atom, Neighbor* neighbor);

typedef struct {
    MD_FLOAT cutforcesq; 
    MD_FLOAT sigma6; 
    MD_FLOAT epsilon;
} LjPair;

static __device__ __forceinline__ LjPair ljPair(
    const NbForceArgs& nb, const AtomParam& pi, const AtomParam& pj)
{
    LjPair p;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    p.cutforcesq = nb.cutforcesq;
    p.sigma6     = nb.sigma6;
    p.epsilon    = nb.epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    p.cutforcesq = nb.cutforcesq;
    p.sigma6     = pi.sigma3 * pj.sigma3;
    p.epsilon    = pi.sqrt_epsilon * pj.sqrt_epsilon;
#else
    const MD_FLOAT4 t = nb.nbfp[pi.type * nb.ntypes + pj.type];
    p.cutforcesq      = t.x;
    p.sigma6          = t.y;
    p.epsilon         = t.z;
#endif
    return p;
}

// Adds fj summed over the 8 lanes cii to atom j of f; lanes cii = 0, 1, 2 add x, y, z
static __device__ __forceinline__ void addForceJ(MD_FLOAT* f, int j, MD_FLOAT3 fj, int cii)
{
    fj.x += SHFL_DOWN_ROW(fj.x, 1);
    fj.y += SHFL_UP_ROW(fj.y, 1);
    fj.z += SHFL_DOWN_ROW(fj.z, 1);

    if (cii & 1) {
        fj.x = fj.y;
    }

    fj.x += SHFL_DOWN_ROW(fj.x, 2);
    fj.z += SHFL_UP_ROW(fj.z, 2);

    if (cii & 2) {
        fj.x = fj.z;
    }

    fj.x += SHFL_DOWN_ROW(fj.x, 4);

    if (cii < 3) {
        atomicAdd(&f[j * 3 + cii], fj.x);
    }
}

// Adds fi summed over the 8 lanes cjj to atom i of f
static __device__ __forceinline__ void addForceI(MD_FLOAT* f, int i, MD_FLOAT3 fi, int cjj)
{
    fi.x += SHFL_DOWN(fi.x, CLUSTER_M);
    fi.y += SHFL_UP(fi.y, CLUSTER_M);
    fi.z += SHFL_DOWN(fi.z, CLUSTER_M);

    if (cjj & 1) {
        fi.x = fi.y;
    }

    fi.x += SHFL_DOWN(fi.x, 2 * CLUSTER_M);
    fi.z += SHFL_UP(fi.z, 2 * CLUSTER_M);

    if (cjj & 2) {
        fi.x = fi.z;
    }

#ifdef __HIP_PLATFORM_AMD__
    // AMD wave64: add the cjj 4-7 half to cjj 0-3 within the wavefront
    if (warpSize == 64) {
        fi.x += SHFL_DOWN(fi.x, 4 * CLUSTER_M);
        if (cjj < 3) {
            atomicAdd(&f[i * 3 + cjj], fi.x);
        }
        return;
    }
#endif
    // One atomic per 32-lane half (NVIDIA, AMD wave32)
    if ((cjj & 3) < 3) {
        atomicAdd(&f[i * 3 + (cjj & 3)], fi.x);
    }
}

extern "C" void launchPackAtomsCuda(Atom* atom, int pack_par);
extern "C" int chooseJSlicesK(
    const Atom* atom, const Neighbor* neighbor, double k_slice, const void* kernel);

extern "C" void updatePbcCUDASimple(Atom* atom, Parameter* param);
extern "C" void updatePbcCUDASup(Atom* atom, Parameter* param);
extern "C" void reverseGhostForcesCUDASimple(Atom* atom, Parameter* param);
extern "C" void reverseGhostForcesCUDASup(Atom* atom, Parameter* param);
extern "C" void pruneNeighborCUDASimple(Parameter* param, Atom* atom, Neighbor* neighbor);
extern "C" void pruneNeighborCUDASup(Parameter* param, Atom* atom, Neighbor* neighbor);
extern "C" void launchPackNeighborsSup(Atom* atom, Neighbor* neighbor);
extern "C" void freeSupData(void);

#endif
