/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

#include <device.h>

#ifdef __HIP_PLATFORM_AMD__
#define FULL_WARP_MASK 0xffffffffffffffffULL
#else
#define FULL_WARP_MASK 0xffffffffU
#endif

extern "C" {
#include <stdio.h>
//---
#include <likwid-marker.h>
//---
#include <allocate.h>
#include <atom.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>
}

#include <nbdata_gpu.h>

extern __global__ void cudaInitialIntegrateSup_warp(MD_FLOAT* cuda_cl_x,
    MD_FLOAT* cuda_cl_v,
    MD_FLOAT* cuda_cl_f,
    int Nclusters_local,
    MD_FLOAT dtforce,
    MD_FLOAT dt);

extern __global__ void cudaFinalIntegrateSup_warp(
    MD_FLOAT* cuda_cl_v, MD_FLOAT* cuda_cl_f, int Nclusters_local, MD_FLOAT dtforce);

extern __global__ void cudaUpdatePbcSup_warp(MD_FLOAT* cuda_cl_x,
    int* cuda_border_map,
    int* cuda_jclusters_natoms,
    int* cuda_PBCx,
    int* cuda_PBCy,
    int* cuda_PBCz,
    int Nclusters_local,
    int Nclusters_ghost,
    MD_FLOAT param_xprd,
    MD_FLOAT param_yprd,
    MD_FLOAT param_zprd);

extern __global__ void cudaReverseGhostForcesSup_warp(MD_FLOAT* cuda_cl_f,
    int* cuda_border_map,
    int Nclusters_local,
    int Nclusters_ghost);

// Packed neighbor list of the super-cluster kernels (force_lj_cuda_sup.cu)
extern "C" void packSupData(Parameter* param, Atom* atom, Neighbor* neighbor);
extern "C" void freeSupData(void);

extern "C" {
MD_FLOAT* cuda_cl_x;
MD_FLOAT* cuda_cl_v;
MD_FLOAT* cuda_cl_f;
int* cuda_neighbors;
unsigned int* cuda_neighbors_imask;
int* cuda_numneigh;
int* cuda_numneigh_inner;
int* cuda_natoms;
int* natoms;
int* ngatoms;
int* cuda_border_map;
int* cuda_jclusters_natoms;
int *cuda_PBCx, *cuda_PBCy, *cuda_PBCz;

// Unpacked per-atom data (NULL when unused)
int* cuda_cl_t;
MD_FLOAT* cuda_cl_sqrt_epsilon;
MD_FLOAT* cuda_cl_sigma3;
MD_FLOAT* cuda_cl_q;
int* cuda_cl_mol;
// Packed data read by the force kernels, see nbdata_gpu.h
MD_FLOAT4* cuda_cl_xq;
AtomParam* cuda_cl_par;
MD_FLOAT4* cuda_nbfp;
}

// Packs atom slot a into cuda_cl_xq and, with pack_par, cuda_cl_par
__global__ void packAtomsCuda(const MD_FLOAT* cuda_cl_x,
    const MD_FLOAT* cuda_cl_q,
    const MD_FLOAT* cuda_cl_sqrt_epsilon,
    const MD_FLOAT* cuda_cl_sigma3,
    const int* cuda_cl_mol,
    const int* cuda_cl_t,
    MD_FLOAT4* xq,
    AtomParam* par,
    int nslots,
    int pack_par)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nslots) {
        return;
    }

    const MD_FLOAT* x = &cuda_cl_x[CJ_VECTOR_BASE_INDEX(a / CLUSTER_N)];
    const int i       = a % CLUSTER_N;
    MD_FLOAT4 v;
    v.x   = x[CL_X_INDEX(i)];
    v.y   = x[CL_Y_INDEX(i)];
    v.z   = x[CL_Z_INDEX(i)];
    v.w   = chargeCuda(cuda_cl_q, a);
    xq[a] = v;

    if (pack_par) {
        AtomParam p = {}; // unused fields stay 0
#if LJ_COMB_RULE == LJ_COMB_GEOM
        p.sqrt_epsilon = cuda_cl_sqrt_epsilon[a];
        p.sigma3       = cuda_cl_sigma3[a];
#elif LJ_COMB_RULE == LJ_COMB_FULL
        p.type = cuda_cl_t[a];
#endif
        p.mol  = moleculeCuda(cuda_cl_mol, a);
        par[a] = p;
    }
}

// Packs all local and ghost atom slots
extern "C" void launchPackAtomsCuda(Atom* atom, int pack_par)
{
    const int nslots = (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) *
                       CLUSTER_N;
    if (nslots == 0) {
        return;
    }
    packAtomsCuda<<<(nslots + 255) / 256, 256>>>(cuda_cl_x,
        cuda_cl_q,
        cuda_cl_sqrt_epsilon,
        cuda_cl_sigma3,
        cuda_cl_mol,
        cuda_cl_t,
        cuda_cl_xq,
        cuda_cl_par,
        nslots,
        pack_par);
    cuda_assert("packAtomsCuda", cudaPeekAtLastError());
}

// Arguments of the force kernels
NbForceArgs makeNbForceArgs(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    NbForceArgs a;
    a.atoms.xq              = cuda_cl_xq;
    a.atoms.par             = cuda_cl_par;
    a.atoms.f               = cuda_cl_f;
    a.atoms.nclusters_local = atom->Nclusters_local;
    a.nbl.numneigh          = cuda_numneigh_inner;
    a.nbl.cj                = cuda_neighbors;
    a.nbl.cj_mask           = NULL; // set by computeForceLJCudaSup()
    a.nbl.maxneighs         = neighbor->maxneighs;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    a.lj.cutforcesq = param->cutforce * param->cutforce;
    a.lj.sigma6     = param->sigma6;
    a.lj.epsilon    = param->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    a.lj.cutforcesq = param->cutforce * param->cutforce;
#else
    a.lj.nbfp   = cuda_nbfp;
    a.lj.ntypes = atom->ntypes;
#endif
    a.cs = coulombScalarInit(param);
    return a;
}

extern "C" void initDevice(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    cuda_cl_x = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * ATOM_DIM * sizeof(MD_FLOAT));
    cuda_cl_v = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
    cuda_cl_f = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
#if LJ_COMB_RULE == LJ_COMB_FULL
    cuda_cl_t = (int*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
    cuda_cl_sqrt_epsilon = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
    cuda_cl_sigma3 = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
#endif
#if LJ_COMB_RULE == LJ_COMB_FULL
    // {cutforcesq, sigma6, epsilon} per type pair
    const int ntypes2 = atom->ntypes * atom->ntypes;
    MD_FLOAT4* nbfp   = (MD_FLOAT4*)malloc(ntypes2 * sizeof(MD_FLOAT4));
    for (int t = 0; t < ntypes2; t++) {
        nbfp[t].x = atom->cutforcesq[t];
        nbfp[t].y = atom->sigma6[t];
        nbfp[t].z = atom->epsilon[t];
        nbfp[t].w = (MD_FLOAT)0.0;
    }
    cuda_nbfp = (MD_FLOAT4*)allocateGPU(ntypes2 * sizeof(MD_FLOAT4));
    memcpyToGPU(cuda_nbfp, nbfp, ntypes2 * sizeof(MD_FLOAT4));
    free(nbfp);
#endif
#if COULOMB_RULE != COULOMB_NONE
    cuda_cl_q = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
#endif
#if USE_EXCLUSIONS
    cuda_cl_mol = (int*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
#endif
    cuda_cl_xq = (MD_FLOAT4*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT4));
    cuda_cl_par = (AtomParam*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(AtomParam));
    cuda_natoms = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_jclusters_natoms = (int*)allocateGPU(
        atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_border_map = (int*)allocateGPU(
        atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_PBCx     = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_PBCy     = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_PBCz     = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_numneigh = (int*)allocateGPU(atom->Nclusters_max * sizeof(int));
    cuda_numneigh_inner = (int*)allocateGPU(atom->Nclusters_max * sizeof(int));
    cuda_neighbors      = (int*)allocateGPU(
        atom->Nclusters_max * neighbor->maxneighs * sizeof(int));
    cuda_neighbors_imask = (unsigned int*)allocateGPU(
        atom->Nclusters_max * neighbor->maxneighs * sizeof(unsigned int));
    // Pinned staging buffers, refilled on every list upload
    natoms = (int*)allocateHostPinned(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    ngatoms = (int*)allocateHostPinned(
        atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
}

extern "C" void copyDataToCUDADevice(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    memcpyToGPU(cuda_cl_x,
        atom->cl_x,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost + 1) * CLUSTER_M * ATOM_DIM * sizeof(MD_FLOAT));
    memcpyToGPU(cuda_cl_v,
        atom->cl_v,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));

#if LJ_COMB_RULE == LJ_COMB_FULL
    memcpyToGPU(cuda_cl_t,
        atom->cl_t,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
            sizeof(int));
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
    memcpyToGPU(cuda_cl_sqrt_epsilon,
        atom->cl_sqrt_epsilon,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
            sizeof(MD_FLOAT));
    memcpyToGPU(cuda_cl_sigma3,
        atom->cl_sigma3,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
            sizeof(MD_FLOAT));
#endif
#if COULOMB_RULE != COULOMB_NONE
    memcpyToGPU(cuda_cl_q,
        atom->cl_q,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
            sizeof(MD_FLOAT));
#endif
#if USE_EXCLUSIONS
    memcpyToGPU(cuda_cl_mol,
        atom->cl_mol,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
            sizeof(int));
#endif
    launchPackAtomsCuda(atom, 1);

    for (int ci = 0; ci < atom->Nclusters_local * SCLUSTER_SIZE; ci++) {
        natoms[ci] = atom->iclusters[ci].natoms;
    }

    // int ncj  = get_ncj_from_nci(atom->Nclusters_local);
    int ncj = CJ0_FROM_CI(atom->Nclusters_local);
    for (int cg = 0; cg < atom->Nclusters_ghost; cg++) {
        const int cj = ncj + cg;
        ngatoms[cg]  = atom->jclusters[cj].natoms;
    }

    memcpyToGPU(cuda_natoms, natoms, atom->Nclusters_local * SCLUSTER_SIZE * sizeof(int));
    memcpyToGPU(cuda_jclusters_natoms, ngatoms, atom->Nclusters_ghost * sizeof(int));
#ifndef _MPI
    memcpyToGPU(cuda_border_map, atom->border_map, atom->Nclusters_ghost * sizeof(int));
    memcpyToGPU(cuda_PBCx, atom->PBCx, atom->Nclusters_ghost * sizeof(int));
    memcpyToGPU(cuda_PBCy, atom->PBCy, atom->Nclusters_ghost * sizeof(int));
    memcpyToGPU(cuda_PBCz, atom->PBCz, atom->Nclusters_ghost * sizeof(int));
#endif

    memcpyToGPU(cuda_numneigh, neighbor->numneigh, atom->Nclusters_local * sizeof(int));
    memcpyToGPU(cuda_numneigh_inner,
        neighbor->numneigh_inner,
        atom->Nclusters_local * sizeof(int));
    memcpyToGPU(cuda_neighbors,
        neighbor->neighbors,
        atom->Nclusters_local * neighbor->maxneighs * sizeof(int));

    if (param->super_clustering) {
        // The host build computed the sub-cluster masks
        memcpyToGPU(cuda_neighbors_imask,
            neighbor->neighbors_imask,
            atom->Nclusters_local * neighbor->maxneighs * sizeof(unsigned int));
        packSupData(param, atom, neighbor);
    } else {
        // Unused by the non-super GPU force kernel.
        memsetGPU(cuda_neighbors_imask,
            0xff,
            atom->Nclusters_local * neighbor->maxneighs * sizeof(unsigned int));
    }
}

extern "C" void copyDataFromCUDADevice(Parameter* param, Atom* atom)
{
    memcpyFromGPU(atom->cl_x,
        cuda_cl_x,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost + 1) * CLUSTER_M * ATOM_DIM * sizeof(MD_FLOAT));
    memcpyFromGPU(atom->cl_v,
        cuda_cl_v,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));
}

extern "C" void cudaDeviceFree(Parameter* param)
{
    GPUfree(cuda_cl_x);
    GPUfree(cuda_cl_v);
    GPUfree(cuda_cl_f);
#if LJ_COMB_RULE == LJ_COMB_FULL
    GPUfree(cuda_cl_t);
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
    GPUfree(cuda_cl_sqrt_epsilon);
    GPUfree(cuda_cl_sigma3);
#endif
#if LJ_COMB_RULE == LJ_COMB_FULL
    GPUfree(cuda_nbfp);
#endif
#if COULOMB_RULE != COULOMB_NONE
    GPUfree(cuda_cl_q);
#endif
#if USE_EXCLUSIONS
    GPUfree(cuda_cl_mol);
#endif
    GPUfree(cuda_cl_xq);
    GPUfree(cuda_cl_par);
    GPUfree(cuda_numneigh);
    GPUfree(cuda_numneigh_inner);
    GPUfree(cuda_neighbors);
    GPUfree(cuda_neighbors_imask);
    GPUfree(cuda_natoms);
    GPUfree(cuda_border_map);
    GPUfree(cuda_jclusters_natoms);
    GPUfree(cuda_PBCx);
    GPUfree(cuda_PBCy);
    GPUfree(cuda_PBCz);

    freeHostPinned(natoms);
    freeHostPinned(ngatoms);
    freeSupData();
}

// gpusimple kernels: one block per i-cluster and j-slice, one pair per thread
__global__ void computeForceLJCudaFullNeigh(const NbForceArgs a)
{
    int ci = blockIdx.x;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    int cii = threadIdx.y;
    int cjj = threadIdx.x;
#else
    int cii = threadIdx.x;
    int cjj = threadIdx.y;
#endif
    int ci_cj0     = CJ0_FROM_CI(ci);
    MD_FLOAT* ci_f = &a.atoms.f[CI_VECTOR3_BASE_INDEX(ci)];
    int numneighs  = a.nbl.numneigh[ci];
    MD_FLOAT fix   = (MD_FLOAT)0.0;
    MD_FLOAT fiy   = (MD_FLOAT)0.0;
    MD_FLOAT fiz   = (MD_FLOAT)0.0;

    // i-atom (charge scaled by epsfac)
    const int slot_i   = CI_SCALAR_BASE_INDEX(ci) + cii;
    const MD_FLOAT4 xi = a.atoms.xq[slot_i];
    const AtomParam pi = a.atoms.par[slot_i];
    const MD_FLOAT qi  = a.cs.epsfac * xi.w;

    // One j-cluster of the neighbor list
    auto jcluster = [&](int k) {
        const int cj =
            neighs_gpu(a.nbl.cj, ci, k, a.atoms.nclusters_local, a.nbl.maxneighs);

        int cond;
#if CLUSTER_M == CLUSTER_N
        cond = ci_cj0 != cj || cii != cjj;
#elif CLUSTER_M < CLUSTER_N
        cond = ci_cj0 != cj || cii + CLUSTER_M * (ci & 0x1) != cjj;
#endif
        if (cond) {
            const int slot_j   = CJ_SCALAR_BASE_INDEX(cj) + cjj;
            const MD_FLOAT4 xj = a.atoms.xq[slot_j];
            const AtomParam pj = a.atoms.par[slot_j];
            MD_FLOAT delx      = xi.x - xj.x;
            MD_FLOAT dely      = xi.y - xj.y;
            MD_FLOAT delz      = xi.z - xj.z;
            MD_FLOAT rsq       = delx * delx + dely * dely + delz * delz;
            const LjPair p     = ljPair(a.lj, pi, pj);

            if (rsq < p.cutforcesq) {
                MD_FLOAT force = pairForceScalar(rsq,
                    p.sigma6,
                    p.epsilon,
                    qi,
                    xj.w,
                    pi.mol,
                    pj.mol,
                    a.cs);

                fix += delx * force;
                fiy += dely * force;
                fiz += delz * force;
            }
        }
    };

    // j-slices: block (ci, s) handles neighbors s, s+S, ... (unit stride without slices)
    if (gridDim.y == 1) {
        for (int k = 0; k < numneighs; k++) {
            jcluster(k);
        }
    } else {
#pragma unroll 2
        for (int k = blockIdx.y; k < numneighs; k += gridDim.y) {
            jcluster(k);
        }
    }

    // Reduce the i-forces with warp shuffles
#if CLUSTER_M <= 32
    auto mask = FULL_WARP_MASK;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    for (int offset = CLUSTER_M / 2; offset > 0; offset /= 2) {
#ifdef CUDA_TARGET
        fix += __shfl_down_sync(mask, fix, offset);
        fiy += __shfl_down_sync(mask, fiy, offset);
        fiz += __shfl_down_sync(mask, fiz, offset);
#else
        fix += __shfl_down(fix, offset);
        fiy += __shfl_down(fiy, offset);
        fiz += __shfl_down(fiz, offset);
#endif
    }

    if (cjj == 0) {
        if (gridDim.y == 1) {
            ci_f[CL_X_INDEX_3D(cii)] = fix;
            ci_f[CL_Y_INDEX_3D(cii)] = fiy;
            ci_f[CL_Z_INDEX_3D(cii)] = fiz;
        } else {
            // several j-slices add to the same i-cluster
            atomicAdd(&ci_f[CL_X_INDEX_3D(cii)], fix);
            atomicAdd(&ci_f[CL_Y_INDEX_3D(cii)], fiy);
            atomicAdd(&ci_f[CL_Z_INDEX_3D(cii)], fiz);
        }
    }
#else
    fix += __shfl_down_sync(mask, fix, CLUSTER_M);
    fiy += __shfl_up_sync(mask, fiy, CLUSTER_M);
    fiz += __shfl_down_sync(mask, fiz, CLUSTER_M);

    if (cjj & 1) {
        fix = fiy;
    }

    fix += __shfl_down_sync(mask, fix, 2 * CLUSTER_M);
    fiz += __shfl_up_sync(mask, fiz, 2 * CLUSTER_M);
    if (cjj & 2) {
        fix = fiz;
    }

    /* Threads 0,1,2 and 4,5,6 increment x,y,z for their warp */
    if ((cjj & 3) < 3) {
        atomicAdd(&ci_f[cii * 3 + (cjj & 3)], fix);
    }

#endif
#else
    atomicAdd(&ci_f[CL_X_INDEX_3D(cii)], fix);
    atomicAdd(&ci_f[CL_Y_INDEX_3D(cii)], fiy);
    atomicAdd(&ci_f[CL_Z_INDEX_3D(cii)], fiz);
#endif
}

__global__ void computeForceLJCudaHalfNeigh(const NbForceArgs a)
{
    int ci = blockIdx.x;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    int cii = threadIdx.y;
    int cjj = threadIdx.x;
#else
    int cii = threadIdx.x;
    int cjj = threadIdx.y;
    MD_FLOAT3 fcj_buf; // j-force of this thread's j-atom from one j-cluster
#endif
    int ci_cj0     = CJ0_FROM_CI(ci);
    MD_FLOAT* ci_f = &a.atoms.f[CI_VECTOR3_BASE_INDEX(ci)];
    int numneighs  = a.nbl.numneigh[ci];
    MD_FLOAT fix   = 0;
    MD_FLOAT fiy   = 0;
    MD_FLOAT fiz   = 0;

    // i-atom (charge scaled by epsfac)
    const int slot_i   = CI_SCALAR_BASE_INDEX(ci) + cii;
    const MD_FLOAT4 xi = a.atoms.xq[slot_i];
    const AtomParam pi = a.atoms.par[slot_i];
    const MD_FLOAT qi  = a.cs.epsfac * xi.w;

    // j-slices: block (ci, s) handles neighbors s, s+S, ... (unrolled twice)
#pragma unroll 2
    for (int k = blockIdx.y; k < numneighs; k += gridDim.y) {
        const int cj =
            neighs_gpu(a.nbl.cj, ci, k, a.atoms.nclusters_local, a.nbl.maxneighs);
        MD_FLOAT* cj_f = &a.atoms.f[CJ_VECTOR3_BASE_INDEX(cj)];
#ifndef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        fcj_buf.x = fcj_buf.y = fcj_buf.z = (MD_FLOAT)0.0;
#endif

        int cond;
#if CLUSTER_M == CLUSTER_N
        cond = ci_cj0 != cj || cii < cjj;
#elif CLUSTER_M < CLUSTER_N
        cond = ci_cj0 != cj || cii + CLUSTER_M * (ci & 0x1) < cjj;
#endif
        if (cond) {
            const int slot_j   = CJ_SCALAR_BASE_INDEX(cj) + cjj;
            const MD_FLOAT4 xj = a.atoms.xq[slot_j];
            const AtomParam pj = a.atoms.par[slot_j];
            MD_FLOAT delx      = xi.x - xj.x;
            MD_FLOAT dely      = xi.y - xj.y;
            MD_FLOAT delz      = xi.z - xj.z;
            MD_FLOAT rsq       = delx * delx + dely * dely + delz * delz;
            const LjPair p     = ljPair(a.lj, pi, pj);

            if (rsq < p.cutforcesq) {
                MD_FLOAT force = pairForceScalar(rsq,
                    p.sigma6,
                    p.epsilon,
                    qi,
                    xj.w,
                    pi.mol,
                    pj.mol,
                    a.cs);
                MD_FLOAT partial_force_x = delx * force;
                MD_FLOAT partial_force_y = dely * force;
                MD_FLOAT partial_force_z = delz * force;

#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
                atomicAdd(&cj_f[CL_X_INDEX_3D(cjj)], -partial_force_x);
                atomicAdd(&cj_f[CL_Y_INDEX_3D(cjj)], -partial_force_y);
                atomicAdd(&cj_f[CL_Z_INDEX_3D(cjj)], -partial_force_z);
#else
                fcj_buf.x -= partial_force_x;
                fcj_buf.y -= partial_force_y;
                fcj_buf.z -= partial_force_z;
#endif

                fix += partial_force_x;
                fiy += partial_force_y;
                fiz += partial_force_z;
            }
        }

#ifndef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        // Reduce the j-force over the 8 i-lanes, one atomic per component
        auto jmask = FULL_WARP_MASK;

        fcj_buf.x += __shfl_down_sync(jmask, fcj_buf.x, 1);
        fcj_buf.y += __shfl_up_sync(jmask, fcj_buf.y, 1);
        fcj_buf.z += __shfl_down_sync(jmask, fcj_buf.z, 1);

        if (cii & 1) {
            fcj_buf.x = fcj_buf.y;
        }

        fcj_buf.x += __shfl_down_sync(jmask, fcj_buf.x, 2);
        fcj_buf.z += __shfl_up_sync(jmask, fcj_buf.z, 2);

        if (cii & 2) {
            fcj_buf.x = fcj_buf.z;
        }

        fcj_buf.x += __shfl_down_sync(jmask, fcj_buf.x, 4);

        // Lanes cii = 0, 1, 2 now hold the x, y, z component of the j-force
        if (cii < 3) {
            const int comp = (cii == 0)   ? CL_X_INDEX_3D(cjj)
                             : (cii == 1) ? CL_Y_INDEX_3D(cjj)
                                          : CL_Z_INDEX_3D(cjj);
            atomicAdd(&cj_f[comp], fcj_buf.x);
        }
#endif
    }

    // Reduce the i-forces with warp shuffles
#if CLUSTER_M <= 32
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    atomicAdd(&ci_f[CL_X_INDEX_3D(cii)], fix);
    atomicAdd(&ci_f[CL_Y_INDEX_3D(cii)], fiy);
    atomicAdd(&ci_f[CL_Z_INDEX_3D(cii)], fiz);
#else
    auto mask = FULL_WARP_MASK;
    fix += __shfl_down_sync(mask, fix, CLUSTER_M);
    fiy += __shfl_up_sync(mask, fiy, CLUSTER_M);
    fiz += __shfl_down_sync(mask, fiz, CLUSTER_M);

    if (cjj & 1) {
        fix = fiy;
    }

    fix += __shfl_down_sync(mask, fix, 2 * CLUSTER_M);
    fiz += __shfl_up_sync(mask, fiz, 2 * CLUSTER_M);
    if (cjj & 2) {
        fix = fiz;
    }

    /* Threads 0,1,2 and 4,5,6 increment x,y,z for their warp */
    if ((cjj & 3) < 3) {
        atomicAdd(&ci_f[cii * 3 + (cjj & 3)], fix);
    }

#endif
#else
    atomicAdd(&ci_f[CL_X_INDEX_3D(cii)], fix);
    atomicAdd(&ci_f[CL_Y_INDEX_3D(cii)], fiy);
    atomicAdd(&ci_f[CL_Z_INDEX_3D(cii)], fiz);
#endif
}

__global__ void cudaInitialIntegrate_warp(MD_FLOAT* cuda_cl_x,
    MD_FLOAT* cuda_cl_v,
    MD_FLOAT* cuda_cl_f,
    int* cuda_natoms,
    int Nclusters_local,
    MD_FLOAT dtforce,
    MD_FLOAT dt)
{

    unsigned int ci = blockDim.x * blockIdx.x + threadIdx.x;
    if (ci >= Nclusters_local) {
        return;
    }

    MD_FLOAT* ci_x  = &cuda_cl_x[CI_VECTOR_BASE_INDEX(ci)];
    int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
    MD_FLOAT* ci_v  = &cuda_cl_v[ci_vec_base];
    MD_FLOAT* ci_f  = &cuda_cl_f[ci_vec_base];

    for (int cii = 0; cii < cuda_natoms[ci]; cii++) {
        ci_v[CL_X_INDEX_3D(cii)] += dtforce * ci_f[CL_X_INDEX_3D(cii)];
        ci_v[CL_Y_INDEX_3D(cii)] += dtforce * ci_f[CL_Y_INDEX_3D(cii)];
        ci_v[CL_Z_INDEX_3D(cii)] += dtforce * ci_f[CL_Z_INDEX_3D(cii)];
        ci_x[CL_X_INDEX(cii)] += dt * ci_v[CL_X_INDEX_3D(cii)];
        ci_x[CL_Y_INDEX(cii)] += dt * ci_v[CL_Y_INDEX_3D(cii)];
        ci_x[CL_Z_INDEX(cii)] += dt * ci_v[CL_Z_INDEX_3D(cii)];
    }
}

extern "C" void cudaInitialIntegrate(Parameter* param, Atom* atom)
{
    const int threads_num = 64;
    const int N           = atom->Nclusters_local;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

    cudaInitialIntegrate_warp<<<grid_size, block_size>>>(cuda_cl_x,
        cuda_cl_v,
        cuda_cl_f,
        cuda_natoms,
        atom->Nclusters_local,
        param->dtforce,
        param->dt);

    cuda_assert("cudaInitialIntegrate", cudaPeekAtLastError());
    cuda_assert("cudaInitialIntegrate", cudaDeviceSynchronize());
}

__global__ void cudaFinalIntegrate_warp(MD_FLOAT* cuda_cl_v,
    MD_FLOAT* cuda_cl_f,
    int* cuda_natoms,
    int Nclusters_local,
    MD_FLOAT dtforce)
{

    unsigned int ci = blockDim.x * blockIdx.x + threadIdx.x;
    if (ci >= Nclusters_local) {
        return;
    }

    int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
    MD_FLOAT* ci_v  = &cuda_cl_v[ci_vec_base];
    MD_FLOAT* ci_f  = &cuda_cl_f[ci_vec_base];

    for (int cii = 0; cii < cuda_natoms[ci]; cii++) {
        ci_v[CL_X_INDEX_3D(cii)] += dtforce * ci_f[CL_X_INDEX_3D(cii)];
        ci_v[CL_Y_INDEX_3D(cii)] += dtforce * ci_f[CL_Y_INDEX_3D(cii)];
        ci_v[CL_Z_INDEX_3D(cii)] += dtforce * ci_f[CL_Z_INDEX_3D(cii)];
    }
}

extern "C" void cudaFinalIntegrate(Parameter* param, Atom* atom)
{
    const int threads_num = 64;
    const int N           = atom->Nclusters_local;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

    cudaFinalIntegrate_warp<<<grid_size, block_size>>>(cuda_cl_v,
        cuda_cl_f,
        cuda_natoms,
        atom->Nclusters_local,
        param->dtforce);

    cuda_assert("cudaFinalIntegrate", cudaPeekAtLastError());
    cuda_assert("cudaFinalIntegrate", cudaDeviceSynchronize());
}

__global__ void cudaUpdatePbc_warp(MD_FLOAT* cuda_cl_x,
    int* cuda_border_map,
    int* cuda_jclusters_natoms,
    int* cuda_PBCx,
    int* cuda_PBCy,
    int* cuda_PBCz,
    int Nclusters_local,
    int Nclusters_ghost,
    MD_FLOAT param_xprd,
    MD_FLOAT param_yprd,
    MD_FLOAT param_zprd)
{

    unsigned int cg = blockDim.x * blockIdx.x + threadIdx.x;
    if (cg >= Nclusters_ghost) {
        return;
    }

    // int ncj       = get_ncj_from_nci(Nclusters_local);
    int ncj       = CJ0_FROM_CI(Nclusters_local);
    MD_FLOAT xprd = param_xprd;
    MD_FLOAT yprd = param_yprd;
    MD_FLOAT zprd = param_zprd;

    const int cj      = ncj + cg;
    int cj_vec_base   = CJ_VECTOR_BASE_INDEX(cj);
    int bmap_vec_base = CJ_VECTOR_BASE_INDEX(cuda_border_map[cg]);
    MD_FLOAT* cj_x    = &cuda_cl_x[cj_vec_base];
    MD_FLOAT* bmap_x  = &cuda_cl_x[bmap_vec_base];

    for (int cjj = 0; cjj < cuda_jclusters_natoms[cg]; cjj++) {
        cj_x[CL_X_INDEX(cjj)] = bmap_x[CL_X_INDEX(cjj)] + cuda_PBCx[cg] * xprd;
        cj_x[CL_Y_INDEX(cjj)] = bmap_x[CL_Y_INDEX(cjj)] + cuda_PBCy[cg] * yprd;
        cj_x[CL_Z_INDEX(cjj)] = bmap_x[CL_Z_INDEX(cjj)] + cuda_PBCz[cg] * zprd;
    }
}

/* update coordinates of ghost atoms */
/* uses mapping created in setupPbc */
extern "C" void updatePbcCUDA(Atom* atom, Parameter* param)
{
    const int N           = atom->Nclusters_ghost;
    const int threads_num = 64;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

    if (param->super_clustering) {
        cudaUpdatePbcSup_warp<<<grid_size, block_size>>>(cuda_cl_x,
            cuda_border_map,
            cuda_jclusters_natoms,
            cuda_PBCx,
            cuda_PBCy,
            cuda_PBCz,
            atom->Nclusters_local,
            atom->Nclusters_ghost,
            param->xprd,
            param->yprd,
            param->zprd);
    } else {
        cudaUpdatePbc_warp<<<grid_size, block_size>>>(cuda_cl_x,
            cuda_border_map,
            cuda_jclusters_natoms,
            cuda_PBCx,
            cuda_PBCy,
            cuda_PBCz,
            atom->Nclusters_local,
            atom->Nclusters_ghost,
            param->xprd,
            param->yprd,
            param->zprd);
    }

    cuda_assert("cudaUpdatePbc", cudaPeekAtLastError());
    cuda_assert("cudaUpdatePbc", cudaDeviceSynchronize());
}

// Adds the ghost forces to their real clusters (atomic: several ghosts per cluster)
__global__ void cudaReverseGhostForces_warp(MD_FLOAT* cuda_cl_f,
    int* cuda_border_map,
    int* cuda_jclusters_natoms,
    int Nclusters_local,
    int Nclusters_ghost)
{
    unsigned int cg = blockDim.x * blockIdx.x + threadIdx.x;
    if (cg >= Nclusters_ghost) {
        return;
    }

    int ncj           = CJ0_FROM_CI(Nclusters_local);
    const int cj      = ncj + cg;
    int cj_vec_base   = CJ_VECTOR3_BASE_INDEX(cj);
    int bmap_vec_base = CJ_VECTOR3_BASE_INDEX(cuda_border_map[cg]);
    MD_FLOAT* cj_f     = &cuda_cl_f[cj_vec_base];
    MD_FLOAT* bmap_f   = &cuda_cl_f[bmap_vec_base];

    for (int cjj = 0; cjj < cuda_jclusters_natoms[cg]; cjj++) {
        atomicAdd(&bmap_f[CL_X_INDEX_3D(cjj)], cj_f[CL_X_INDEX_3D(cjj)]);
        atomicAdd(&bmap_f[CL_Y_INDEX_3D(cjj)], cj_f[CL_Y_INDEX_3D(cjj)]);
        atomicAdd(&bmap_f[CL_Z_INDEX_3D(cjj)], cj_f[CL_Z_INDEX_3D(cjj)]);
    }
}

extern "C" void reverseGhostForcesCUDA(Atom* atom, Parameter* param)
{
    const int N           = atom->Nclusters_ghost;
    const int threads_num = 64;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

    if (param->super_clustering) {
        cudaReverseGhostForcesSup_warp<<<grid_size, block_size>>>(cuda_cl_f,
            cuda_border_map,
            atom->Nclusters_local,
            atom->Nclusters_ghost);
    } else {
        cudaReverseGhostForces_warp<<<grid_size, block_size>>>(cuda_cl_f,
            cuda_border_map,
            cuda_jclusters_natoms,
            atom->Nclusters_local,
            atom->Nclusters_ghost);
    }

    cuda_assert("cudaReverseGhostForces", cudaPeekAtLastError());
    cuda_assert("cudaReverseGhostForces", cudaDeviceSynchronize());
}

__global__ void cudaPruneNeighbor(MD_FLOAT* cuda_cl_x,
    int* cuda_numneigh,
    int* cuda_numneigh_inner,
    int* cuda_neighbors,
    int Nclusters_local,
    int maxneighs,
    MD_FLOAT cutsq)
{
    unsigned int ci = blockDim.x * blockIdx.x + threadIdx.x;
    if (ci >= Nclusters_local) {
        return;
    }

    const int numneighs = cuda_numneigh[ci];
    MD_FLOAT* ci_x      = &cuda_cl_x[CI_VECTOR_BASE_INDEX(ci)];
    int lo              = 0;

    for (int hi = 0; hi < numneighs; hi++) {
        int cj         = neighs_gpu(cuda_neighbors, ci, hi, Nclusters_local, maxneighs);
        MD_FLOAT* cj_x = &cuda_cl_x[CJ_VECTOR_BASE_INDEX(cj)];
        int is_inner   = 0;

        for (int cii = 0; cii < CLUSTER_M && !is_inner; cii++) {
            MD_FLOAT xtmp = ci_x[CL_X_INDEX(cii)];
            MD_FLOAT ytmp = ci_x[CL_Y_INDEX(cii)];
            MD_FLOAT ztmp = ci_x[CL_Z_INDEX(cii)];
            for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
                MD_FLOAT delx = xtmp - cj_x[CL_X_INDEX(cjj)];
                MD_FLOAT dely = ytmp - cj_x[CL_Y_INDEX(cjj)];
                MD_FLOAT delz = ztmp - cj_x[CL_Z_INDEX(cjj)];
                if (delx * delx + dely * dely + delz * delz < cutsq) {
                    is_inner = 1;
                    break;
                }
            }
        }

        if (is_inner) {
            if (hi != lo) {
                int t_cj = neighs_gpu(cuda_neighbors, ci, lo, Nclusters_local, maxneighs);
                neighs_gpu(cuda_neighbors, ci, lo, Nclusters_local, maxneighs) = cj;
                neighs_gpu(cuda_neighbors, ci, hi, Nclusters_local, maxneighs) = t_cj;
            }
            lo++;
        }
    }

    cuda_numneigh_inner[ci] = lo;
}

extern "C" void pruneNeighborCUDASup(Parameter*, Atom*, Neighbor*);

extern "C" void pruneNeighborCUDA(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("pruneNeighborCUDA start\n");

    if (param->super_clustering) {
        pruneNeighborCUDASup(param, atom, neighbor);
        return;
    }

    if (param->outer_skin <= 0.0) {
        // No outer list: use the full counts
        memcpyOnGPU(cuda_numneigh_inner,
            cuda_numneigh,
            atom->Nclusters_local * sizeof(int));
        return;
    }

    const MD_FLOAT cut_inner = param->cutforce + param->skin;
    const MD_FLOAT cutsq     = cut_inner * cut_inner;
    const int threads_num    = 64;
    const int N              = atom->Nclusters_local;
    dim3 block_size          = dim3(threads_num, 1, 1);
    dim3 grid_size           = dim3((N + threads_num - 1) / threads_num, 1, 1);

    cudaPruneNeighbor<<<grid_size, block_size>>>(cuda_cl_x,
        cuda_numneigh,
        cuda_numneigh_inner,
        cuda_neighbors,
        atom->Nclusters_local,
        neighbor->maxneighs,
        cutsq);

    cuda_assert("cudaPruneNeighbor", cudaPeekAtLastError());
    cuda_assert("cudaPruneNeighbor", cudaDeviceSynchronize());
    DEBUG_MESSAGE("pruneNeighborCUDA end\n");
}

// j-slices, see force_lj_cuda_sup.cu
extern "C" int chooseJSlicesK(
    const Atom* atom, const Neighbor* neighbor, double k_slice, const void* kernel);

// K of the gpusimple kernels (tuned on an RTX 2080)
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
#ifndef K_SLICE_SIMPLE_HALF
#define K_SLICE_SIMPLE_HALF 2.0
#endif
#ifndef K_SLICE_SIMPLE_FULL
#define K_SLICE_SIMPLE_FULL 2.0
#endif
#else
#ifndef K_SLICE_SIMPLE_HALF
#define K_SLICE_SIMPLE_HALF 3.3
#endif
#ifndef K_SLICE_SIMPLE_FULL
#define K_SLICE_SIMPLE_FULL 2.5
#endif
#endif

extern "C" double computeForceLJCuda(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCuda start\r\n");

    const NbForceArgs args = makeNbForceArgs(param, atom, neighbor);

    // memsetGPU(cuda_cl_f, 0, atom->Nclusters_local * CLUSTER_M * 3 * sizeof(MD_FLOAT));
    memsetGPU(cuda_cl_f,
        0,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));
    dim3 block_size = dim3(CLUSTER_N, CLUSTER_M, 1);
    dim3 grid_size  = dim3(atom->Nclusters_local, 1, 1);
    double S        = getTimeStamp();
    LIKWID_MARKER_START("force");

    // Repack x/y/z/q (the positions moved), inside the timed region
    launchPackAtomsCuda(atom, 0);

    if (neighbor->half_neigh) {
        grid_size.y = chooseJSlicesK(atom,
            neighbor,
            K_SLICE_SIMPLE_HALF,
            (const void*)computeForceLJCudaHalfNeigh);
        computeForceLJCudaHalfNeigh<<<grid_size, block_size>>>(args);
    } else {
        grid_size.y = chooseJSlicesK(atom,
            neighbor,
            K_SLICE_SIMPLE_FULL,
            (const void*)computeForceLJCudaFullNeigh);
        computeForceLJCudaFullNeigh<<<grid_size, block_size>>>(args);
    }

    cuda_assert("computeForceLJ_cuda", cudaPeekAtLastError());
    cuda_assert("computeForceLJ_cuda", cudaDeviceSynchronize());
    LIKWID_MARKER_STOP("force");

#ifdef COMPUTE_STATS
    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        addStat(stats->calculated_forces, 1);
        addStat(stats->num_neighs, neighbor->numneigh[ci]);
    }
#endif

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJCuda stop\r\n");
    return E - S;
}

extern "C" void copyGhostFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->cl_x,
        cuda_cl_x,
        atom->Nclusters_local * CLUSTER_M * SCLUSTER_SIZE * ATOM_DIM * sizeof(MD_FLOAT));
}

extern "C" void copyGhostToGPU(Atom* atom)
{
    // The ghost clusters follow the local ones
    int ncj               = atom->Nclusters_local * SCLUSTER_SIZE;
    int cj_vec_base       = CJ_VECTOR_BASE_INDEX(ncj);
    MD_FLOAT* host_cl_x   = &atom->cl_x[cj_vec_base];
    MD_FLOAT* device_cl_x = &cuda_cl_x[cj_vec_base];
    memcpyToGPU(device_cl_x,
        host_cl_x,
        atom->Nclusters_ghost * CLUSTER_N * ATOM_DIM * sizeof(MD_FLOAT));
}

extern "C" void copyForceFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->cl_f,
        cuda_cl_f,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 * sizeof(MD_FLOAT));
}

extern "C" void copyForceToGPU(Atom* atom)
{
    memcpyToGPU(cuda_cl_f,
        atom->cl_f,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 * sizeof(MD_FLOAT));
}

extern "C" void growClustersCUDA(Atom* atom)
{
    if (cuda_cl_x) {
        cuda_cl_x             = (MD_FLOAT*)reallocateGPU(cuda_cl_x,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * ATOM_DIM *
                sizeof(MD_FLOAT));
        cuda_cl_v             = (MD_FLOAT*)reallocateGPU(cuda_cl_v,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
        cuda_cl_f             = (MD_FLOAT*)reallocateGPU(cuda_cl_f,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
        cuda_natoms           = (int*)reallocateGPU(cuda_natoms,
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
        cuda_jclusters_natoms = (int*)reallocateGPU(cuda_jclusters_natoms,
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
        cuda_numneigh         = (int*)reallocateGPU(cuda_numneigh,
            atom->Nclusters_max * sizeof(int));
        cuda_numneigh_inner   = (int*)reallocateGPU(cuda_numneigh_inner,
            atom->Nclusters_max * sizeof(int));

        freeHostPinned(natoms);
        freeHostPinned(ngatoms);

        natoms = (int*)allocateHostPinned(
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
        ngatoms = (int*)allocateHostPinned(
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
#if LJ_COMB_RULE == LJ_COMB_FULL
        cuda_cl_t = (int*)reallocateGPU(cuda_cl_t,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
        cuda_cl_sqrt_epsilon = (MD_FLOAT*)reallocateGPU(cuda_cl_sqrt_epsilon,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
        cuda_cl_sigma3       = (MD_FLOAT*)reallocateGPU(cuda_cl_sigma3,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
#endif
#if COULOMB_RULE != COULOMB_NONE
        cuda_cl_q = (MD_FLOAT*)reallocateGPU(cuda_cl_q,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
#endif
#if USE_EXCLUSIONS
        cuda_cl_mol = (int*)reallocateGPU(cuda_cl_mol,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
#endif
        cuda_cl_xq  = (MD_FLOAT4*)reallocateGPU(cuda_cl_xq,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT4));
        cuda_cl_par = (AtomParam*)reallocateGPU(cuda_cl_par,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(AtomParam));
    }
}

extern void growNeighborCUDA(Atom* atom, Neighbor* neighbor)
{
    cuda_neighbors       = (int*)reallocateGPU(cuda_neighbors,
        atom->Nclusters_max * neighbor->maxneighs * sizeof(int));
    cuda_neighbors_imask = (unsigned int*)reallocateGPU(cuda_neighbors_imask,
        atom->Nclusters_max * neighbor->maxneighs * sizeof(unsigned int));
}

/*
extern "C" void forwardCommCUDA(Comm* comm, Atom* atom, int iswap)
{
    int nrqst = 0, offset = 0, nsend = 0, nrecv = 0;
    int pbc[3];
    int size    = comm->forwardSize;
    int maxrqst = comm->numneigh;
    int* cuda_sendlist;
    int max_list_size;
    cuda_buf_send =  (MD_FLOAT*)allocateGPU(comm->max_send * sizeof(MD_FLOAT));
    cuda_buf_recv =  (MD_FLOAT*)allocateGPU(comm->max_recv * sizeof(MD_FLOAT));

    //use a single buffer and takes the highes list size to move list of cluster to send
    for (int ineigh = 0; ineigh < comm->numneigh; ineigh){
        max_list_size = comm->maxsendlist[ineigh];
    }
    //allocate the memory for the unique buffer
    cuda_sendlist = (int*)allocateGPU(max_list_size * sizeof(int));

    for (int ineigh = comm->sendfrom[iswap]; ineigh < comm->sendtill[iswap]; ineigh++) {
        offset  = comm->off_atom_send[ineigh];
        pbc[_x] = comm->pbc_x[ineigh];
        pbc[_y] = comm->pbc_y[ineigh];
        pbc[_z] = comm->pbc_z[ineigh];
        //copy lists into the buffer
        memcpyToGPU(cuda_sendlist, comm->sendlist[ineigh], comm->atom_send[ineigh] *
sizeof(int));

        const int threads_num = 64;
        dim3 block_size       = dim3(threads_num, 1, 1);
        dim3 grid_size = dim3(( comm->atom_send[ineigh] + threads_num - 1) / threads_num,
1, 1); packForwardCuda<<<grid_size, block_size>>>( cuda_cl_x, //MD_FLOAT*  -->need to be
in tye device cuda_jclusters_natoms,       //int*       -->need to be in tye device
                                            comm->atom_send[ineigh],     //int
                                            cuda_sendlist,               //int* -->need to
be in tye device cuda_buf_send[offset*size],  //MD_FLOAT*  -->need to be in tye device
                                            pbc[_x],                     //int
                                            pbc[_y],                     //int
                                            pbc[_z],                     //int
                                            atom->xprd,                  //MD_FLOAT
                                            atom->yrpd,                  //MD_FLOAT
                                            atom->zprd);                 //MD_FLOAT
    }


#ifdef _MPI
    MPI_Request requests[maxrqst];
    // Receives elements
    if (comm->othersend[iswap])
        for (int ineigh = comm->recvfrom[iswap]; ineigh < comm->recvtill[iswap];ineigh++)
{ offset = comm->off_atom_recv[ineigh] * size; nrecv  = comm->atom_recv[ineigh] * size;
            MPI_Irecv(&cuda_buf_recv[offset],
                    nrecv,
                    type,
                    comm->nrecv[ineigh],
                    0,
                    world,
                    &requests[nrqst++]);
        }

    // Send elements
    if (comm->othersend[iswap])
        for (int ineigh = comm->sendfrom[iswap]; ineigh < comm->sendtill[iswap]; ineigh++)
{ offset = comm->off_atom_send[ineigh] * size; nsend  = comm->atom_send[ineigh] * size;
            MPI_Send(&cuda_buf_send[offset], nsend, type, comm->nsend[ineigh], 0, world);
        }

    if (comm->othersend[iswap]) MPI_Waitall(nrqst, requests, MPI_STATUS_IGNORE);
#endif


    //unpack buffer
    for (int ineigh = comm->recvfrom[iswap]; ineigh < comm->recvtill[iswap]; ineigh++) {
        offset = comm->off_atom_recv[ineigh];
        MD_FLOAT *buf = (comm->othersend[iswap]) ? cuda_buf_recv : cuda_buf_send;
        const int threads_num = 64;
        dim3 block_size       = dim3(threads_num, 1, 1);
        dim3 grid_size = dim3(( comm->atom_recv[ineigh] + threads_num - 1) / threads_num,
1, 1);

        unpackForwardCuda<<<grid_size, block_size>>>(
                                        cuda_cl_x,                          //MD_FLOAT*
--> need to be in the device cuda_jclusters_natoms,              //int*      --> need to
be in the device comm->atom_recv[ineigh],            //int comm->firstrecv[iswap] +
offset,    //int &buf[offset * size]);               //MD_FLOAT* --> need to be in the
devic
    }
    GPUfree(cuda_sendlist);
    GPUfree(cuda_buf_recv);
    GPUfree(cuda_buf_send);
}

__global__ void packForwardCuda(MD_FLOAT* cuda_cl_x, int* cuda_jclusters_natoms, int nc,
int* cuda_list, MD_FLOAT* cuda_buf, int PBCx, int PBCY, int PBCz, MD_FLOAT xprd, MD_FLOAT
yrpd, MD_FLOAT zprd)
{
    unsigned int j = blockDim.x * blockIdx.x + threadIdx.x;
    if (j >= nc) return;

    int cj          = cuda_list[j];
    int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
    MD_FLOAT* cj_x  = &cuda_cl_x[cj_vec_base];
    int displ       = j * CLUSTER_N;

    for (int cjj = 0; cjj < cuda_jclusters_natoms[cj]; cjj++) {
        cuda_buf[3 * (displ + cjj) + 0] = cj_x[CL_X_INDEX(cjj)] +
                                    PBCx * xprd;
        cuda_buf[3 * (displ + cjj) + 1] = cj_x[CL_Y_INDEX(cjj)] +
                                    PBCy * yprd;
        cuda_buf[3 * (displ + cjj) + 2] = cj_x[CL_Z_INDEX(cjj)] +
                                    PBCz * zprd;
    }

    for (int cjj = cuda_jclusters_natoms[cj]; cjj < CLUSTER_N; cjj++) {
        cuda_buf[3 * (displ + cjj) + 0] = -1; // x
        cuda_buf[3 * (displ + cjj) + 1] = -1; // y
        cuda_buf[3 * (displ + cjj) + 2] = -1; // z
    }

}
__global__ void unpackForwardCuda(MD_FLOAT* cuda_cl_x, int* cuda_jclusters_natoms, int nc,
int c0, MD_FLOAT* cuda_buf)
{
    unsigned int j = blockDim.x * blockIdx.x + threadIdx.x;
    if (j >= nc) return;
    int cj = c0 + j;
    int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
    MD_FLOAT* cj_x  = &cuda_cl_x[cj_vec_base];
    int displ       = j * CLUSTER_N;

    for (int cjj = 0; cjj < cuda_jclusters_natoms[cj]; cjj++) {
        if (cj_x[CL_X_INDEX(cjj)] < INF)
            cj_x[CL_X_INDEX(cjj)] = cuda_buf[3 * (displ + cjj) + 0];
        if (cj_x[CL_Y_INDEX(cjj)] < INF)
            cj_x[CL_Y_INDEX(cjj)] = cuda_buf[3 * (displ + cjj) + 1];
        if (cj_x[CL_Z_INDEX(cjj)] < INF)
            cj_x[CL_Z_INDEX(cjj)] = cuda_buf[3 * (displ + cjj) + 2];
    }
}
*/
