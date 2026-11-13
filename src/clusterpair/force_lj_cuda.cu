/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

// gpusimple GPU kernels: one block per i-cluster and j-slice

#include <device.h>

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

#include <force_lj_cuda_common.h>

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

// Updates the ghost positions (mapping from setupPbc)
extern "C" void updatePbcCUDASimple(Atom* atom, Parameter* param)
{
    const int N           = atom->Nclusters_ghost;
    const int threads_num = 64;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

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

// Adds the ghost forces to their real clusters
extern "C" void reverseGhostForcesCUDASimple(Atom* atom, Parameter* param)
{
    const int N           = atom->Nclusters_ghost;
    const int threads_num = 64;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

    cudaReverseGhostForces_warp<<<grid_size, block_size>>>(cuda_cl_f,
        cuda_border_map,
        cuda_jclusters_natoms,
        atom->Nclusters_local,
        atom->Nclusters_ghost);
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

extern "C" void pruneNeighborCUDASimple(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("pruneNeighborCUDA start\n");

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

    // K of the j-slices per list type and thread mapping, tuned on an RTX 2080
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    const double k_half = 2.0, k_full = 2.0;
#else
    const double k_half = 3.3, k_full = 2.5;
#endif
    if (neighbor->half_neigh) {
        grid_size.y = chooseJSlicesK(atom,
            neighbor,
            k_half,
            (const void*)computeForceLJCudaHalfNeigh);
        computeForceLJCudaHalfNeigh<<<grid_size, block_size>>>(args);
    } else {
        grid_size.y = chooseJSlicesK(atom,
            neighbor,
            k_full,
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
