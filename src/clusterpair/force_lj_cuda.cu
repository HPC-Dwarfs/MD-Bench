/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

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

__global__ void computeForceLJCudaFullNeigh(const NbForceArgs nb)
{
    int ci         = blockIdx.x;
    int cii        = threadIdx.x;
    int cjj        = threadIdx.y;
    int ci_cj0     = CJ0_FROM_CI(ci);
    MD_FLOAT* ci_f = &nb.f[CI_VECTOR3_BASE_INDEX(ci)];
    int numneighs  = nb.numneigh[ci];
    MD_FLOAT3 fi   = { (MD_FLOAT)0.0, (MD_FLOAT)0.0, (MD_FLOAT)0.0 };

    const int slot_i   = CI_SCALAR_BASE_INDEX(ci) + cii;
    const MD_FLOAT4 xi = nb.xq[slot_i];
    const AtomParam pi = nb.par[slot_i];
    const MD_FLOAT qi  = nb.cs.epsfac * xi.w;

#pragma unroll 2
    for (int k = blockIdx.y; k < numneighs; k += gridDim.y) {
        const int cj = neighs_gpu(nb.cj, ci, k, nb.nclusters_local, nb.maxneighs);

        int cond;
#if CLUSTER_M == CLUSTER_N
        cond = ci_cj0 != cj || cii != cjj;
#elif CLUSTER_M < CLUSTER_N
        cond = ci_cj0 != cj || cii + CLUSTER_M * (ci & 0x1) != cjj;
#endif
        if (cond) {
            const int slot_j   = CJ_SCALAR_BASE_INDEX(cj) + cjj;
            const MD_FLOAT4 xj = nb.xq[slot_j];
            const AtomParam pj = nb.par[slot_j];
            MD_FLOAT delx      = xi.x - xj.x;
            MD_FLOAT dely      = xi.y - xj.y;
            MD_FLOAT delz      = xi.z - xj.z;
            MD_FLOAT rsq       = delx * delx + dely * dely + delz * delz;
            const LjPair p     = ljPair(nb, pi, pj);

            if (rsq < p.cutforcesq) {
                MD_FLOAT force = pairForceScalar(rsq,
                    p.sigma6,
                    p.epsilon,
                    qi,
                    xj.w,
                    pi.mol,
                    pj.mol,
                    nb.cs);

                fi.x += delx * force;
                fi.y += dely * force;
                fi.z += delz * force;
            }
        }
    }

    // Reduce the i-forces with warp shuffles
    addForceI(ci_f, cii, fi, cjj);
}

__global__ void computeForceLJCudaHalfNeigh(const NbForceArgs nb)
{
    int ci  = blockIdx.x;
    int cii = threadIdx.x;
    int cjj = threadIdx.y;
    MD_FLOAT3 fcj_buf;
    int ci_cj0     = CJ0_FROM_CI(ci);
    MD_FLOAT* ci_f = &nb.f[CI_VECTOR3_BASE_INDEX(ci)];
    int numneighs  = nb.numneigh[ci];
    MD_FLOAT3 fi   = { (MD_FLOAT)0.0, (MD_FLOAT)0.0, (MD_FLOAT)0.0 };

    const int slot_i   = CI_SCALAR_BASE_INDEX(ci) + cii;
    const MD_FLOAT4 xi = nb.xq[slot_i];
    const AtomParam pi = nb.par[slot_i];
    const MD_FLOAT qi  = nb.cs.epsfac * xi.w;

#pragma unroll 2
    for (int k = blockIdx.y; k < numneighs; k += gridDim.y) {
        const int cj = neighs_gpu(nb.cj, ci, k, nb.nclusters_local, nb.maxneighs);
        MD_FLOAT* cj_f = &nb.f[CJ_VECTOR3_BASE_INDEX(cj)];
        fcj_buf        = { 0.0, 0.0, 0.0 };

        int cond;
#if CLUSTER_M == CLUSTER_N
        cond = ci_cj0 != cj || cii < cjj;
#elif CLUSTER_M < CLUSTER_N
        cond = ci_cj0 != cj || cii + CLUSTER_M * (ci & 0x1) < cjj;
#endif
        if (cond) {
            const int slot_j   = CJ_SCALAR_BASE_INDEX(cj) + cjj;
            const MD_FLOAT4 xj = nb.xq[slot_j];
            const AtomParam pj = nb.par[slot_j];
            MD_FLOAT delx      = xi.x - xj.x;
            MD_FLOAT dely      = xi.y - xj.y;
            MD_FLOAT delz      = xi.z - xj.z;
            MD_FLOAT rsq       = delx * delx + dely * dely + delz * delz;
            const LjPair p     = ljPair(nb, pi, pj);

            if (rsq < p.cutforcesq) {
                MD_FLOAT force = pairForceScalar(rsq,
                    p.sigma6,
                    p.epsilon,
                    qi,
                    xj.w,
                    pi.mol,
                    pj.mol,
                    nb.cs);
                MD_FLOAT partial_force_x = delx * force;
                MD_FLOAT partial_force_y = dely * force;
                MD_FLOAT partial_force_z = delz * force;

                fcj_buf.x -= partial_force_x;
                fcj_buf.y -= partial_force_y;
                fcj_buf.z -= partial_force_z;

                fi.x += partial_force_x;
                fi.y += partial_force_y;
                fi.z += partial_force_z;
            }
        }

        addForceJ(cj_f, cjj, fcj_buf, cii);
    }

    // Reduce the i-forces with warp shuffles
    addForceI(ci_f, cii, fi, cjj);
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

static int chooseJSlices(const Atom* atom, const Neighbor* neighbor)
{
#if GPU_J_SLICES > 0
    return GPU_J_SLICES;
#else
    const double k_half = 3.3, k_full = 2.5;
    return neighbor->half_neigh
               ? chooseJSlicesK(atom,
                     neighbor,
                     k_half,
                     (const void*)computeForceLJCudaHalfNeigh)
               : chooseJSlicesK(atom,
                     neighbor,
                     k_full,
                     (const void*)computeForceLJCudaFullNeigh);
#endif
}

extern "C" double computeForceLJCuda(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCuda start\r\n");

    const NbForceArgs nb = makeNbForceArgs(param, atom, neighbor);

    // memsetGPU(cuda_cl_f, 0, atom->Nclusters_local * CLUSTER_M * 3 * sizeof(MD_FLOAT));
    memsetGPU(cuda_cl_f,
        0,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));
    dim3 block_size = dim3(CLUSTER_N, CLUSTER_M, 1);
    dim3 grid_size  = dim3(atom->Nclusters_local, 1, 1);
    double S        = getTimeStamp();
    LIKWID_MARKER_START("force");

    launchPackAtomsCuda(atom, 0);

    grid_size.y = chooseJSlices(atom, neighbor);
    if (neighbor->half_neigh) {
        computeForceLJCudaHalfNeigh<<<grid_size, block_size>>>(nb);
    } else {
        computeForceLJCudaFullNeigh<<<grid_size, block_size>>>(nb);
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

// ---------------------------------------------------------------------------------------
// Reference force kernel, selected at runtime with -ref.
// ---------------------------------------------------------------------------------------

__global__ void computeForceLJCudaRef_warp(MD_FLOAT* cuda_cl_x,
    MD_FLOAT* cuda_cl_f,
    int Nclusters_local,
    int* cuda_numneigh,
    int* cuda_neighs,
    int maxneighs,
    int half,
    MD_FLOAT* cuda_cl_q,
    int* cuda_cl_mol,
    CoulombScalar cs,
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq,
    MD_FLOAT* cuda_cl_sqrt_epsilon,
    MD_FLOAT* cuda_cl_sigma3
#else
    int* cuda_cl_t,
    MD_FLOAT* atom_cutforcesq,
    MD_FLOAT* atom_sigma6,
    MD_FLOAT* atom_epsilon,
    int ntypes
#endif
)
{
    int ci         = blockIdx.x;
    int cii        = threadIdx.y;
    int cjj        = threadIdx.x;
    int ci_cj0     = CJ0_FROM_CI(ci);
    MD_FLOAT* ci_x = &cuda_cl_x[CI_VECTOR_BASE_INDEX(ci)];
    MD_FLOAT* ci_f = &cuda_cl_f[CI_VECTOR3_BASE_INDEX(ci)];
    int numneighs  = cuda_numneigh[ci];
    const int i    = CI_SCALAR_BASE_INDEX(ci) + cii; 
    MD_FLOAT xtmp  = ci_x[CL_X_INDEX(cii)];
    MD_FLOAT ytmp  = ci_x[CL_Y_INDEX(cii)];
    MD_FLOAT ztmp  = ci_x[CL_Z_INDEX(cii)];
    MD_FLOAT3 fi   = { 0.0, 0.0, 0.0 };

    for (int k = 0; k < numneighs; k++) {
        int cj         = neighs_gpu(cuda_neighs, ci, k, Nclusters_local, maxneighs);
        MD_FLOAT* cj_x = &cuda_cl_x[CJ_VECTOR_BASE_INDEX(cj)];
        MD_FLOAT* cj_f = &cuda_cl_f[CJ_VECTOR3_BASE_INDEX(cj)];
        const int j    = CJ_SCALAR_BASE_INDEX(cj) + cjj; 

        bool skip;
        if (half) {
            skip = (ci_cj0 == cj && cii >= cjj);
        } else {
            skip = (ci_cj0 == cj && cii == cjj);
        }

        if (!skip) {
            MD_FLOAT delx = xtmp - cj_x[CL_X_INDEX(cjj)];
            MD_FLOAT dely = ytmp - cj_x[CL_Y_INDEX(cjj)];
            MD_FLOAT delz = ztmp - cj_x[CL_Z_INDEX(cjj)];
            MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
            MD_FLOAT sigma6  = cuda_cl_sigma3[i] * cuda_cl_sigma3[j];
            MD_FLOAT epsilon = cuda_cl_sqrt_epsilon[i] * cuda_cl_sqrt_epsilon[j];
#elif LJ_COMB_RULE == LJ_COMB_FULL
            int type_index      = cuda_cl_t[i] * ntypes + cuda_cl_t[j];
            MD_FLOAT cutforcesq = atom_cutforcesq[type_index];
            MD_FLOAT sigma6     = atom_sigma6[type_index];
            MD_FLOAT epsilon    = atom_epsilon[type_index];
#endif

            if (rsq < cutforcesq) {
                MD_FLOAT force = pairForceScalar(rsq,
                    sigma6,
                    epsilon,
                    cs.epsfac * chargeCuda(cuda_cl_q, i),
                    chargeCuda(cuda_cl_q, j),
                    moleculeCuda(cuda_cl_mol, i),
                    moleculeCuda(cuda_cl_mol, j),
                    cs);
                MD_FLOAT fx    = delx * force;
                MD_FLOAT fy    = dely * force;
                MD_FLOAT fz    = delz * force;

                fi.x += fx;
                fi.y += fy;
                fi.z += fz;

                if (half) {
                    // Newton's third law: the opposite force on atom j
                    atomicAdd(&cj_f[CL_X_INDEX_3D(cjj)], -fx);
                    atomicAdd(&cj_f[CL_Y_INDEX_3D(cjj)], -fy);
                    atomicAdd(&cj_f[CL_Z_INDEX_3D(cjj)], -fz);
                }
            }
        }
    }

    atomicAdd(&ci_f[CL_X_INDEX_3D(cii)], fi.x);
    atomicAdd(&ci_f[CL_Y_INDEX_3D(cii)], fi.y);
    atomicAdd(&ci_f[CL_Z_INDEX_3D(cii)], fi.z);
}

extern "C" double computeForceLJCudaRef(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCudaRef start\r\n");

#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif
    const CoulombScalar cs = coulombScalarInit(param);

    memsetGPU(cuda_cl_f,
        0,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));
    dim3 block_size = dim3(CLUSTER_N, CLUSTER_M, 1);
    dim3 grid_size  = dim3(atom->Nclusters_local, 1, 1);
    double S        = getTimeStamp();
    LIKWID_MARKER_START("force");

    computeForceLJCudaRef_warp<<<grid_size, block_size>>>(cuda_cl_x,
        cuda_cl_f,
        atom->Nclusters_local,
        cuda_numneigh_inner,
        cuda_neighbors,
        neighbor->maxneighs,
        neighbor->half_neigh,
        cuda_cl_q,
        cuda_cl_mol,
        cs,
#if LJ_COMB_RULE == LJ_COMB_SINGLE
        cutforcesq,
        sigma6,
        epsilon
#elif LJ_COMB_RULE == LJ_COMB_GEOM
        cutforcesq,
        cuda_cl_sqrt_epsilon,
        cuda_cl_sigma3
#else
        cuda_cl_t,
        cuda_cutforcesq,
        cuda_sigma6,
        cuda_epsilon,
        atom->ntypes
#endif
    );

    cuda_assert("computeForceLJCudaRef", cudaPeekAtLastError());
    cuda_assert("computeForceLJCudaRef", cudaDeviceSynchronize());
    LIKWID_MARKER_STOP("force");

#ifdef COMPUTE_STATS
    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        addStat(stats->calculated_forces, 1);
        addStat(stats->num_neighs, neighbor->numneigh[ci]);
    }
#endif

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJCudaRef stop\r\n");
    return E - S;
}
