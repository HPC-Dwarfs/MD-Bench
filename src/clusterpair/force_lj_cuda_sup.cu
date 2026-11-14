#ifdef NBLIST_CSR
#error "Super-clustering CUDA kernels are not supported with NBLIST_DATA_LAYOUT=CSR"
#endif



extern "C" {

#include <stdio.h>
//---
#include <likwid-marker.h>
//---
#include <atom.h>
#include <device.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>
}

#include <force_lj_cuda_common.h>

static int2* sup_nbl = NULL; // {cj, sub-cluster mask}

// Packs sup_nbl as {cj, mask}: low bits = sub-clusters to visit, high bits = diag_ci + 1
__global__ void packNeighborsSup(const int* numneigh,
    const int* neighbors,
    const unsigned int* imasks,
    const int* border_map,
    int2* nbl,
    int nsc,
    int maxneighs,
    int half)
{
    const int sci       = blockIdx.x;
    const int ncj_local = nsc * SCLUSTER_SIZE; 
    const int n         = numneigh[sci];

    for (int k = threadIdx.x; k < n; k += blockDim.x) {
        const int cj       = neighs_gpu(neighbors, sci, k, nsc, maxneighs);
        unsigned int imask = neighs_gpu(imasks, sci, k, nsc, maxneighs);

        // half: get cj from border_map if it is a ghost cluster
        const int cj_real = (half && cj >= ncj_local) ? border_map[cj - ncj_local] : cj;

        int diag_ci = -1;
        for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
            const int ci = sci * SCLUSTER_SIZE + sci_ci;
            if (half && ci > cj_real) imask &= ~(1u << sci_ci); // only ci <= cj
            if (ci == cj_real) diag_ci = sci_ci;
        }

        imask &= (1u << SCLUSTER_SIZE) - 1u;
        const unsigned int mask = imask | ((unsigned int)(diag_ci + 1) << SCLUSTER_SIZE);
        neighs_gpu(nbl, sci, k, nsc, maxneighs) = make_int2(cj, (int)mask);
    }
}

// Packs sup_nbl from cuda_neighbors and cuda_neighbors_imask
extern "C" void launchPackNeighborsSup(Atom* atom, Neighbor* neighbor)
{
    static size_t cap = 0; // entries allocated in sup_nbl
    const size_t nnbl = (size_t)atom->Nclusters_local * neighbor->maxneighs;

    if (nnbl > 0) {
        if (sup_nbl == NULL || nnbl > cap) {
            sup_nbl = (int2*)reallocateGPU(sup_nbl, nnbl * sizeof(int2));
            cap     = nnbl;
        }
        packNeighborsSup<<<atom->Nclusters_local, 128>>>(cuda_numneigh,
            cuda_neighbors,
            cuda_neighbors_imask,
            cuda_border_map,
            sup_nbl,
            atom->Nclusters_local,
            neighbor->maxneighs,
            neighbor->half_neigh);
        cuda_assert("packNeighborsSup", cudaPeekAtLastError());
    }
    cuda_assert("packNeighborsSup", cudaDeviceSynchronize());
}

extern "C" void freeSupData(void)
{
    if (sup_nbl) {
        GPUfree(sup_nbl);
    }
    sup_nbl = NULL;
}

__global__ void cudaInitialIntegrateSup_warp(MD_FLOAT* cuda_cl_x,
    MD_FLOAT* cuda_cl_v,
    MD_FLOAT* cuda_cl_f,
    int Nclusters_local,
    MD_FLOAT dtforce,
    MD_FLOAT dt)
{

    int sci = blockIdx.x;
    int ci  = threadIdx.x;
    int cii = threadIdx.y;

    if (sci >= Nclusters_local) {
        return;
    }

    int sci_vec_base = SCI_VECTOR3_BASE_INDEX(sci);
    int i            = ci * CLUSTER_M + cii;
    MD_FLOAT* ci_x   = &cuda_cl_x[SCI_VECTOR_BASE_INDEX(sci)];
    MD_FLOAT* ci_v   = &cuda_cl_v[sci_vec_base];
    MD_FLOAT* ci_f   = &cuda_cl_f[sci_vec_base];

    ci_v[CL_X_INDEX_3D(i)] += dtforce * ci_f[CL_X_INDEX_3D(i)];
    ci_v[CL_Y_INDEX_3D(i)] += dtforce * ci_f[CL_Y_INDEX_3D(i)];
    ci_v[CL_Z_INDEX_3D(i)] += dtforce * ci_f[CL_Z_INDEX_3D(i)];
    ci_x[CL_X_INDEX(i)] += dt * ci_v[CL_X_INDEX_3D(i)];
    ci_x[CL_Y_INDEX(i)] += dt * ci_v[CL_Y_INDEX_3D(i)];
    ci_x[CL_Z_INDEX(i)] += dt * ci_v[CL_Z_INDEX_3D(i)];
}

extern "C" void cudaInitialIntegrateSup(Parameter* param, Atom* atom)
{
    dim3 block_size = dim3(SCLUSTER_SIZE, CLUSTER_M, 1);
    dim3 grid_size  = dim3(atom->Nclusters_local, 1, 1);

    cudaInitialIntegrateSup_warp<<<grid_size, block_size>>>(cuda_cl_x,
        cuda_cl_v,
        cuda_cl_f,
        atom->Nclusters_local,
        param->dtforce,
        param->dt);

    cuda_assert("cudaInitialIntegrateSup", cudaPeekAtLastError());
    cuda_assert("cudaInitialIntegrateSup", cudaDeviceSynchronize());
}

__global__ void cudaFinalIntegrateSup_warp(
    MD_FLOAT* cuda_cl_v, MD_FLOAT* cuda_cl_f, int Nclusters_local, MD_FLOAT dtforce)
{

    int sci = blockIdx.x;
    int ci  = threadIdx.x;
    int cii = threadIdx.y;

    if (sci >= Nclusters_local) {
        return;
    }

    int sci_vec_base = SCI_VECTOR3_BASE_INDEX(sci);
    int i            = ci * CLUSTER_M + cii;
    MD_FLOAT* ci_v   = &cuda_cl_v[sci_vec_base];
    MD_FLOAT* ci_f   = &cuda_cl_f[sci_vec_base];

    ci_v[CL_X_INDEX_3D(i)] += dtforce * ci_f[CL_X_INDEX_3D(i)];
    ci_v[CL_Y_INDEX_3D(i)] += dtforce * ci_f[CL_Y_INDEX_3D(i)];
    ci_v[CL_Z_INDEX_3D(i)] += dtforce * ci_f[CL_Z_INDEX_3D(i)];
}

extern "C" void cudaFinalIntegrateSup(Parameter* param, Atom* atom)
{
    dim3 block_size = dim3(SCLUSTER_SIZE, CLUSTER_M, 1);
    dim3 grid_size  = dim3(atom->Nclusters_local, 1, 1);

    cudaFinalIntegrateSup_warp<<<grid_size, block_size>>>(cuda_cl_v,
        cuda_cl_f,
        atom->Nclusters_local,
        param->dtforce);

    cuda_assert("cudaFinalIntegrateSup", cudaPeekAtLastError());
    cuda_assert("cudaFinalIntegrateSup", cudaDeviceSynchronize());
}

static __device__ __forceinline__ void loadSciToShared(
    const NbForceArgs& nb, int sci, MD_FLOAT4* sh_xq, AtomParam* sh_par)
{
    const int slot0 = SCI_SCALAR_BASE_INDEX(sci);
    const int tid   = threadIdx.y * CLUSTER_N + threadIdx.x; 
    for (int idx = tid; idx < SCLUSTER_SIZE * CLUSTER_M; idx += blockDim.x * blockDim.y) {
        MD_FLOAT4 xq = nb.atoms.xq[slot0 + idx];
        xq.w *= nb.cs.epsfac;
        sh_xq[idx]  = xq;
        sh_par[idx] = nb.atoms.par[slot0 + idx];
    }
    __syncthreads();
}

__global__ void computeForceLJCudaSup_halfwarp(const NbForceArgs nb)
{
    __shared__ MD_FLOAT4 sh_xq[SCLUSTER_SIZE * CLUSTER_M];
    __shared__ AtomParam sh_par[SCLUSTER_SIZE * CLUSTER_M];
    const int sci     = blockIdx.x;
    const int slice   = blockIdx.y;
    const int nslices = gridDim.y;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    const int cii = threadIdx.y;
    const int cjj = threadIdx.x;
#else
    const int cii = threadIdx.x;
    const int cjj = threadIdx.y;
    MD_FLOAT3 fcj_buf;
#endif
    MD_FLOAT* sci_f = &nb.atoms.f[SCI_VECTOR3_BASE_INDEX(sci)];
    MD_FLOAT3 fbuf[SCLUSTER_SIZE];

#pragma unroll
    for (int i = 0; i < SCLUSTER_SIZE; i++) {
        fbuf[i].x = (MD_FLOAT)0.0;
        fbuf[i].y = (MD_FLOAT)0.0;
        fbuf[i].z = (MD_FLOAT)0.0;
    }

    loadSciToShared(nb, sci, sh_xq, sh_par);

    // only cii < cjj interact
    const bool lowerPair   = cii < cjj;
    const int numneigh_sci = nb.nbl.numneigh[sci];
    for (int k = slice; k < numneigh_sci; k += nslices) {
        const int2 cj_mask = neighs_gpu(nb.nbl.cj_mask, sci, k, nb.atoms.nclusters_local, nb.nbl.maxneighs);
        const int cj = cj_mask.x;
        // imask: sub-clusters to compute; diag_ci: the sub-cluster that is cj (-1: none)
        const unsigned int imask = (unsigned int)cj_mask.y & ((1u << SCLUSTER_SIZE) - 1u);
        const int diag_ci        = (int)((unsigned int)cj_mask.y >> SCLUSTER_SIZE) - 1;
        const int slot_j         = CJ_SCALAR_BASE_INDEX(cj) + cjj;
        const MD_FLOAT4 xj       = nb.atoms.xq[slot_j];
        const AtomParam pj       = nb.atoms.par[slot_j];

#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        MD_FLOAT* cj_f = &nb.atoms.f[CJ_VECTOR3_BASE_INDEX(cj)];
#else
        fcj_buf = float3 { 0.0f, 0.0f, 0.0f };
#endif

#pragma unroll
        for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
            if (imask & (1u << sci_ci)) {
                const int ai       = sci_ci * CLUSTER_M + cii;
                const MD_FLOAT4 xi = sh_xq[ai];
                const AtomParam pi = sh_par[ai];
                MD_FLOAT delx      = xi.x - xj.x;
                MD_FLOAT dely      = xi.y - xj.y;
                MD_FLOAT delz      = xi.z - xj.z;
                MD_FLOAT rsq       = delx * delx + dely * dely + delz * delz;
                const LjPair p     = ljPair(nb.lj, pi, pj);

                if ((rsq < p.cutforcesq) & ((sci_ci != diag_ci) | lowerPair)) {
                    MD_FLOAT force = pairForceScalar(rsq,
                        p.sigma6,
                        p.epsilon,
                        xi.w,
                        xj.w,
                        pi.mol,
                        pj.mol,
                        nb.cs);
                    MD_FLOAT fx    = delx * force;
                    MD_FLOAT fy    = dely * force;
                    MD_FLOAT fz    = delz * force;

                    fbuf[sci_ci].x += fx;
                    fbuf[sci_ci].y += fy;
                    fbuf[sci_ci].z += fz;

#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
                    atomicAdd(&cj_f[CL_X_INDEX_3D(cjj)], -fx);
                    atomicAdd(&cj_f[CL_Y_INDEX_3D(cjj)], -fy);
                    atomicAdd(&cj_f[CL_Z_INDEX_3D(cjj)], -fz);
#else
                    fcj_buf.x -= fx;
                    fcj_buf.y -= fy;
                    fcj_buf.z -= fz;
#endif
                }
            }
        }

#ifndef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        auto mask = FULL_WARP_MASK;

        fcj_buf.x += __shfl_down_sync(mask, fcj_buf.x, 1);
        fcj_buf.y += __shfl_up_sync(mask, fcj_buf.y, 1);
        fcj_buf.z += __shfl_down_sync(mask, fcj_buf.z, 1);

        if (cii & 1) {
            fcj_buf.x = fcj_buf.y;
        }

        fcj_buf.x += __shfl_down_sync(mask, fcj_buf.x, 2);
        fcj_buf.z += __shfl_up_sync(mask, fcj_buf.z, 2);

        if (cii & 2) {
            fcj_buf.x = fcj_buf.z;
        }

        fcj_buf.x += __shfl_down_sync(mask, fcj_buf.x, 4);

        if (cii < 3) {
            atomicAdd(&nb.atoms.f[(cj * CLUSTER_N + cjj) * 3 + cii], fcj_buf.x);
        }
#endif
    }

#pragma unroll
    for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
        const int ai = sci_ci * CLUSTER_M + cii;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        atomicAdd(&sci_f[CL_X_INDEX_3D(ai)], fbuf[sci_ci].x);
        atomicAdd(&sci_f[CL_Y_INDEX_3D(ai)], fbuf[sci_ci].y);
        atomicAdd(&sci_f[CL_Z_INDEX_3D(ai)], fbuf[sci_ci].z);
#else
        MD_FLOAT fix = fbuf[sci_ci].x;
        MD_FLOAT fiy = fbuf[sci_ci].y;
        MD_FLOAT fiz = fbuf[sci_ci].z;
        auto mask    = FULL_WARP_MASK;
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

        if ((cjj & 3) < 3) {
            atomicAdd(&sci_f[ai * 3 + (cjj & 3)], fix);
        }
#endif
    }
}

// Full-list super-cluster kernel: as the half-list one, without j-forces
__global__ void computeForceLJCudaSup_fullwarp(const NbForceArgs nb)
{
    __shared__ MD_FLOAT4 sh_xq[SCLUSTER_SIZE * CLUSTER_M];
    __shared__ AtomParam sh_par[SCLUSTER_SIZE * CLUSTER_M];
    const int sci     = blockIdx.x;
    const int slice   = blockIdx.y;
    const int nslices = gridDim.y;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    const int cii = threadIdx.y;
    const int cjj = threadIdx.x;
#else
    const int cii = threadIdx.x;
    const int cjj = threadIdx.y;
#endif
    MD_FLOAT* sci_f = &nb.atoms.f[SCI_VECTOR3_BASE_INDEX(sci)];
    MD_FLOAT3 fbuf[SCLUSTER_SIZE];

#pragma unroll
    for (int i = 0; i < SCLUSTER_SIZE; i++) {
        fbuf[i].x = (MD_FLOAT)0.0;
        fbuf[i].y = (MD_FLOAT)0.0;
        fbuf[i].z = (MD_FLOAT)0.0;
    }

    loadSciToShared(nb, sci, sh_xq, sh_par);

    // cii == cjj is excluded
    const bool offDiag     = cii != cjj;
    const int numneigh_sci = nb.nbl.numneigh[sci];
    for (int k = slice; k < numneigh_sci; k += nslices) {
        const int2 cj_mask = neighs_gpu(nb.nbl.cj_mask, sci, k, nb.atoms.nclusters_local, nb.nbl.maxneighs);
        const int cj             = cj_mask.x;
        const unsigned int imask = (unsigned int)cj_mask.y & ((1u << SCLUSTER_SIZE) - 1u);
        const int diag_ci        = (int)((unsigned int)cj_mask.y >> SCLUSTER_SIZE) - 1;
        const int slot_j         = CJ_SCALAR_BASE_INDEX(cj) + cjj;
        const MD_FLOAT4 xj       = nb.atoms.xq[slot_j];
        const AtomParam pj       = nb.atoms.par[slot_j];

#pragma unroll
        for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
            if (imask & (1u << sci_ci)) {
                const int ai       = sci_ci * CLUSTER_M + cii;
                const MD_FLOAT4 xi = sh_xq[ai];
                const AtomParam pi = sh_par[ai];
                MD_FLOAT delx      = xi.x - xj.x;
                MD_FLOAT dely      = xi.y - xj.y;
                MD_FLOAT delz      = xi.z - xj.z;
                MD_FLOAT rsq       = delx * delx + dely * dely + delz * delz;
                const LjPair p     = ljPair(nb.lj, pi, pj);

                if ((rsq < p.cutforcesq) & ((sci_ci != diag_ci) | offDiag)) {
                    MD_FLOAT force = pairForceScalar(rsq,
                        p.sigma6,
                        p.epsilon,
                        xi.w,
                        xj.w,
                        pi.mol,
                        pj.mol,
                        nb.cs);

                    fbuf[sci_ci].x += delx * force;
                    fbuf[sci_ci].y += dely * force;
                    fbuf[sci_ci].z += delz * force;
                }
            }
        }
    }

#pragma unroll
    for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
        const int ai = sci_ci * CLUSTER_M + cii;
        MD_FLOAT fix = fbuf[sci_ci].x;
        MD_FLOAT fiy = fbuf[sci_ci].y;
        MD_FLOAT fiz = fbuf[sci_ci].z;
        auto mask    = FULL_WARP_MASK;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        
        for (int offset = CLUSTER_N / 2; offset > 0; offset /= 2) {
            fix += __shfl_down_sync(mask, fix, offset);
            fiy += __shfl_down_sync(mask, fiy, offset);
            fiz += __shfl_down_sync(mask, fiz, offset);
        }

        if (cjj == 0) {
#if GPU_J_SLICES == 1
            sci_f[CL_X_INDEX_3D(ai)] = fix;
            sci_f[CL_Y_INDEX_3D(ai)] = fiy;
            sci_f[CL_Z_INDEX_3D(ai)] = fiz;
#else
            atomicAdd(&sci_f[CL_X_INDEX_3D(ai)], fix);
            atomicAdd(&sci_f[CL_Y_INDEX_3D(ai)], fiy);
            atomicAdd(&sci_f[CL_Z_INDEX_3D(ai)], fiz);
#endif
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

        if ((cjj & 3) < 3) {
            atomicAdd(&sci_f[ai * 3 + (cjj & 3)], fix);
        }
#endif
    }
}

__global__ void cudaUpdatePbcSup_warp(MD_FLOAT* cuda_cl_x,
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

    int cg = blockDim.x * blockIdx.x + threadIdx.x;
    if (cg >= Nclusters_ghost) {
        return;
    }

    int ncj           = Nclusters_local * SCLUSTER_SIZE;
    int cj            = ncj + cg;
    int cj_vec_base   = CJ_VECTOR_BASE_INDEX(cj);
    int bmap_vec_base = CJ_VECTOR_BASE_INDEX(cuda_border_map[cg]);
    MD_FLOAT* cj_x    = &cuda_cl_x[cj_vec_base];
    MD_FLOAT* bmap_x  = &cuda_cl_x[bmap_vec_base];

    for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
        cj_x[CL_X_INDEX(cjj)] = bmap_x[CL_X_INDEX(cjj)] + cuda_PBCx[cg] * param_xprd;
        cj_x[CL_Y_INDEX(cjj)] = bmap_x[CL_Y_INDEX(cjj)] + cuda_PBCy[cg] * param_yprd;
        cj_x[CL_Z_INDEX(cjj)] = bmap_x[CL_Z_INDEX(cjj)] + cuda_PBCz[cg] * param_zprd;
    }
}

// Updates the ghost positions (mapping from setupPbc)
extern "C" void updatePbcCUDASup(Atom* atom, Parameter* param)
{
    const int N           = atom->Nclusters_ghost;
    const int threads_num = 64;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

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
    cuda_assert("cudaUpdatePbcSup", cudaPeekAtLastError());
    cuda_assert("cudaUpdatePbcSup", cudaDeviceSynchronize());
}

__global__ void cudaReverseGhostForcesSup_warp(MD_FLOAT* cuda_cl_f,
    int* cuda_border_map,
    int Nclusters_local,
    int Nclusters_ghost)
{
    int cg = blockDim.x * blockIdx.x + threadIdx.x;
    if (cg >= Nclusters_ghost) {
        return;
    }

    int ncj           = Nclusters_local * SCLUSTER_SIZE;
    int cj            = ncj + cg;
    int cj_vec_base   = CJ_VECTOR3_BASE_INDEX(cj);
    int bmap_vec_base = CJ_VECTOR3_BASE_INDEX(cuda_border_map[cg]);
    MD_FLOAT* cj_f    = &cuda_cl_f[cj_vec_base];
    MD_FLOAT* bmap_f  = &cuda_cl_f[bmap_vec_base];

    for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
        atomicAdd(&bmap_f[CL_X_INDEX_3D(cjj)], cj_f[CL_X_INDEX_3D(cjj)]);
        atomicAdd(&bmap_f[CL_Y_INDEX_3D(cjj)], cj_f[CL_Y_INDEX_3D(cjj)]);
        atomicAdd(&bmap_f[CL_Z_INDEX_3D(cjj)], cj_f[CL_Z_INDEX_3D(cjj)]);
    }
}


extern "C" void reverseGhostForcesCUDASup(Atom* atom, Parameter* param)
{
    const int N           = atom->Nclusters_ghost;
    const int threads_num = 64;
    dim3 block_size       = dim3(threads_num, 1, 1);
    dim3 grid_size        = dim3((N + threads_num - 1) / threads_num, 1, 1);

    cudaReverseGhostForcesSup_warp<<<grid_size, block_size>>>(cuda_cl_f,
        cuda_border_map,
        atom->Nclusters_local,
        atom->Nclusters_ghost);
    cuda_assert("cudaReverseGhostForcesSup", cudaPeekAtLastError());
    cuda_assert("cudaReverseGhostForcesSup", cudaDeviceSynchronize());
}

// selects the number of slices 
static int chooseJSlices(const Atom* atom, const Neighbor* neighbor)
{
#if GPU_J_SLICES > 0
    return GPU_J_SLICES;
#else
    
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    const double k_half = 0.6, k_full = 1.5;
#else
    const double k_half = 3.05, k_full = 2.5;
#endif
    return neighbor->half_neigh
               ? chooseJSlicesK(atom,
                     neighbor,
                     k_half,
                     (const void*)computeForceLJCudaSup_halfwarp)
               : chooseJSlicesK(atom,
                     neighbor,
                     k_full,
                     (const void*)computeForceLJCudaSup_fullwarp);
#endif
}

extern "C" double computeForceLJCudaSup(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCudaSup start\r\n");

    NbForceArgs nb = makeNbForceArgs(param, atom, neighbor);
    nb.nbl.cj_mask = sup_nbl;

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
        computeForceLJCudaSup_halfwarp<<<grid_size, block_size>>>(nb);
    } else {
        computeForceLJCudaSup_fullwarp<<<grid_size, block_size>>>(nb);
    }

    cuda_assert("computeForceLJCudaSup", cudaPeekAtLastError());
    cuda_assert("computeForceLJCudaSup", cudaDeviceSynchronize());

    LIKWID_MARKER_STOP("force");

#ifdef COMPUTE_STATS
    for (int sci = 0; sci < atom->Nclusters_local; sci++) {
        addStat(stats->calculated_forces, 1);
        addStat(stats->num_neighs, neighbor->numneigh[sci]);
    }
#endif

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJCudaSup stop\r\n");
    return E - S;
}

__global__ void cudaPruneNeighborSup(MD_FLOAT* cuda_cl_x,
    int* cuda_numneigh,
    int* cuda_numneigh_inner,
    int* cuda_neighbors,
    unsigned int* cuda_neighbors_imask,
    int Nclusters_local,
    int maxneighs,
    MD_FLOAT cutsq)
{
    __shared__ MD_FLOAT4 sh_sci_x[SCLUSTER_SIZE * CLUSTER_M];
    __shared__ unsigned int sh_mask;
    __shared__ int sh_lo;

    int sci = blockIdx.x;
    if (sci >= Nclusters_local) {
        return;
    }

    int cii = threadIdx.y;
    int cjj = threadIdx.x;
    int tid = cii * CLUSTER_N + cjj;

    MD_FLOAT* sci_x = &cuda_cl_x[SCI_VECTOR_BASE_INDEX(sci)];

    for (int idx = tid; idx < SCLUSTER_SIZE * CLUSTER_M; idx += blockDim.x * blockDim.y) {
        sh_sci_x[idx].x = sci_x[CL_X_INDEX(idx)];
        sh_sci_x[idx].y = sci_x[CL_Y_INDEX(idx)];
        sh_sci_x[idx].z = sci_x[CL_Z_INDEX(idx)];
    }

    if (tid == 0) {
        sh_lo = 0;
    }
    __syncthreads();

    const int numneighs = cuda_numneigh[sci];

    for (int hi = 0; hi < numneighs; hi++) {
        if (tid == 0) {
            sh_mask = 0;
        }
        __syncthreads();

        int cj         = neighs_gpu(cuda_neighbors, sci, hi, Nclusters_local, maxneighs);
        MD_FLOAT* cj_x = &cuda_cl_x[CJ_VECTOR_BASE_INDEX(cj)];
        MD_FLOAT xjtmp = cj_x[CL_X_INDEX(cjj)];
        MD_FLOAT yjtmp = cj_x[CL_Y_INDEX(cjj)];
        MD_FLOAT zjtmp = cj_x[CL_Z_INDEX(cjj)];

        unsigned int local_mask = 0;
#pragma unroll
        for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
            int ai        = sci_ci * CLUSTER_M + cii;
            MD_FLOAT delx = sh_sci_x[ai].x - xjtmp;
            MD_FLOAT dely = sh_sci_x[ai].y - yjtmp;
            MD_FLOAT delz = sh_sci_x[ai].z - zjtmp;
            if (delx * delx + dely * dely + delz * delz < cutsq) {
                local_mask |= (1u << sci_ci);
            }
        }

        if (local_mask) {
            atomicOr(&sh_mask, local_mask);
        }
        __syncthreads();

        if (tid == 0) {
            unsigned int imask = sh_mask;
            neighs_gpu(cuda_neighbors_imask, sci, hi, Nclusters_local, maxneighs) = imask;

            if (imask != 0) {
                if (hi != sh_lo) {
                    int t_cj = neighs_gpu(cuda_neighbors,
                        sci,
                        sh_lo,
                        Nclusters_local,
                        maxneighs);
                    unsigned int t_im = neighs_gpu(cuda_neighbors_imask,
                        sci,
                        sh_lo,
                        Nclusters_local,
                        maxneighs);
                    neighs_gpu(cuda_neighbors, sci, sh_lo, Nclusters_local, maxneighs) =
                        cj;
                    neighs_gpu(cuda_neighbors, sci, hi, Nclusters_local, maxneighs) =
                        t_cj;
                    neighs_gpu(cuda_neighbors_imask, sci, sh_lo, Nclusters_local, maxneighs) =
                        imask;
                    neighs_gpu(cuda_neighbors_imask, sci, hi, Nclusters_local, maxneighs) =
                        t_im;
                }
                sh_lo++;
            }
        }
        __syncthreads();
    }

    if (tid == 0) {
        cuda_numneigh_inner[sci] = sh_lo;
    }
}

extern "C" void pruneNeighborCUDASup(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("pruneNeighborCUDASup start\n");

    if (param->outer_skin <= 0.0) {
        // The build-time per-sub-cluster imask (copied with the list) stays valid.
        memcpyOnGPU(cuda_numneigh_inner,
            cuda_numneigh,
            atom->Nclusters_local * sizeof(int));
        return;
    }

    const MD_FLOAT cut_inner = param->cutforce + param->skin;
    const MD_FLOAT cutsq     = cut_inner * cut_inner;
    dim3 block_size          = dim3(CLUSTER_N, CLUSTER_M, 1);
    dim3 grid_size           = dim3(atom->Nclusters_local, 1, 1);

    cudaPruneNeighborSup<<<grid_size, block_size>>>(cuda_cl_x,
        cuda_numneigh,
        cuda_numneigh_inner,
        cuda_neighbors,
        cuda_neighbors_imask,
        atom->Nclusters_local,
        neighbor->maxneighs,
        cutsq);

    cuda_assert("cudaPruneNeighborSup", cudaPeekAtLastError());
    cuda_assert("cudaPruneNeighborSup", cudaDeviceSynchronize());
    // The pruning reordered cuda_neighbors/cuda_neighbors_imask in place
    launchPackNeighborsSup(atom, neighbor);
    DEBUG_MESSAGE("pruneNeighborCUDASup end\n");
}
