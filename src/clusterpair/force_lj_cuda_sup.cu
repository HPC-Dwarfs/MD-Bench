#ifdef NBLIST_CSR
#error "Super-clustering CUDA kernels are not supported with NBLIST_DATA_LAYOUT=CSR"
#endif

#ifdef __HIP_PLATFORM_AMD__
#define FULL_WARP_MASK 0xffffffffffffffffULL
#else
#define FULL_WARP_MASK 0xffffffffU
#endif

// j-slices per i-supercluster (grid.y): 0 = chosen per force call, > 0 = fixed
#ifndef GPU_J_SLICES
#define GPU_J_SLICES 0
#endif

// K of S = K * sqrt(N_SM * n_tot) / N_sc per list type and thread mapping (GH200)
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
#ifndef K_SLICE_HALF
#define K_SLICE_HALF 0.6
#endif
#ifndef K_SLICE_FULL
#define K_SLICE_FULL 1.5
#endif
#else
#ifndef K_SLICE_HALF
#define K_SLICE_HALF 3.05
#endif
#ifndef K_SLICE_FULL
#define K_SLICE_FULL 2.5
#endif
#endif

// At most this many waves of resident blocks
#ifndef SLICE_MAX_WAVES
#define SLICE_MAX_WAVES 10
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

#include <nbdata_gpu.h>

// Packs cuda_cl_xq (and cuda_cl_par with pack_par), see force_lj_cuda.cu
extern "C" void launchPackAtomsCuda(Atom* atom, int pack_par);

// Neighbor list of the super-cluster kernels, {cj, sub-cluster mask} per entry
static int2* sup_nbl      = NULL;
static size_t sup_nbl_cap = 0;

/* Packs the list into {cj, mask}: mask bits 0..SCLUSTER_SIZE-1 are the sub-clusters to
 * visit, the bits above 1 + the sub-cluster that is cj itself (0 if none) */
__global__ void packNeighborsSup(const int* cuda_numneigh,
    const int* cuda_neighbors,
    const unsigned int* cuda_neighbors_imask,
    const int* cuda_border_map,
    int2* nbl,
    int Nclusters_local,
    int maxneighs,
    int half)
{
    const int sci       = blockIdx.x;
    const int ncj_local = Nclusters_local * SCLUSTER_SIZE;
    const int n         = cuda_numneigh[sci];

    for (int k = threadIdx.x; k < n; k += blockDim.x) {
        const int cj       = neighs_gpu(cuda_neighbors, sci, k, Nclusters_local, maxneighs);
        unsigned int imask = neighs_gpu(cuda_neighbors_imask,
                                 sci,
                                 k,
                                 Nclusters_local,
                                 maxneighs) &
                             ((1u << SCLUSTER_SIZE) - 1u);
        if (!half) {
            // A ghost cj is never one of the local sub-clusters
            const int d = cj - sci * SCLUSTER_SIZE;
            const int diag = (d >= 0 && d < SCLUSTER_SIZE) ? d : -1;
            neighs_gpu(nbl, sci, k, Nclusters_local, maxneighs) = make_int2(cj,
                (int)(imask | ((unsigned int)(diag + 1) << SCLUSTER_SIZE)));
            continue;
        }
        // Order a ghost cj by its real cluster
        const int order_cj = (cj < ncj_local) ? cj : cuda_border_map[cj - ncj_local];
        // Highest sub-cluster with ci <= order_cj
        const int last = order_cj - sci * SCLUSTER_SIZE;
        int diag       = -1;
        if (last < 0) {
            imask = 0;
        } else if (last < SCLUSTER_SIZE) {
            imask &= (2u << last) - 1u;
            diag = last;
        }
        neighs_gpu(nbl, sci, k, Nclusters_local, maxneighs) = make_int2(cj,
            (int)(imask | ((unsigned int)(diag + 1) << SCLUSTER_SIZE)));
    }
}

// Packs sup_nbl from cuda_neighbors and cuda_neighbors_imask
static void launchPackNeighborsSup(Atom* atom, Neighbor* neighbor)
{
    const size_t nnbl = (size_t)atom->Nclusters_local * neighbor->maxneighs;

    if (nnbl > sup_nbl_cap) {
        sup_nbl     = (int2*)reallocateGPU(sup_nbl, nnbl * sizeof(int2));
        sup_nbl_cap = nnbl;
    }
    if (nnbl > 0) {
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

extern "C" void packSupData(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    launchPackNeighborsSup(atom, neighbor);
}

extern "C" void freeSupData(void)
{
    if (sup_nbl) {
        GPUfree(sup_nbl);
    }
    sup_nbl     = NULL;
    sup_nbl_cap = 0;
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
    const NbForceArgs& a, int sci, MD_FLOAT4* sh_xq, AtomParam* sh_par)
{
    const int slot0 = SCI_SCALAR_BASE_INDEX(sci);
    const int tid   = threadIdx.y * CLUSTER_N + threadIdx.x; 
    for (int idx = tid; idx < SCLUSTER_SIZE * CLUSTER_M; idx += blockDim.x * blockDim.y) {
        MD_FLOAT4 xq = a.atoms.xq[slot0 + idx];
        xq.w *= a.cs.epsfac;
        sh_xq[idx]  = xq;
        sh_par[idx] = a.atoms.par[slot0 + idx];
    }
    __syncthreads();
}


__global__ void computeForceLJCudaSup_halfwarp(const NbForceArgs a)
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
    MD_FLOAT* sci_f = &a.atoms.f[SCI_VECTOR3_BASE_INDEX(sci)];
    MD_FLOAT3 fbuf[SCLUSTER_SIZE];

#pragma unroll
    for (int i = 0; i < SCLUSTER_SIZE; i++) {
        fbuf[i].x = (MD_FLOAT)0.0;
        fbuf[i].y = (MD_FLOAT)0.0;
        fbuf[i].z = (MD_FLOAT)0.0;
    }

    loadSciToShared(a, sci, sh_xq, sh_par);

    // In the sub-cluster that is cj only pairs cii < cjj interact
    const bool lowerPair   = cii < cjj;
    const int numneigh_sci = a.nbl.numneigh[sci];
    for (int k = slice; k < numneigh_sci; k += nslices) {
        const int2 nb =
            neighs_gpu(a.nbl.cj_mask, sci, k, a.atoms.nclusters_local, a.nbl.maxneighs);
        const int cj = nb.x;
        // Sub-clusters to visit and the one that is cj itself (-1 if none)
        const unsigned int imask = (unsigned int)nb.y & ((1u << SCLUSTER_SIZE) - 1u);
        const int diag_ci        = (int)((unsigned int)nb.y >> SCLUSTER_SIZE) - 1;
        const int slot_j         = CJ_SCALAR_BASE_INDEX(cj) + cjj;
        const MD_FLOAT4 xj       = a.atoms.xq[slot_j];
        const AtomParam pj       = a.atoms.par[slot_j];

#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        MD_FLOAT* cj_f = &a.atoms.f[CJ_VECTOR3_BASE_INDEX(cj)];
#else
        fcj_buf.x = fcj_buf.y = fcj_buf.z = (MD_FLOAT)0.0;
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
                const LjPair p     = ljPair(a.lj, pi, pj);

                if ((rsq < p.cutforcesq) & ((sci_ci != diag_ci) | lowerPair)) {
                    MD_FLOAT force = pairForceScalar(rsq,
                        p.sigma6,
                        p.epsilon,
                        xi.w,
                        xj.w,
                        pi.mol,
                        pj.mol,
                        a.cs);
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
        // Reduce the j-forces over the i-threads with shuffles
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
            atomicAdd(&a.atoms.f[(cj * CLUSTER_N + cjj) * 3 + cii], fcj_buf.x);
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
__global__ void computeForceLJCudaSup_fullwarp(const NbForceArgs a)
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
    MD_FLOAT* sci_f = &a.atoms.f[SCI_VECTOR3_BASE_INDEX(sci)];
    MD_FLOAT3 fbuf[SCLUSTER_SIZE];

#pragma unroll
    for (int i = 0; i < SCLUSTER_SIZE; i++) {
        fbuf[i].x = (MD_FLOAT)0.0;
        fbuf[i].y = (MD_FLOAT)0.0;
        fbuf[i].z = (MD_FLOAT)0.0;
    }

    loadSciToShared(a, sci, sh_xq, sh_par);

    // In the sub-cluster that is cj only the pair cii == cjj is excluded
    const bool offDiag     = cii != cjj;
    const int numneigh_sci = a.nbl.numneigh[sci];
    for (int k = slice; k < numneigh_sci; k += nslices) {
        const int2 nb =
            neighs_gpu(a.nbl.cj_mask, sci, k, a.atoms.nclusters_local, a.nbl.maxneighs);
        const int cj             = nb.x;
        const unsigned int imask = (unsigned int)nb.y & ((1u << SCLUSTER_SIZE) - 1u);
        const int diag_ci        = (int)((unsigned int)nb.y >> SCLUSTER_SIZE) - 1;
        const int slot_j         = CJ_SCALAR_BASE_INDEX(cj) + cjj;
        const MD_FLOAT4 xj       = a.atoms.xq[slot_j];
        const AtomParam pj       = a.atoms.par[slot_j];

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
                const LjPair p     = ljPair(a.lj, pi, pj);

                if ((rsq < p.cutforcesq) & ((sci_ci != diag_ci) | offDiag)) {
                    MD_FLOAT force = pairForceScalar(rsq,
                        p.sigma6,
                        p.epsilon,
                        xi.w,
                        xj.w,
                        pi.mol,
                        pj.mol,
                        a.cs);

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
        // cjj is the lane index within each group of CLUSTER_N lanes
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

// Adds the ghost forces to their real clusters (atomic: several ghosts per cluster)
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

// Resident blocks per SM of a force kernel (one wave = N_SM of these)
static int residentBlocksPerSM(const void* kernel)
{
    static const void* kernels[8];
    static int blocks[8];
    static int nkernels = 0;
    for (int i = 0; i < nkernels; i++) {
        if (kernels[i] == kernel) return blocks[i];
    }

    int n = 0;
    cuda_assert("cudaOccupancyMaxActiveBlocksPerMultiprocessor",
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n, kernel, CLUSTER_N * CLUSTER_M, 0));
    if (n < 1) n = 1;
    if (nkernels < 8) {
        kernels[nkernels] = kernel;
        blocks[nkernels]  = n;
        nkernels++;
    }
    return n;
}

// j-slices S = K * sqrt(N_SM * n_tot) / N_c, capped at SLICE_MAX_WAVES waves
extern "C" int chooseJSlicesK(
    const Atom* atom, const Neighbor* neighbor, double k_slice, const void* kernel)
{
    static int nsm = 0;
    if (nsm == 0) {                                  
        int dev = 0;
        cudaGetDevice(&dev);
        cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev);
    }
    const int blocks_per_sm = residentBlocksPerSM(kernel);

    const int nc = atom->Nclusters_local;             
    if (nc <= 0) return 1;

    long n_tot = 0;
    int  n_max = 1;
    for (int ci = 0; ci < nc; ci++) {
        const int n = neighbor->numneigh[ci];      
        n_tot += n;
        if (n > n_max) n_max = n;
    }

    int S = (int)lround(k_slice * sqrt((double)nsm) * sqrt((double)n_tot) / (double)nc);
    // At most SLICE_MAX_WAVES waves
    const long max_blocks = (long)SLICE_MAX_WAVES * nsm * blocks_per_sm;
    if (S > max_blocks / nc) S = (int)(max_blocks / nc);
    if (S < 1)     S = 1;
    if (S > n_max) S = n_max;
    if (S > 65535) S = 65535;                   
    return S;
}

// j-slices (grid.y) of this super-cluster force call, see GPU_J_SLICES
static int chooseJSlices(const Atom* atom, const Neighbor* neighbor)
{
#if GPU_J_SLICES > 0
    return GPU_J_SLICES;
#else
    return neighbor->half_neigh
               ? chooseJSlicesK(atom,
                     neighbor,
                     K_SLICE_HALF,
                     (const void*)computeForceLJCudaSup_halfwarp)
               : chooseJSlicesK(atom,
                     neighbor,
                     K_SLICE_FULL,
                     (const void*)computeForceLJCudaSup_fullwarp);
#endif
}

extern "C" double computeForceLJCudaSup(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCudaSup start\r\n");

    NbForceArgs args = makeNbForceArgs(param, atom, neighbor);
    args.nbl.cj_mask = sup_nbl;

    memsetGPU(cuda_cl_f,
        0,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));
    dim3 block_size = dim3(CLUSTER_N, CLUSTER_M, 1);
    dim3 grid_size  = dim3(atom->Nclusters_local, 1, 1);
    double S        = getTimeStamp();
    LIKWID_MARKER_START("force");

    // Repack x/y/z/q 
    launchPackAtomsCuda(atom, 0);

    grid_size.y = chooseJSlices(atom, neighbor);
    if (neighbor->half_neigh) {
        computeForceLJCudaSup_halfwarp<<<grid_size, block_size>>>(args);
    } else {
        computeForceLJCudaSup_fullwarp<<<grid_size, block_size>>>(args);
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

// Prunes the list: one block per super-cluster, thread (cii, cjj) tests one pair
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
