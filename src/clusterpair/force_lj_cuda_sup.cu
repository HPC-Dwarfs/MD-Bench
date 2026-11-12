#ifdef NBLIST_CSR
#error "Super-clustering CUDA kernels are not supported with NBLIST_DATA_LAYOUT=CSR"
#endif

#ifdef __HIP_PLATFORM_AMD__
#define FULL_WARP_MASK 0xffffffffffffffffULL
#else
#define FULL_WARP_MASK 0xffffffffU
#endif

/* j-slices per i-supercluster (grid.y) of the force kernels, config.mk's SUPERCLUSTER_J_SLICES
 * (-DGPU_J_SLICES=...; this fallback only applies outside that build system). 0: chosen per
 * force call from the neighbor list as
 *     S = K * sqrt(N_SM * n_tot) / N_sc    (n_tot: j-entries of all N_sc superclusters),
 * at most SLICE_MAX_WAVES waves of resident blocks; > 0: fixed S for both kernels. */
#ifndef GPU_J_SLICES
#define GPU_J_SLICES 0
#endif

/* K per list type and thread mapping, tuned on GH200 (argon 1k-1M atoms, r 0.9-2.8: mean loss
 * against the best S of each case 0.9-1.5%). With the inverse mapping the half kernel adds every
 * j-force atomically, which wants few slices. */
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

/* A wave is the number of blocks all SMs hold at once. Beyond about 10 waves the tail is already
 * small and further slices only add per-block overhead (half list, 1M atoms: one block per
 * supercluster, about 7 waves, is 4-6% faster than the 13-34 waves of the formula alone) */
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

#include <pair_cuda.h>

extern "C" {
extern MD_FLOAT* cuda_cl_x;
extern MD_FLOAT* cuda_cl_v;
extern MD_FLOAT* cuda_cl_f;
extern int* cuda_neighbors;
extern unsigned int* cuda_neighbors_imask;
extern int* cuda_numneigh;
extern int* cuda_numneigh_inner;
extern int* cuda_natoms;
extern int* natoms;
extern int* ngatoms;
extern int* cuda_border_map;
extern int* cuda_jclusters_natoms;
extern MD_FLOAT *cuda_bbminx, *cuda_bbmaxx;
extern MD_FLOAT *cuda_bbminy, *cuda_bbmaxy;
extern MD_FLOAT *cuda_bbminz, *cuda_bbmaxz;
extern int *cuda_PBCx, *cuda_PBCy, *cuda_PBCz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
extern int* cuda_cl_t;
extern MD_FLOAT* cuda_cutforcesq;
extern MD_FLOAT* cuda_sigma6;
extern MD_FLOAT* cuda_epsilon;
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
extern MD_FLOAT* cuda_cl_sqrt_epsilon;
extern MD_FLOAT* cuda_cl_sigma3;
#endif
extern MD_FLOAT* cuda_cl_q;
extern int* cuda_cl_mol;
}

/* Packed per-atom data of the super-cluster kernels (computeForceLJCudaSup_halfwarp,
 * computeForceLJCudaSup_fullwarp), like the xq/ljComb arrays of the GROMACS NBNxM GPU
 * kernels: one 16-byte load per atom instead of one scalar load per field, and cj with its
 * sub-cluster mask in one 8-byte load. */
/* Indexed by the scalar cluster slot CJ_SCALAR_BASE_INDEX(cj) + cjj:
 *   sup_xq:  x, y, z, q (repacked every step, the positions move on the device)
 *   sup_par: LJ parameters, molecule index (int bits), 0 (packed on list upload):
 *            geometric: sqrt(epsilon), sigma^3; full: atom type (int bits), 0; single: 0, 0
 * sup_nbl has the layout of cuda_neighbors, {cj, imask} per entry (packed on list upload
 * and after pruning on the GPU) */
static MD_FLOAT4* sup_xq   = NULL;
static MD_FLOAT4* sup_par  = NULL;
static int2* sup_nbl       = NULL;
static size_t sup_atom_cap = 0;
static size_t sup_nbl_cap  = 0;

/* The molecule index travels in a float slot of sup_par, bit-exact */
static __device__ __forceinline__ MD_FLOAT packMoleculeCuda(int mol)
{
#if PRECISION == 1
    return __int_as_float(mol);
#else
    return __longlong_as_double((long long)mol);
#endif
}

static __device__ __forceinline__ int unpackMoleculeCuda(MD_FLOAT m)
{
#if PRECISION == 1
    return __float_as_int(m);
#else
    return (int)__double_as_longlong(m);
#endif
}

__global__ void packXqSup(
    const MD_FLOAT* cuda_cl_x, const MD_FLOAT* cuda_cl_q, MD_FLOAT4* xq, int nclusters)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nclusters * CLUSTER_N) {
        return;
    }

    const int c       = a / CLUSTER_N;
    const int i       = a % CLUSTER_N;
    const MD_FLOAT* x = &cuda_cl_x[CJ_VECTOR_BASE_INDEX(c)];
    const int sca     = CJ_SCALAR_BASE_INDEX(c) + i;
    MD_FLOAT4 v;
    v.x     = x[CL_X_INDEX(i)];
    v.y     = x[CL_Y_INDEX(i)];
    v.z     = x[CL_Z_INDEX(i)];
    v.w     = chargeCuda(cuda_cl_q, sca);
    xq[sca] = v;
}

/* The atom type travels in a float slot of sup_par, bit-exact, like the molecule */
#define packTypeCuda   packMoleculeCuda
#define unpackTypeCuda unpackMoleculeCuda

__global__ void packParamsSup(
#if LJ_COMB_RULE == LJ_COMB_GEOM
    const MD_FLOAT* cuda_cl_sqrt_epsilon,
    const MD_FLOAT* cuda_cl_sigma3,
#elif LJ_COMB_RULE == LJ_COMB_FULL
    const int* cuda_cl_t,
#endif
    const int* cuda_cl_mol,
    MD_FLOAT4* par,
    int nslots)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nslots) {
        return;
    }

    MD_FLOAT4 v;
#if LJ_COMB_RULE == LJ_COMB_GEOM
    v.x = cuda_cl_sqrt_epsilon[a];
    v.y = cuda_cl_sigma3[a];
#elif LJ_COMB_RULE == LJ_COMB_FULL
    v.x = packTypeCuda(cuda_cl_t[a]);
    v.y = (MD_FLOAT)0.0;
#else
    v.x = (MD_FLOAT)0.0;
    v.y = (MD_FLOAT)0.0;
#endif
    v.z    = packMoleculeCuda(moleculeCuda(cuda_cl_mol, a));
    v.w    = (MD_FLOAT)0.0;
    par[a] = v;
}

/* One block per i-super-cluster. Besides packing {cj, imask}, this prepares the
 * self-interaction handling once per list instead of in every force call, so that the
 * kernels only test the warp-uniform imask bit per sub-cluster (as GROMACS does).
 * Half list (half != 0):
 *  - clears the sub-clusters with ci > order_cj (the CPU builder already drops them, the
 *    GPU pruning does not),
 *  - stores 1 + the index of the sub-cluster ci == order_cj in the bits above the mask
 *    (0 if cj is not one of the sub-clusters or a periodic image of one), where the
 *    kernel keeps only the pairs cii < cjj.
 * Full list: the mask is kept, the bits above it store 1 + the index of the sub-cluster
 * ci == cj (0 if none), where the kernel drops the pairs cii == cjj. */
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
        /* Resolve a ghost cj to its real mirror before ordering, else a
         * periodic-boundary pair gets counted from both sides. */
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

static int numSupClusters(const Atom* atom)
{
    return atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost;
}

/* Packs sup_par (with_params) and sup_nbl from the arrays copyDataToCUDADevice() uploaded,
 * growing the packed buffers when the clusters or the list grew */
static void packSup(Atom* atom, Neighbor* neighbor, int with_params)
{
    const int threads    = 256;
    const size_t nslots  = (size_t)numSupClusters(atom) * CLUSTER_N;
    const size_t nnbl    = (size_t)atom->Nclusters_local * neighbor->maxneighs;

    if (nslots > sup_atom_cap) {
        if (sup_xq) {
            GPUfree(sup_xq);
            GPUfree(sup_par);
        }
        sup_xq       = (MD_FLOAT4*)allocateGPU(nslots * sizeof(MD_FLOAT4));
        sup_par      = (MD_FLOAT4*)allocateGPU(nslots * sizeof(MD_FLOAT4));
        sup_atom_cap = nslots;
        with_params  = 1;
    }
    if (nnbl > sup_nbl_cap) {
        if (sup_nbl) {
            GPUfree(sup_nbl);
        }
        sup_nbl     = (int2*)allocateGPU(nnbl * sizeof(int2));
        sup_nbl_cap = nnbl;
    }

    if (with_params && nslots > 0) {
        packParamsSup<<<(nslots + threads - 1) / threads, threads>>>(
#if LJ_COMB_RULE == LJ_COMB_GEOM
            cuda_cl_sqrt_epsilon,
            cuda_cl_sigma3,
#elif LJ_COMB_RULE == LJ_COMB_FULL
            cuda_cl_t,
#endif
            cuda_cl_mol,
            sup_par,
            (int)nslots);
        cuda_assert("packParamsSup", cudaPeekAtLastError());
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
    cuda_assert("packSup", cudaDeviceSynchronize());
}

/* Called by copyDataToCUDADevice() after every list upload */
extern "C" void packSupData(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    if (param->super_clustering) {
        packSup(atom, neighbor, 1);
    }
}

extern "C" void freeSupData(void)
{
    if (sup_xq) {
        GPUfree(sup_xq);
        GPUfree(sup_par);
    }
    if (sup_nbl) {
        GPUfree(sup_nbl);
    }
    sup_xq = sup_par = NULL;
    sup_nbl          = NULL;
    sup_atom_cap = sup_nbl_cap = 0;
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

/* LJ parameters of the super-cluster force kernels: the cut-off (and for the single rule
 * sigma^6/epsilon) by value, the per-type-pair tables for the full rule. The geometric rule
 * takes sqrt(epsilon) and sigma^3 from the packed per-atom data. */
typedef struct {
#if LJ_COMB_RULE == LJ_COMB_FULL
    const MD_FLOAT* cutforcesq;
    const MD_FLOAT* sigma6;
    const MD_FLOAT* epsilon;
    int ntypes;
#else
    MD_FLOAT cutforcesq;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT sigma6;
    MD_FLOAT epsilon;
#endif
#endif
} SupLJ;

/* Cut-off of the pair (i-atom ai in shared memory, j-atom of type type_j). The full rule
 * looks it up per type pair and returns the table index in *tidx for supLJParams() */
static __device__ __forceinline__ MD_FLOAT supCutoffSq(
    const SupLJ& lj, const MD_FLOAT4* sh_par, int ai, int type_j, int* tidx)
{
#if LJ_COMB_RULE == LJ_COMB_FULL
    *tidx = unpackTypeCuda(sh_par[ai].x) * lj.ntypes + type_j;
    return lj.cutforcesq[*tidx];
#else
    (void)sh_par;
    (void)ai;
    (void)type_j;
    *tidx = 0;
    return lj.cutforcesq;
#endif
}

static __device__ __forceinline__ void supLJParams(const SupLJ& lj,
    const MD_FLOAT4& pari,
    const MD_FLOAT4& parj,
    int tidx,
    MD_FLOAT* sigma6,
    MD_FLOAT* epsilon)
{
#if LJ_COMB_RULE == LJ_COMB_GEOM
    (void)lj;
    (void)tidx;
    *sigma6  = pari.y * parj.y;
    *epsilon = pari.x * parj.x;
#elif LJ_COMB_RULE == LJ_COMB_FULL
    (void)pari;
    (void)parj;
    *sigma6  = lj.sigma6[tidx];
    *epsilon = lj.epsilon[tidx];
#else
    (void)pari;
    (void)parj;
    (void)tidx;
    *sigma6  = lj.sigma6;
    *epsilon = lj.epsilon;
#endif
}

/* Half-list super-cluster kernel. Per j-cluster one {cj, imask} load and two 16-byte atom
 * loads (x/y/z/q and LJ parameters/molecule); the i-atom data is loaded into shared memory
 * once per block, a pair reads it with 16-byte shared-memory loads. The sub-cluster loop
 * only tests the warp-uniform imask bit, the half-list ordering is prepared by
 * packNeighborsSup() and the diagonal sub-cluster is handled in the cut-off test.
 * With the GROMACS thread mapping (cii = threadIdx.x) the j-forces are reduced over the
 * i-threads with shuffles, with SUPERCLUSTER_INVERSE_THREAD_MAPPING every pair adds its
 * j-force atomically. */
__global__ void computeForceLJCudaSup_halfwarp(const MD_FLOAT4* cl_xq,
    const MD_FLOAT4* cl_par,
    MD_FLOAT* cuda_cl_f,
    int Nclusters_local,
    const int* cuda_numneigh,
    const int2* cuda_nbl,
    int maxneighs,
    CoulombScalar cs,
    SupLJ lj)
{
    /* i-atoms: x/y/z/q, then the LJ parameters and molecule (like ljcpib in the GROMACS
     * kernel). One array, so both are addressed from one base register */
    __shared__ MD_FLOAT4 sh_sci[2 * SCLUSTER_SIZE * CLUSTER_M];
    MD_FLOAT4* sh_sci_x   = sh_sci;
    MD_FLOAT4* sh_sci_par = sh_sci + SCLUSTER_SIZE * CLUSTER_M;
    const int sci         = blockIdx.x;
    const int slice       = blockIdx.y;
    const int nslices     = gridDim.y;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    const int cii = threadIdx.y;
    const int cjj = threadIdx.x;
#else
    const int cii = threadIdx.x;
    const int cjj = threadIdx.y;
    MD_FLOAT3 fcj_buf;
#endif
    const int tid          = threadIdx.y * CLUSTER_N + threadIdx.x; // blockDim.x == CLUSTER_N
    const int sci_sca_base = SCI_SCALAR_BASE_INDEX(sci);
    MD_FLOAT* sci_f        = &cuda_cl_f[SCI_VECTOR3_BASE_INDEX(sci)];
    MD_FLOAT3 fbuf[SCLUSTER_SIZE];

#pragma unroll
    for (int i = 0; i < SCLUSTER_SIZE; i++) {
        fbuf[i].x = (MD_FLOAT)0.0;
        fbuf[i].y = (MD_FLOAT)0.0;
        fbuf[i].z = (MD_FLOAT)0.0;
    }

    // i-atoms: position and charge (scaled by epsfac once here) and the parameters,
    // invariant over the whole j-loop
    for (int idx = tid; idx < SCLUSTER_SIZE * CLUSTER_M; idx += blockDim.x * blockDim.y) {
        MD_FLOAT4 xqi = cl_xq[sci_sca_base + idx];
        xqi.w *= cs.epsfac;
        sh_sci_x[idx]   = xqi;
        sh_sci_par[idx] = cl_par[sci_sca_base + idx];
    }

    __syncthreads();

    // In the sub-cluster that is cj itself only the pairs cii < cjj interact (half list)
    const bool lowerPair   = cii < cjj;
    const int numneigh_sci = cuda_numneigh[sci];
    for (int k = slice; k < numneigh_sci; k += nslices) {
        const int2 nb = neighs_gpu(cuda_nbl, sci, k, Nclusters_local, maxneighs);
        const int cj  = nb.x;
        /* Sub-clusters to visit (packNeighborsSup() already dropped ci > cj) and the one
         * that is cj itself, -1 if none: both are warp-uniform */
        const unsigned int imask = (unsigned int)nb.y & ((1u << SCLUSTER_SIZE) - 1u);
        const int diag_ci        = (int)((unsigned int)nb.y >> SCLUSTER_SIZE) - 1;
        const int aj             = CJ_SCALAR_BASE_INDEX(cj) + cjj;
        const MD_FLOAT4 xqj      = cl_xq[aj];
        const MD_FLOAT4 parj     = cl_par[aj];
        const int mol_j          = unpackMoleculeCuda(parj.z);
        const int type_j         = unpackTypeCuda(parj.x); // used by the full rule only

#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
        MD_FLOAT* cj_f = &cuda_cl_f[CJ_VECTOR3_BASE_INDEX(cj)];
#else
        fcj_buf.x = fcj_buf.y = fcj_buf.z = (MD_FLOAT)0.0;
#endif

#pragma unroll
        for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
            if (imask & (1u << sci_ci)) {
                const int ai        = sci_ci * CLUSTER_M + cii;
                const MD_FLOAT4 xqi = sh_sci_x[ai];
                MD_FLOAT delx       = xqi.x - xqj.x;
                MD_FLOAT dely       = xqi.y - xqj.y;
                MD_FLOAT delz       = xqi.z - xqj.z;
                MD_FLOAT rsq        = delx * delx + dely * dely + delz * delz;
                int tidx;
                const MD_FLOAT cutforcesq = supCutoffSq(lj, sh_sci_par, ai, type_j, &tidx);

                if ((rsq < cutforcesq) & ((sci_ci != diag_ci) | lowerPair)) {
                    const MD_FLOAT4 pari = sh_sci_par[ai];
                    MD_FLOAT sigma6, epsilon;
                    supLJParams(lj, pari, parj, tidx, &sigma6, &epsilon);
                    MD_FLOAT force = pairForceScalar(rsq,
                        sigma6,
                        epsilon,
                        xqi.w,
                        xqj.w,
                        unpackMoleculeCuda(pari.z),
                        mol_j,
                        cs);
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
        // Reduce the j-forces over the i-threads (as reduce_force_j_warp_shfl in GROMACS)
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
            atomicAdd(&cuda_cl_f[(cj * CLUSTER_N + cjj) * 3 + cii], fcj_buf.x);
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

/* Full-list counterpart of computeForceLJCudaSup_halfwarp: the same loads and sub-cluster
 * loop, the self pair (ci == cj, cii == cjj) dropped in the cut-off test with the
 * sub-cluster packNeighborsSup() marked. No j-forces. */
__global__ void computeForceLJCudaSup_fullwarp(const MD_FLOAT4* cl_xq,
    const MD_FLOAT4* cl_par,
    MD_FLOAT* cuda_cl_f,
    int Nclusters_local,
    const int* cuda_numneigh,
    const int2* cuda_nbl,
    int maxneighs,
    CoulombScalar cs,
    SupLJ lj)
{
    __shared__ MD_FLOAT4 sh_sci[2 * SCLUSTER_SIZE * CLUSTER_M];
    MD_FLOAT4* sh_sci_x   = sh_sci;
    MD_FLOAT4* sh_sci_par = sh_sci + SCLUSTER_SIZE * CLUSTER_M;
    const int sci         = blockIdx.x;
    const int slice       = blockIdx.y;
    const int nslices     = gridDim.y;
#ifdef SUPERCLUSTER_INVERSE_THREAD_MAPPING
    const int cii = threadIdx.y;
    const int cjj = threadIdx.x;
#else
    const int cii = threadIdx.x;
    const int cjj = threadIdx.y;
#endif
    const int tid          = threadIdx.y * CLUSTER_N + threadIdx.x; // blockDim.x == CLUSTER_N
    const int sci_sca_base = SCI_SCALAR_BASE_INDEX(sci);
    MD_FLOAT* sci_f        = &cuda_cl_f[SCI_VECTOR3_BASE_INDEX(sci)];
    MD_FLOAT3 fbuf[SCLUSTER_SIZE];

#pragma unroll
    for (int i = 0; i < SCLUSTER_SIZE; i++) {
        fbuf[i].x = (MD_FLOAT)0.0;
        fbuf[i].y = (MD_FLOAT)0.0;
        fbuf[i].z = (MD_FLOAT)0.0;
    }

    for (int idx = tid; idx < SCLUSTER_SIZE * CLUSTER_M; idx += blockDim.x * blockDim.y) {
        MD_FLOAT4 xqi = cl_xq[sci_sca_base + idx];
        xqi.w *= cs.epsfac;
        sh_sci_x[idx]   = xqi;
        sh_sci_par[idx] = cl_par[sci_sca_base + idx];
    }

    __syncthreads();

    // In the sub-cluster that is cj itself only the pair cii == cjj is excluded
    const bool offDiag     = cii != cjj;
    const int numneigh_sci = cuda_numneigh[sci];
    for (int k = slice; k < numneigh_sci; k += nslices) {
        const int2 nb            = neighs_gpu(cuda_nbl, sci, k, Nclusters_local, maxneighs);
        const int cj             = nb.x;
        const unsigned int imask = (unsigned int)nb.y & ((1u << SCLUSTER_SIZE) - 1u);
        const int diag_ci        = (int)((unsigned int)nb.y >> SCLUSTER_SIZE) - 1;
        const int aj             = CJ_SCALAR_BASE_INDEX(cj) + cjj;
        const MD_FLOAT4 xqj      = cl_xq[aj];
        const MD_FLOAT4 parj     = cl_par[aj];
        const int mol_j          = unpackMoleculeCuda(parj.z);
        const int type_j         = unpackTypeCuda(parj.x); // used by the full rule only

#pragma unroll
        for (int sci_ci = 0; sci_ci < SCLUSTER_SIZE; sci_ci++) {
            if (imask & (1u << sci_ci)) {
                const int ai        = sci_ci * CLUSTER_M + cii;
                const MD_FLOAT4 xqi = sh_sci_x[ai];
                MD_FLOAT delx       = xqi.x - xqj.x;
                MD_FLOAT dely       = xqi.y - xqj.y;
                MD_FLOAT delz       = xqi.z - xqj.z;
                MD_FLOAT rsq        = delx * delx + dely * dely + delz * delz;
                int tidx;
                const MD_FLOAT cutforcesq = supCutoffSq(lj, sh_sci_par, ai, type_j, &tidx);

                if ((rsq < cutforcesq) & ((sci_ci != diag_ci) | offDiag)) {
                    const MD_FLOAT4 pari = sh_sci_par[ai];
                    MD_FLOAT sigma6, epsilon;
                    supLJParams(lj, pari, parj, tidx, &sigma6, &epsilon);
                    MD_FLOAT force = pairForceScalar(rsq,
                        sigma6,
                        epsilon,
                        xqi.w,
                        xqj.w,
                        unpackMoleculeCuda(pari.z),
                        mol_j,
                        cs);

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

/* fold ghost-cluster reaction forces back onto the real cluster they mirror,
 * entirely on-device; a real cluster can be mirrored by several ghost
 * clusters (one per PBC image), so the accumulation must be atomic */
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

// Resident blocks per SM of a force kernel (registers, shared memory): one wave is N_SM of
// these. Queried once per kernel, so it follows the kernel's register count.
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

/* j-slices (grid.y) S = K * sqrt(N_SM * n_tot) / N_c of a force kernel with one block of
 * CLUSTER_N x CLUSTER_M threads per i-(super)cluster (N_c of them, n_tot j-entries in all), at
 * most SLICE_MAX_WAVES waves of the kernel's resident blocks. Also used by the gpusimple
 * kernels (force_lj_cuda.cu). */
extern "C" int chooseJSlicesK(
    const Atom* atom, const Neighbor* neighbor, double k_slice, const void* kernel)
{
    static int nsm = 0;
    if (nsm == 0) {                                   // N_SM
        int dev = 0;
        cudaGetDevice(&dev);
        cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev);
    }
    const int blocks_per_sm = residentBlocksPerSM(kernel);

    const int nc = atom->Nclusters_local;             // N_c
    if (nc <= 0) return 1;

    long n_tot = 0;
    int  n_max = 1;
    for (int ci = 0; ci < nc; ci++) {
        const int n = neighbor->numneigh[ci];         // total number of neighbors for this cluster
        n_tot += n;
        if (n > n_max) n_max = n;
    }

    int S = (int)lround(k_slice * sqrt((double)nsm) * sqrt((double)n_tot) / (double)nc);
    // At most SLICE_MAX_WAVES waves: the formula keeps adding slices as the system grows
    const long max_blocks = (long)SLICE_MAX_WAVES * nsm * blocks_per_sm;
    if (S > max_blocks / nc) S = (int)(max_blocks / nc);
    if (S < 1)     S = 1;
    if (S > n_max) S = n_max;
    if (S > 65535) S = 65535;                         // gridDim.y limit
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

    SupLJ lj;
#if LJ_COMB_RULE == LJ_COMB_FULL
    lj.cutforcesq = cuda_cutforcesq;
    lj.sigma6     = cuda_sigma6;
    lj.epsilon    = cuda_epsilon;
    lj.ntypes     = atom->ntypes;
#else
    lj.cutforcesq = param->cutforce * param->cutforce;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    lj.sigma6  = param->sigma6;
    lj.epsilon = param->epsilon;
#endif
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

    {
        // The positions moved since the last step: repack x/y/z/q (as the x_to_nbat
        // kernel of GROMACS), inside the timed region so FORCE includes it
        const int nclusters = numSupClusters(atom);
        const int nslots    = nclusters * CLUSTER_N;
        packXqSup<<<(nslots + 255) / 256, 256>>>(cuda_cl_x, cuda_cl_q, sup_xq, nclusters);
        cuda_assert("packXqSup", cudaPeekAtLastError());
    }

    grid_size.y = chooseJSlices(atom, neighbor);
    if (neighbor->half_neigh) {
        computeForceLJCudaSup_halfwarp<<<grid_size, block_size>>>(sup_xq,
            sup_par,
            cuda_cl_f,
            atom->Nclusters_local,
            cuda_numneigh_inner,
            sup_nbl,
            neighbor->maxneighs,
            cs,
            lj);
    } else {
        computeForceLJCudaSup_fullwarp<<<grid_size, block_size>>>(sup_xq,
            sup_par,
            cuda_cl_f,
            atom->Nclusters_local,
            cuda_numneigh_inner,
            sup_nbl,
            neighbor->maxneighs,
            cs,
            lj);
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

// One block per supercluster, CLUSTER_N x CLUSTER_M threads -- same shape as
// computeForceLJCudaSup_*() -- instead of one thread per supercluster
// serially checking all SCLUSTER_SIZE sub-clusters. Each thread owns one
// (cii, cjj) pair and loops over sub-clusters; imask is reduced via shared
// memory. Compaction stays single-threaded (tid 0): it's cheap bookkeeping,
// not the part this rewrite parallelizes.
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
    packSup(atom, neighbor, 0);
    DEBUG_MESSAGE("pruneNeighborCUDASup end\n");
}
