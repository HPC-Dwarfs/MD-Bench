/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
//---
#include <likwid-marker.h>

#include <device.h>
#include <gpu_profiler.h>

extern "C" {
#include <allocate.h>
#include <atom.h>
#include <comm.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <timing.h>
#include <util.h>
}

#include <energy_gpu.h>
#include <kernel_select_gpu.h>
#include <pair_scalar.h>

#ifdef __HIP_PLATFORM_AMD__
#define VL_WARP_SIZE 64
#define VL_LANE_MASK unsigned long long
// Lane reduction step of the force kernels on AMD, see "Reduction" above the kernels.
// DPP row_shl:n: lane l reads lane l+n of its 16-lane row (0 past the row end). It is a
// modifier of the add itself, much cheaper than __shfl_down (ds_bpermute). The builtin is
// int-only (a float argument would be value-converted!), hence the bit casts; a double
// is moved as two 32-bit halves.
#define VL_DPP_I32(v, n) __builtin_amdgcn_update_dpp(0, v, 0x100 + (n), 0xf, 0xf, false)
#if PRECISION == 1
#define VL_ROW_DOWN(x, n) __int_as_float(VL_DPP_I32(__float_as_int(x), n))
#else
#define VL_ROW_DOWN(x, n)                                                                \
    __hiloint2double(VL_DPP_I32(__double2hiint(x), n), VL_DPP_I32(__double2loint(x), n))
#endif

#define VL_ROW_STEP(x, y, z, n)                                                          \
    do {                                                                                 \
        x += VL_ROW_DOWN(x, n);                                                          \
        y += VL_ROW_DOWN(y, n);                                                          \
        z += VL_ROW_DOWN(z, n);                                                          \
    } while (0)
#else
#define VL_WARP_SIZE 32
#define VL_LANE_MASK unsigned int
#endif

// par[] is only read with geometric LJ or exclusions (the matrix rule reads mol directly)
template <int RULE, bool EXCL>
static __host__ __device__ constexpr bool vlPackPar()
{
    return RULE == LJ_COMB_GEOM || (EXCL && RULE != LJ_COMB_MATRIX);
}

static bool vlPackPar(int rule, bool excl)
{
    return rule == LJ_COMB_GEOM || (excl && rule != LJ_COMB_MATRIX);
}

#if PRECISION == 1
#define make_md_float4(x, y, z, w) make_float4((x), (y), (z), (w))
#define int_as_md_float(i)         __int_as_float(i)
#define md_float_as_int(f)         __float_as_int(f)
#else
#define make_md_float4(x, y, z, w) make_double4((x), (y), (z), (w))
#define int_as_md_float(i)         __longlong_as_double((long long)(i))
#define md_float_as_int(f)         ((int)__double_as_longlong(f))
#endif

// Read-only 16/32-byte load (there is no __ldg overload for double4)
static __device__ __forceinline__ MD_FLOAT4 ldg_md_float4(const MD_FLOAT4* p)
{
#if PRECISION == 1
    return __ldg(p);
#else
    const double2 lo = __ldg((const double2*)p);
    const double2 hi = __ldg((const double2*)p + 1);
    return make_double4(lo.x, lo.y, hi.x, hi.y);
#endif
}

// xq = {x, y, z, type (matrix rule) or charge}, par = {sqrt_epsilon, sigma3, mol, 0}
__global__ void packAtomsVL(DeviceAtom a,
    int n,
    MD_FLOAT4* __restrict__ xq,
    MD_FLOAT4* __restrict__ par,
    int pack_par,
    int rule,
    int excl)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) {
        return;
    }
    DeviceAtom* atom = &a;
#if COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT q = atom->q[i];
#else
    const MD_FLOAT q = (MD_FLOAT)0.0;
#endif
    const MD_FLOAT w = rule == LJ_COMB_MATRIX ? int_as_md_float(atom->type[i]) : q;
    xq[i]            = make_md_float4(atom_x(i), atom_y(i), atom_z(i), w);
    if (pack_par) {
        const bool geom   = rule == LJ_COMB_GEOM;
        const MD_FLOAT p0 = geom ? atom->sqrt_epsilon[i] : (MD_FLOAT)0.0;
        const MD_FLOAT p1 = geom ? atom->sigma3[i] : (MD_FLOAT)0.0;
        const MD_FLOAT p2 = excl ? int_as_md_float(atom->mol[i]) : (MD_FLOAT)0.0;
        par[i]            = make_md_float4(p0, p1, p2, (MD_FLOAT)0.0);
    }
}

// Arguments of the force kernels
typedef struct {
    DeviceAtom a;
    DeviceNeighbor neigh;
    int Nlocal;
    int ntypes;
    MD_FLOAT cutforcesq, sigma6, epsilon; // single rule (cutforcesq also geometric)
    CoulombScalar cs;
    const MD_FLOAT4* xq;
    const MD_FLOAT4* par;
    double* energies; // -energy sums, NULL without
} VLForceArgs;

/* Force kernels: T = 2^LGT lanes per atom (TPA, as in the LAMMPS GPU package), T = 16
 * with full and 8 with half lists. Lane t of atom i handles neighbors t, t+T, t+2T, ...;
 * the list is laid out so that these T lanes read consecutive entries (neighs()). The
 * host keeps the T lanes of an atom consecutive in one warp.
 *
 * Reduction: a tree over log2(T) steps into lane 0 (lane t += lane t+off). Only lane 0's
 * result is used, and it only reads lanes of its own atom.
 *  - NVIDIA: shuffle mask = this atom's T lanes only (atoms i >= Nlocal have returned)
 *  - AMD: DPP row shifts, which stay inside a 16-lane row, hence T <= 16
 *
 * HALF: half lists (Newton's third law for local j), RULE: LJ combination rule, ENERGY:
 * also sum the pair energies (-energy), EXCL: molecular exclusions, ACCURATE: IEEE
 * reciprocals (--accurate-recip)
 */
template <bool HALF, int RULE, bool ENERGY, bool EXCL, bool ACCURATE>
static __device__ __forceinline__ void forceLJVL(const VLForceArgs& args)
{
    constexpr int LGT  = HALF ? VL_TPA_LOG2_HALF : VL_TPA_LOG2_FULL;
    constexpr int T    = 1 << LGT;
    constexpr bool PAR = vlPackPar<RULE, EXCL>();
    const int gtid     = blockIdx.x * blockDim.x + threadIdx.x;
    const int i        = gtid >> LGT;
    const int t        = gtid & (T - 1);
    const int Nlocal   = args.Nlocal;
    const int ntypes   = args.ntypes;
    extern __shared__ MD_FLOAT4 lj_table[];
    if constexpr (RULE == LJ_COMB_MATRIX) {
        for (int k = threadIdx.x; k < ntypes * ntypes; k += blockDim.x) {
            lj_table[k] = make_md_float4(args.a.cutforcesq[k],
                args.a.sigma6[k],
                args.a.epsilon[k],
                (MD_FLOAT)0.0);
        }
        __syncthreads();
    }
    if (i >= Nlocal) {
        return;
    }
    // The host set the same value; as a constant, neighs() indexes with immediate shifts
    DeviceNeighbor neigh = args.neigh;
    neigh.tpa_log2       = LGT;

    const DeviceAtom* atom      = &args.a;
    const int* neigh_neighbors  = neigh.neighbors;
    const int numneighs         = neigh.numneigh[i];
    const MD_FLOAT4* __restrict__ xq  = args.xq;
    const MD_FLOAT4* __restrict__ par = args.par;

    const MD_FLOAT4 xqi = xq[i];
    MD_FLOAT fix        = 0;
    MD_FLOAT fiy        = 0;
    MD_FLOAT fiz        = 0;
    double vlj          = 0.0;
    double vc           = 0.0;

    MD_FLOAT4 pi = make_md_float4((MD_FLOAT)0.0, (MD_FLOAT)0.0, (MD_FLOAT)0.0, (MD_FLOAT)0.0);
    if constexpr (PAR) {
        pi = par[i];
    }
    // Row of i in lj_table (matrix rule)
    const int row_i = RULE == LJ_COMB_MATRIX ? md_float_as_int(xqi.w) * ntypes : 0;
#if COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT qi = args.cs.epsfac * (RULE == LJ_COMB_MATRIX ? atom->q[i] : xqi.w);
#else
    const MD_FLOAT qi = (MD_FLOAT)0.0;
#endif
    int mol_i = 0;
    if constexpr (EXCL) {
        mol_i = RULE == LJ_COMB_MATRIX ? atom->mol[i] : md_float_as_int(pi.z);
    }

    auto pair = [&](int k) {
        const int j         = __ldg(&neighs(neigh_neighbors, i, k, Nlocal, &neigh));
        const MD_FLOAT4 xqj = ldg_md_float4(&xq[j]);
        MD_FLOAT delx       = xqi.x - xqj.x;
        MD_FLOAT dely       = xqi.y - xqj.y;
        MD_FLOAT delz       = xqi.z - xqj.z;
        MD_FLOAT rsq        = delx * delx + dely * dely + delz * delz;

        MD_FLOAT4 pj = make_md_float4((MD_FLOAT)0.0, (MD_FLOAT)0.0, (MD_FLOAT)0.0, (MD_FLOAT)0.0);
        if constexpr (PAR) {
            pj = ldg_md_float4(&par[j]);
        }
        MD_FLOAT cutforcesq = args.cutforcesq;
        MD_FLOAT sigma6     = args.sigma6;
        MD_FLOAT epsilon    = args.epsilon;
        if constexpr (RULE == LJ_COMB_GEOM) {
            sigma6  = pi.y * pj.y;
            epsilon = pi.x * pj.x;
        } else if constexpr (RULE == LJ_COMB_MATRIX) {
            const MD_FLOAT4 lj = lj_table[row_i + md_float_as_int(xqj.w)];
            cutforcesq         = lj.x;
            sigma6             = lj.y;
            epsilon            = lj.z;
        }

        if (rsq < cutforcesq) {
#if COULOMB_RULE != COULOMB_NONE
            const MD_FLOAT qj = RULE == LJ_COMB_MATRIX ? __ldg(&atom->q[j]) : xqj.w;
#else
            const MD_FLOAT qj = (MD_FLOAT)0.0;
#endif
            int mol_j = 0;
            if constexpr (EXCL) {
                mol_j = RULE == LJ_COMB_MATRIX ? __ldg(&atom->mol[j])
                                               : md_float_as_int(pj.z);
            }
            MD_FLOAT force = pairForceScalarAcc(
                rsq, sigma6, epsilon, qi, qj, mol_i, mol_j, EXCL, ACCURATE, args.cs);
            if constexpr (ENERGY) {
                MD_FLOAT e_lj, e_c;
                pairEnergyScalarAcc(rsq,
                    sigma6,
                    epsilon,
                    qi,
                    qj,
                    mol_i,
                    mol_j,
                    EXCL,
                    ACCURATE,
                    args.cs,
                    &e_lj,
                    &e_c);
                // Half lists: ghost pairs are computed from both sides (no ghost force
                // fold-back)
                const double w = (!HALF || j < Nlocal) ? 1.0 : 0.5;
                vlj += w * e_lj;
                vc += w * e_c;
            }
            MD_FLOAT partial_force_x = delx * force;
            MD_FLOAT partial_force_y = dely * force;
            MD_FLOAT partial_force_z = delz * force;

            if constexpr (HALF) {
                if (j < Nlocal) {
                    atomicAdd(&atom_fx(j), -partial_force_x);
                    atomicAdd(&atom_fy(j), -partial_force_y);
                    atomicAdd(&atom_fz(j), -partial_force_z);
                }
            }

            fix += partial_force_x;
            fiy += partial_force_y;
            fiz += partial_force_z;
        }
    };

    // ceil(numneighs / T) chunks of T neighbors
    const int nchunks = (numneighs + T - 1) >> LGT;
    if (gridDim.y == 1) {
        for (int c = 0; c < nchunks; c++) {
            const int k = (c << LGT) + t;
            if (T == 1 || k < numneighs) {
                pair(k);
            }
        }
    } else {
#pragma unroll 2
        for (int c = blockIdx.y; c < nchunks; c += gridDim.y) {
            const int k = (c << LGT) + t;
            if (T == 1 || k < numneighs) {
                pair(k);
            }
        }
    }

    // Reduction of the T partial forces into lane 0, see above
    if (T > 1) {
#ifdef __HIP_PLATFORM_AMD__
        static_assert(T <= 16, "the lanes of an atom must lie in one DPP row");
        if (T > 8) VL_ROW_STEP(fix, fiy, fiz, 8);
        if (T > 4) VL_ROW_STEP(fix, fiy, fiz, 4);
        if (T > 2) VL_ROW_STEP(fix, fiy, fiz, 2);
        VL_ROW_STEP(fix, fiy, fiz, 1);
#else
        // mask = T bits at the atom's first lane (T % VL_WARP_SIZE: no shift by the type
        // width when T == VL_WARP_SIZE, where that branch is not taken)
        const unsigned int lane  = threadIdx.x & (VL_WARP_SIZE - 1);
        const VL_LANE_MASK group = T >= VL_WARP_SIZE
                                       ? ~(VL_LANE_MASK)0
                                       : (((VL_LANE_MASK)1 << (T % VL_WARP_SIZE)) - 1);
        const VL_LANE_MASK mask  = group << (lane & ~(unsigned int)(T - 1));
#pragma unroll
        for (int off = T / 2; off > 0; off /= 2) {
            fix += __shfl_down_sync(mask, fix, off, T);
            fiy += __shfl_down_sync(mask, fiy, off, T);
            fiz += __shfl_down_sync(mask, fiz, off, T);
        }
#endif
    }
    if (t == 0) {
        if (!HALF && gridDim.y == 1) {
            atom_fx(i) = fix;
            atom_fy(i) = fiy;
            atom_fz(i) = fiz;
        } else {
            atomicAdd(&atom_fx(i), fix);
            atomicAdd(&atom_fy(i), fiy);
            atomicAdd(&atom_fz(i), fiz);
        }
    }
    if constexpr (ENERGY) {
        warpEnergiesAdd(vlj, vc, args.energies);
    }
}

template <int RULE, bool ENERGY, bool EXCL, bool ACCURATE>
__global__ void computeForceLJCudaFullNeigh(const VLForceArgs args)
{
    forceLJVL<false, RULE, ENERGY, EXCL, ACCURATE>(args);
}

template <int RULE, bool ENERGY, bool EXCL, bool ACCURATE>
__global__ void computeForceLJCudaHalfNeigh(const VLForceArgs args)
{
    forceLJVL<true, RULE, ENERGY, EXCL, ACCURATE>(args);
}

__global__ void kernel_initial_integrate(
    MD_FLOAT dtforce, MD_FLOAT dt, int Nlocal, DeviceAtom a)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom = &a;

    atom_vx(i) += dtforce * atom_fx(i);
    atom_vy(i) += dtforce * atom_fy(i);
    atom_vz(i) += dtforce * atom_fz(i);
    atom_x(i) = atom_x(i) + dt * atom_vx(i);
    atom_y(i) = atom_y(i) + dt * atom_vy(i);
    atom_z(i) = atom_z(i) + dt * atom_vz(i);
}

__global__ void kernel_final_integrate(MD_FLOAT dtforce, int Nlocal, DeviceAtom a)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom = &a;

    atom_vx(i) += dtforce * atom_fx(i);
    atom_vy(i) += dtforce * atom_fy(i);
    atom_vz(i) += dtforce * atom_fz(i);
}

extern "C" {

void finalIntegrateCUDA(bool reneigh, Parameter* param, Atom* atom)
{
    DEBUG_MESSAGE("finalIntegrateCUDA begin\n");

    const int Nlocal                = atom->Nlocal;
    const int num_threads_per_block = get_cuda_num_threads();
    const int num_blocks            = ceil((float)Nlocal / (float)num_threads_per_block);

    kernel_final_integrate<<<num_blocks, num_threads_per_block>>>(param->dtforce,
        Nlocal,
        atom->d_atom);
    cuda_assert("kernel_final_integrate", cudaPeekAtLastError());
    cuda_assert("kernel_final_integrate", cudaDeviceSynchronize());

    if (reneigh) {
        memcpyFromGPU(atom->vx, atom->d_atom.vx, sizeof(MD_FLOAT) * atom->Nlocal * 3);
    }

    DEBUG_MESSAGE("finalIntegrateCUDA end\n");
}

void initialIntegrateCUDA(bool reneigh, Parameter* param, Atom* atom)
{
    DEBUG_MESSAGE("initialIntegrateCUDA begin\n");

    const int Nlocal                = atom->Nlocal;
    const int num_threads_per_block = get_cuda_num_threads();
    const int num_blocks            = ceil((float)Nlocal / (float)num_threads_per_block);

    kernel_initial_integrate<<<num_blocks, num_threads_per_block>>>(param->dtforce,
        param->dt,
        Nlocal,
        atom->d_atom);
    cuda_assert("kernel_initial_integrate", cudaPeekAtLastError());
    cuda_assert("kernel_initial_integrate", cudaDeviceSynchronize());

    if (reneigh) {
        memcpyFromGPU(atom->vx, atom->d_atom.vx, sizeof(MD_FLOAT) * atom->Nlocal * 3);
    }

    DEBUG_MESSAGE("initialIntegrateCUDA end\n");
}

// Neighbor entries of all local atoms after the last list build (neighborCuda.cu)
extern "C" long numNeighborsTotalCUDA(void);
extern "C" int neighborBuildCountCUDA(void);

/* Packs xq every call and par after each list build (atoms are only reordered, exchanged and
 * given new ghosts together with a build), for all local and ghost atoms. */
static void packAtomsForForceVL(
    Parameter* param, Atom* atom, MD_FLOAT4** xq, MD_FLOAT4** par)
{
    static MD_FLOAT4* d_xq  = NULL;
    static MD_FLOAT4* d_par = NULL;
    static int nalloc       = 0;
    static int packed_build = -1;
    const int rule          = param->lj_comb_rule;
    const bool excl         = atom->molsize > 1;
    const bool use_par      = vlPackPar(rule, excl);
    const int n             = atom->Nlocal + atom->Nghost;
    int pack_par            = use_par && neighborBuildCountCUDA() != packed_build;
    if (n > nalloc) {
        nalloc = atom->Nmax > n ? atom->Nmax : n;
        d_xq   = (MD_FLOAT4*)reallocateGPU(d_xq, sizeof(MD_FLOAT4) * nalloc);
        if (use_par) {
            d_par    = (MD_FLOAT4*)reallocateGPU(d_par, sizeof(MD_FLOAT4) * nalloc);
            pack_par = 1;
        }
    }
    if (n > 0) {
        const int threads = 256;
        packAtomsVL<<<(n + threads - 1) / threads, threads>>>(atom->d_atom,
            n,
            d_xq,
            d_par,
            pack_par,
            rule,
            excl);
        cuda_assert("packAtomsVL", cudaPeekAtLastError());
    }
    packed_build = neighborBuildCountCUDA();
    *xq          = d_xq;
    *par         = d_par;
}

// Resident blocks per SM of a force kernel: one wave is N_SM of these. Queried once per kernel.
static int residentBlocksPerSMVL(const void* kernel, int threads, size_t shmem)
{
    static const void* kernels[16];
    static int blocks[16];
    static int nkernels = 0;
    for (int i = 0; i < nkernels; i++) {
        if (kernels[i] == kernel) return blocks[i];
    }
    int n = 0;
    cuda_assert("cudaOccupancyMaxActiveBlocksPerMultiprocessor",
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n, kernel, threads, shmem));
    if (n < 1) n = 1;
    if (nkernels < 16) {
        kernels[nkernels] = kernel;
        blocks[nkernels]  = n;
        nkernels++;
    }
    return n;
}

// j-slices (grid.y) of this force call, T threads per atom and slice, see --gpu-j-slices
static int chooseJSlicesVL(const Parameter* param,
    const Neighbor* neighbor,
    const void* kernel,
    int nblocks,
    int threads,
    size_t shmem)
{
    if (param->gpu_j_slices > 0) {
        return param->gpu_j_slices;
    }
    // S = K * sqrt(N_SM * n_tot) / N_b if K > 0, else W waves of resident blocks;
    // K and W per list layout and type, tuned on GH200 (CSR uses the AOS values)
#if defined(NBLIST_AOS) || defined(NBLIST_CSR)
    const double k_half = 0.0, w_half = 0.35, k_full = 0.0, w_full = 0.23;
#else
    const double k_half = 0.13, w_half = 0.0, k_full = 0.0, w_full = 0.3;
#endif
    static int nsm = 0;
    if (nsm == 0) {
        int dev = 0;
        cuda_assert("cudaGetDevice", cudaGetDevice(&dev));
        cuda_assert("cudaDeviceGetAttribute",
            cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev));
    }
    if (nblocks <= 0) return 1;
    const int bps       = residentBlocksPerSMVL(kernel, threads, shmem);
    const double k      = neighbor->half_neigh ? k_half : k_full;
    const double w      = neighbor->half_neigh ? w_half : w_full;
    const double s      = k > 0.0 ? k * sqrt((double)nsm) * sqrt((double)numNeighborsTotalCUDA()) /
                                   (double)nblocks
                                  : w * nsm * bps / (double)nblocks;
    int S               = (int)lround(s);
    const long max_blocks = 10L * nsm * bps; // 10 waves
    if (S > max_blocks / nblocks) S = (int)(max_blocks / nblocks);
    if (S < 1) S = 1;
    if (S > neighbor->maxneighs) S = neighbor->maxneighs > 0 ? neighbor->maxneighs : 1;
    if (S > 65535) S = 65535;                         // gridDim.y limit
    return S;
}

double computeForceLJCUDA(Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCUDA begin\n");

    const int num_threads_per_block = get_cuda_num_threads();
    int Nlocal                      = atom->Nlocal;
    int Nmax                        = atom->Nmax;

    VLForceArgs args;
    args.a          = atom->d_atom;
    args.neigh      = neighbor->d_neighbor;
    args.Nlocal     = Nlocal;
    args.ntypes     = atom->ntypes;
    args.cutforcesq = param->cutforce * param->cutforce;
    args.sigma6     = param->sigma6;
    args.epsilon    = param->epsilon;
    args.cs         = coulombScalarInit(param);
    args.energies   = deviceEnergies(param);

    const bool excl     = atom->molsize > 1; // molecular exclusions
    const int rule      = param->lj_comb_rule;
    const bool energy   = param->compute_energy;
    const bool accurate = param->accurate_reciprocal;

    // HINT: Run with cuda-memcheck ./MDBench-NVCC in case of error

    GPU_PROFILE_START("force_lj");
    double S = getTimeStamp();
    LIKWID_MARKER_START("force");

    MD_FLOAT4 *xq, *par;
    packAtomsForForceVL(param, atom, &xq, &par);
    args.xq  = xq;
    args.par = par;

    // Lanes per atom T, see the comment above forceLJVL
    const int half  = neighbor->half_neigh;
    const int lanes = 1 << (half ? VL_TPA_LOG2_HALF : VL_TPA_LOG2_FULL);
    auto kernel     = half
                          ? SELECT_GPU_KERNEL(
                            computeForceLJCudaHalfNeigh, rule, energy, excl, accurate)
                          : SELECT_GPU_KERNEL(
                            computeForceLJCudaFullNeigh, rule, energy, excl, accurate);
    // the lanes of an atom must stay in one block (and so in one warp) for the reduction
    if (num_threads_per_block % lanes != 0) {
        fprintf(stderr, "Error: NUM_THREADS must be a multiple of %d!\n", lanes);
        exit(EXIT_FAILURE);
    }
    const size_t shmem = rule == LJ_COMB_MATRIX // the type-pair table
                             ? (size_t)atom->ntypes * atom->ntypes * sizeof(MD_FLOAT4)
                             : 0;
    if (shmem > 48 * 1024) {
        fprintf(stderr,
            "Error: the table of %d atom types does not fit in shared memory!\n",
            atom->ntypes);
        exit(EXIT_FAILURE);
    }
    // grid.x: T threads per local atom; grid.y: j-slices
    const int bx = (int)(((long)Nlocal * lanes + num_threads_per_block - 1) /
                         num_threads_per_block);
    const dim3 grid(bx,
        chooseJSlicesVL(param,
            neighbor,
            (const void*)kernel,
            bx,
            num_threads_per_block,
            shmem),
        1);
    if (half || grid.y > 1) { // the kernel accumulates with atomicAdd
#ifdef ATOM_POSITION_AOS
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax * 3);
#else
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fy, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fz, 0, sizeof(MD_FLOAT) * Nmax);
#endif
    }
    if (bx > 0) {
        kernel<<<grid, num_threads_per_block, shmem>>>(args);
    }

    cuda_assert("computeForceLJCuda", cudaPeekAtLastError());
    cuda_assert("computeForceLJCuda", cudaDeviceSynchronize());
    GPU_PROFILE_STOP();

    LIKWID_MARKER_STOP("force");
    double E = getTimeStamp();
    // A full list has every pair twice (the half list kernel weights its pairs itself)
    storeDeviceEnergies(param,
        Nlocal,
        atom->type,
        atom->q_per_type,
        args.energies,
        half ? 1.0 : 0.5);
    DEBUG_MESSAGE("computeForceLJCUDA end\n");
    return E - S;
}
}

// Reference force kernel (-ref): one thread per atom, no TPA, packing or slicing
template <int RULE, bool ENERGY, bool EXCL, bool ACCURATE>
__global__ void computeForceLJCudaRef_warp(DeviceAtom a,
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    int Nlocal,
    DeviceNeighbor neigh,
    int* neigh_neighbors,
    int* neigh_numneigh,
    int ntypes,
    int half,
    CoulombScalar cs,
    double* energies)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom    = &a;
    const int numneighs = neigh_numneigh[i];
    MD_FLOAT xtmp       = atom_x(i);
    MD_FLOAT ytmp       = atom_y(i);
    MD_FLOAT ztmp       = atom_z(i);
    MD_FLOAT3 fi        = { 0.0, 0.0, 0.0 };
    double vlj          = 0.0;
    double vc           = 0.0;
#if COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT qi = cs.epsfac * atom->q[i];
#else
    const MD_FLOAT qi = (MD_FLOAT)0.0;
#endif
    const int mol_i = EXCL ? atom->mol[i] : 0;

    for (int k = 0; k < numneighs; k++) {
        // neighs() hides the list layout (AOS, SOA with lanes per atom)
        int j         = neighs(neigh_neighbors, i, k, Nlocal, &neigh);
        MD_FLOAT delx = xtmp - atom_x(j);
        MD_FLOAT dely = ytmp - atom_y(j);
        MD_FLOAT delz = ztmp - atom_z(j);
        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

        MD_FLOAT cutforcesq_ij = cutforcesq;
        MD_FLOAT sigma6_ij     = sigma6;
        MD_FLOAT epsilon_ij    = epsilon;
        if constexpr (RULE == LJ_COMB_GEOM) {
            sigma6_ij  = atom->sigma3[i] * atom->sigma3[j];
            epsilon_ij = atom->sqrt_epsilon[i] * atom->sqrt_epsilon[j];
        } else if constexpr (RULE == LJ_COMB_MATRIX) {
            int type_ij    = atom->type[i] * ntypes + atom->type[j];
            cutforcesq_ij  = atom->cutforcesq[type_ij];
            sigma6_ij      = atom->sigma6[type_ij];
            epsilon_ij     = atom->epsilon[type_ij];
        }

        if (rsq < cutforcesq_ij) {
#if COULOMB_RULE != COULOMB_NONE
            const MD_FLOAT qj = atom->q[j];
#else
            const MD_FLOAT qj = (MD_FLOAT)0.0;
#endif
            const int mol_j = EXCL ? atom->mol[j] : 0;
            MD_FLOAT force  = pairForceScalarAcc(rsq,
                sigma6_ij,
                epsilon_ij,
                qi,
                qj,
                mol_i,
                mol_j,
                EXCL,
                ACCURATE,
                cs);
            if constexpr (ENERGY) {
                MD_FLOAT e_lj, e_c;
                pairEnergyScalarAcc(rsq,
                    sigma6_ij,
                    epsilon_ij,
                    qi,
                    qj,
                    mol_i,
                    mol_j,
                    EXCL,
                    ACCURATE,
                    cs,
                    &e_lj,
                    &e_c);
                // Half lists: ghost pairs are computed from both sides
                const double w = (!half || j < Nlocal) ? 1.0 : 0.5;
                vlj += w * e_lj;
                vc += w * e_c;
            }
            MD_FLOAT fx = delx * force;
            MD_FLOAT fy = dely * force;
            MD_FLOAT fz = delz * force;

            fi.x += fx;
            fi.y += fy;
            fi.z += fz;

            if (half && j < Nlocal) {
                atomicAdd(&atom_fx(j), -fx);
                atomicAdd(&atom_fy(j), -fy);
                atomicAdd(&atom_fz(j), -fz);
            }
        }
    }

    atomicAdd(&atom_fx(i), fi.x);
    atomicAdd(&atom_fy(i), fi.y);
    atomicAdd(&atom_fz(i), fi.z);
    if constexpr (ENERGY) {
        warpEnergiesAdd(vlj, vc, energies);
    }
}

extern "C" double computeForceLJCUDARef(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCUDARef begin\n");

    const int num_threads_per_block = get_cuda_num_threads();
    int Nlocal                      = atom->Nlocal;
    int Nmax                        = atom->Nmax;
    const CoulombScalar cs          = coulombScalarInit(param);
    double* energies                = deviceEnergies(param);
    const bool excl                 = atom->molsize > 1; // molecular exclusions
    const int rule                  = param->lj_comb_rule;
    const bool energy               = param->compute_energy;
    const bool accurate             = param->accurate_reciprocal;
    auto kernel                     = SELECT_GPU_KERNEL(
        computeForceLJCudaRef_warp, rule, energy, excl, accurate);

    GPU_PROFILE_START("force_lj");
    const int num_blocks = (Nlocal + num_threads_per_block - 1) / num_threads_per_block;
    double S             = getTimeStamp();
    LIKWID_MARKER_START("force");

#ifdef ATOM_POSITION_AOS
    memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax * 3);
#else
    memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax);
    memsetGPU(atom->d_atom.fy, 0, sizeof(MD_FLOAT) * Nmax);
    memsetGPU(atom->d_atom.fz, 0, sizeof(MD_FLOAT) * Nmax);
#endif
    if (num_blocks > 0) {
        kernel<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
            param->cutforce * param->cutforce,
            param->sigma6,
            param->epsilon,
            Nlocal,
            neighbor->d_neighbor,
            neighbor->d_neighbor.neighbors,
            neighbor->d_neighbor.numneigh,
            atom->ntypes,
            neighbor->half_neigh,
            cs,
            energies);
    }

    cuda_assert("computeForceLJCudaRef", cudaPeekAtLastError());
    cuda_assert("computeForceLJCudaRef", cudaDeviceSynchronize());
    GPU_PROFILE_STOP();

    LIKWID_MARKER_STOP("force");
    double E = getTimeStamp();
    storeDeviceEnergies(param,
        Nlocal,
        atom->type,
        atom->q_per_type,
        energies,
        neighbor->half_neigh ? 1.0 : 0.5);
    DEBUG_MESSAGE("computeForceLJCUDARef end\n");
    return E - S;
}

extern "C" void copyGhostFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->x, atom->d_atom.x, atom->Nlocal * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyGhostToGPU(Atom* atom)
{
    memcpyToGPU(&atom->d_atom.x[atom->Nlocal * 3],
        &atom->x[atom->Nlocal * 3],
        atom->Nghost * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyForceFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->fx, atom->d_atom.fx, atom->Nmax * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyForceToGPU(Atom* atom)
{
    memcpyToGPU(atom->d_atom.fx, atom->fx, atom->Nmax * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyDataFromCUDADevice(Parameter* param, Atom* atom)
{
    memcpyFromGPU(atom->x, atom->d_atom.x, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyFromGPU(atom->vx, atom->d_atom.vx, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyFromGPU(atom->type, atom->d_atom.type, atom->Nmax * sizeof(int));
}

extern "C" void copyDataToCUDADevice(Parameter* param, Atom* atom)
{
    memcpyToGPU(atom->d_atom.x, atom->x, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyToGPU(atom->d_atom.vx, atom->vx, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyToGPU(atom->d_atom.type, atom->type, atom->Nmax * sizeof(int));
}

/*
__global__ void pack_forward_Cuda(DeviceAtom a,
    int n,
    int* cuda_list,
    MD_FLOAT* buf,
    int PBCx,
    int PBCy,
    int PBCz,
    MD_FLOAT xprd,
    MD_FLOAT yprd,
    MD_FLOAT zprd)
{
    DeviceAtom* atom    = &a;
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    int j        = cuda_list[i];
    buf_x(i) = atom_x(j) + PBCx * xprd;
    buf_y(i) = atom_y(j) + PBCy * yprd;
    buf_z(i) = atom_z(j) + PBCz * zprd;
}

__global__ void unpack_forward_Cuda(DeviceAtom a,
    int n,
    int first,
    MD_FLOAT* buf)
{
    DeviceAtom* atom    = &a;
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    atom_x((first + i)) = buf_x(i);
    atom_y((first + i)) = buf_y(i);
    atom_z((first + i)) = buf_z(i);
}

extern "C" void forwardCommCUDA(Comm* comm, Atom* atom, int iswap)
{
    int nrqst = 0, offset = 0, nsend = 0, nrecv = 0;
    int pbc[3];
    int size    = comm->forwardSize;
    int maxrqst = comm->numneigh;
    int* cuda_sendlist;
    int max_list_size;
    cuda_buf_send =  (MD_FLOAT*)allocateGPU(comm->maxsend * sizeof(MD_FLOAT));
    cuda_buf_recv =  (MD_FLOAT*)allocateGPU(comm->maxrecv * sizeof(MD_FLOAT));

    //use a single buffer and takes the highes list size to move list of cluster to send
    for (int ineigh = 0; ineigh < comm->numneigh; ineigh++){
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

        const int num_threads_per_block = get_cuda_num_threads();
        const int num_blocks = ceil((float) comm->atom_send[ineigh] /
(float)num_threads_per_block);

        pack_forward_Cuda<<<num_blocks, num_threads_per_block>>>(
                                            atom->d_atom,                   //MD_FLOAT*
-->need to be in tye device comm->atom_send[ineigh],        //int cuda_sendlist, //int*
-->need to be in tye device &cuda_buf_send[offset*size],    //MD_FLOAT*  -->need to be in
tye device pbc[_x],                        //int pbc[_y],                        //int
                                            pbc[_z],                        //int
                                            atom->mybox.xprd, //MD_FLOAT atom->mybox.yprd,
//MD_FLOAT atom->mybox.zprd);                    //MD_FLOAT cudaDeviceSynchronize();
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

    // unpack buffer
    for (int ineigh = comm->recvfrom[iswap]; ineigh < comm->recvtill[iswap]; ineigh++) {
        offset = comm->off_atom_recv[ineigh];
        MD_FLOAT *buf = (comm->othersend[iswap]) ? cuda_buf_recv : cuda_buf_send;

        const int num_threads_per_block = get_cuda_num_threads();
        const int num_blocks = ceil((float) comm->atom_send[ineigh] /
(float)num_threads_per_block);

        unpack_forward_Cuda<<<num_blocks, num_threads_per_block>>>(
                                        atom->d_atom,                       //MD_FLOAT*
--> need to be in the device comm->atom_recv[ineigh],            //int
                                        comm->firstrecv[iswap] + offset,    //int
                                        &buf[offset * size]);               //MD_FLOAT*
--> need to be in the devic cudaDeviceSynchronize();
        }
    cuda_assert("cudaDeviceFree", cudaFree(cuda_sendlist));
    cuda_assert("cudaDeviceFree", cudaFree(cuda_buf_recv));
    cuda_assert("cudaDeviceFree", cudaFree(cuda_buf_send));
}
*/
