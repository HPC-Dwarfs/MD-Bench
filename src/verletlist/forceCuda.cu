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

#include <pair_scalar.h>

/* Packed per-atom data of the force kernels, after LAMMPS' GPU package (x and type in one
 * float4): a neighbor costs one 16-byte load per array instead of up to seven scalar gathers
 * from five arrays (x, y, z, sqrt_epsilon, sigma3, q, mol).
 *   xq[i]  = {x, y, z, q}                                         every force call
 *   par[i] = {sqrt_epsilon, sigma3, mol, 0} (geometric) or
 *            {type, 0, mol, 0} (full matrix)                       after each list build
 * Integers (type, mol) are stored bitwise. With LJ_COMB_SINGLE and no exclusions the kernels
 * need no per-atom parameters and par is not used. */
#if LJ_COMB_RULE == LJ_COMB_GEOM || LJ_COMB_RULE == LJ_COMB_FULL || USE_EXCLUSIONS
#define VL_PACK_PAR 1
#else
#define VL_PACK_PAR 0
#endif

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

__global__ void packAtomsVL(DeviceAtom a,
    int n,
    MD_FLOAT4* __restrict__ xq,
    MD_FLOAT4* __restrict__ par,
    int pack_par)
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
    xq[i] = make_md_float4(atom_x(i), atom_y(i), atom_z(i), q);
#if VL_PACK_PAR
    if (pack_par) {
#if LJ_COMB_RULE == LJ_COMB_GEOM
        const MD_FLOAT p0 = atom->sqrt_epsilon[i];
        const MD_FLOAT p1 = atom->sigma3[i];
#elif LJ_COMB_RULE == LJ_COMB_FULL
        const MD_FLOAT p0 = int_as_md_float(atom->type[i]);
        const MD_FLOAT p1 = (MD_FLOAT)0.0;
#else
        const MD_FLOAT p0 = (MD_FLOAT)0.0;
        const MD_FLOAT p1 = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS
        const MD_FLOAT p2 = int_as_md_float(atom->mol[i]);
#else
        const MD_FLOAT p2 = (MD_FLOAT)0.0;
#endif
        par[i] = make_md_float4(p0, p1, p2, (MD_FLOAT)0.0);
    }
#else
    (void)par;
    (void)pack_par;
#endif
}

/* Full-list kernel with T = 2^VL_TPA_LOG2_FULL lanes per atom (TPA, as the half-list kernel
 * below): lane t of atom i handles the neighbors t, t+T, t+2T, ... (chunks of T, interleaved
 * in SOA lists, see neighs), the partial forces of i are summed over the T lanes with warp
 * shuffles and written once: no atomics without j-slices. */
__global__ void computeForceLJCudaFullNeigh(DeviceAtom a,
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    int Nlocal,
    DeviceNeighbor neigh,
    int* neigh_neighbors,
    int* neigh_numneigh,
    int ntypes,
    CoulombScalar cs,
    const MD_FLOAT4* __restrict__ xq,
    const MD_FLOAT4* __restrict__ par)
{

    constexpr int LGT = VL_TPA_LOG2_FULL; // T = 2^LGT lanes per atom
    constexpr int T   = 1 << LGT;
    const int gtid    = blockIdx.x * blockDim.x + threadIdx.x;
    const int i       = gtid >> LGT;
    const int t       = gtid & (T - 1);
    if (i >= Nlocal) {
        return; // all lanes of an atom leave together
    }
    neigh.tpa_log2 = LGT; // the list layout as a constant: index arithmetic with immediates

    DeviceAtom* atom         = &a;
    const int numneighs      = neigh_numneigh[i];

    const MD_FLOAT4 xqi = xq[i];
    MD_FLOAT fix = 0;
    MD_FLOAT fiy = 0;
    MD_FLOAT fiz = 0;

#if VL_PACK_PAR
    const MD_FLOAT4 pi = par[i];
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
    const MD_FLOAT sqrt_eps_i = pi.x;
    const MD_FLOAT sigma3_i   = pi.y;
#elif LJ_COMB_RULE == LJ_COMB_FULL
    const int type_i = md_float_as_int(pi.x);
#endif
#if COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT qi = cs.epsfac * xqi.w;
#else
    const MD_FLOAT qi = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS
    const int mol_i = md_float_as_int(pi.z);
#else
    const int mol_i = 0;
#endif

    // One neighbor of atom i
    auto pair = [&](int k) {
        // __ldg: neighbor indices and packed atom data are read-only here, but the
        // compiler cannot prove they don't alias fx, so it won't use the read-only
        // cache on its own
        const int j        = __ldg(&neighs(neigh_neighbors, i, k, Nlocal, &neigh));
        const MD_FLOAT4 xqj = ldg_md_float4(&xq[j]);
        MD_FLOAT delx       = xqi.x - xqj.x;
        MD_FLOAT dely       = xqi.y - xqj.y;
        MD_FLOAT delz       = xqi.z - xqj.z;
        MD_FLOAT rsq        = delx * delx + dely * dely + delz * delz;

#if VL_PACK_PAR
        const MD_FLOAT4 pj = ldg_md_float4(&par[j]);
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
        const MD_FLOAT sigma6  = sigma3_i * pj.y;
        const MD_FLOAT epsilon = sqrt_eps_i * pj.x;
#elif LJ_COMB_RULE == LJ_COMB_FULL
        const int type_ij         = type_i * ntypes + md_float_as_int(pj.x);
        const MD_FLOAT cutforcesq = __ldg(&atom->cutforcesq[type_ij]);
        const MD_FLOAT sigma6     = __ldg(&atom->sigma6[type_ij]);
        const MD_FLOAT epsilon    = __ldg(&atom->epsilon[type_ij]);
#endif

        if (rsq < cutforcesq) {
#if COULOMB_RULE != COULOMB_NONE
            const MD_FLOAT qj = xqj.w;
#else
            const MD_FLOAT qj = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS
            const int mol_j = md_float_as_int(pj.z);
#else
            const int mol_j = 0;
#endif
            MD_FLOAT force = pairForceScalar(rsq, sigma6, epsilon, qi, qj, mol_i, mol_j, cs);
            fix += delx * force;
            fiy += dely * force;
            fiz += delz * force;
        }
    };

    // Chunks of T neighbors; j-slices: block (bx, s) handles the chunks s, s+S, s+2S, ... of
    // its atoms (S = gridDim.y). Without slices the loop keeps unit stride, which nvcc
    // unrolls; the strided loop is not unrolled by itself and up to 46% slower at S = 1
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

    if (T > 1) {
        // the T lanes of an atom are consecutive lanes of one warp
        const unsigned int lane  = threadIdx.x & 31;
        const unsigned int group = T >= 32 ? 0xffffffffu : ((1u << (T % 32)) - 1u);
        const unsigned int mask  = group << (lane & ~(unsigned int)(T - 1));
#pragma unroll
        for (int off = T / 2; off > 0; off /= 2) {
            fix += __shfl_down_sync(mask, fix, off, T);
            fiy += __shfl_down_sync(mask, fiy, off, T);
            fiz += __shfl_down_sync(mask, fiz, off, T);
        }
    }
    if (t == 0) {
        if (gridDim.y == 1) {
            atom_fx(i) = fix;
            atom_fy(i) = fiy;
            atom_fz(i) = fiz;
        } else {
            // several j-slices add to the same atom
            atomicAdd(&atom_fx(i), fix);
            atomicAdd(&atom_fy(i), fiy);
            atomicAdd(&atom_fz(i), fiz);
        }
    }
}

/* Half-list kernel with T = 2^VL_TPA_LOG2_HALF lanes per atom (TPA, after the LAMMPS GPU
 * package and KOKKOS neigh/thread): lane t of atom i handles the neighbors t, t+T, t+2T, ...
 * (chunks of T, interleaved in SOA lists, see neighs), the partial forces of i are summed over
 * the T lanes with warp shuffles and added once. The reaction on a neighbor j is added
 * atomically, but only for local j: the GPU half list keeps all ghosts (pairs across the
 * boundary are computed by both sides) and ghost forces are not folded back (force.c), as
 * LAMMPS' half list with newton off. */
__global__ void computeForceLJCudaHalfNeigh(DeviceAtom a,
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    int Nlocal,
    DeviceNeighbor neigh,
    int* neigh_neighbors,
    int* neigh_numneigh,
    int ntypes,
    CoulombScalar cs,
    const MD_FLOAT4* __restrict__ xq,
    const MD_FLOAT4* __restrict__ par)
{

    constexpr int LGT = VL_TPA_LOG2_HALF; // T = 2^LGT lanes per atom
    constexpr int T   = 1 << LGT;
    const int gtid    = blockIdx.x * blockDim.x + threadIdx.x;
    const int i       = gtid >> LGT;
    const int t       = gtid & (T - 1);
    if (i >= Nlocal) {
        return; // all lanes of an atom leave together
    }
    neigh.tpa_log2 = LGT; // the list layout as a constant: index arithmetic with immediates

    DeviceAtom* atom         = &a;
    const int numneighs      = neigh_numneigh[i];

    const MD_FLOAT4 xqi = xq[i];
    MD_FLOAT fix = 0;
    MD_FLOAT fiy = 0;
    MD_FLOAT fiz = 0;

#if VL_PACK_PAR
    const MD_FLOAT4 pi = par[i];
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
    const MD_FLOAT sqrt_eps_i = pi.x;
    const MD_FLOAT sigma3_i   = pi.y;
#elif LJ_COMB_RULE == LJ_COMB_FULL
    const int type_i = md_float_as_int(pi.x);
#endif
#if COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT qi = cs.epsfac * xqi.w;
#else
    const MD_FLOAT qi = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS
    const int mol_i = md_float_as_int(pi.z);
#else
    const int mol_i = 0;
#endif

    // One neighbor of atom i
    auto pair = [&](int k) {
        // __ldg: neighbor indices and packed atom data are read-only here, but the
        // compiler cannot prove they don't alias fx, so it won't use the read-only
        // cache on its own
        const int j        = __ldg(&neighs(neigh_neighbors, i, k, Nlocal, &neigh));
        const MD_FLOAT4 xqj = ldg_md_float4(&xq[j]);
        MD_FLOAT delx       = xqi.x - xqj.x;
        MD_FLOAT dely       = xqi.y - xqj.y;
        MD_FLOAT delz       = xqi.z - xqj.z;
        MD_FLOAT rsq        = delx * delx + dely * dely + delz * delz;

#if VL_PACK_PAR
        const MD_FLOAT4 pj = ldg_md_float4(&par[j]);
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
        const MD_FLOAT sigma6  = sigma3_i * pj.y;
        const MD_FLOAT epsilon = sqrt_eps_i * pj.x;
#elif LJ_COMB_RULE == LJ_COMB_FULL
        const int type_ij         = type_i * ntypes + md_float_as_int(pj.x);
        const MD_FLOAT cutforcesq = __ldg(&atom->cutforcesq[type_ij]);
        const MD_FLOAT sigma6     = __ldg(&atom->sigma6[type_ij]);
        const MD_FLOAT epsilon    = __ldg(&atom->epsilon[type_ij]);
#endif

        if (rsq < cutforcesq) {
#if COULOMB_RULE != COULOMB_NONE
            const MD_FLOAT qj = xqj.w;
#else
            const MD_FLOAT qj = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS
            const int mol_j = md_float_as_int(pj.z);
#else
            const int mol_j = 0;
#endif
            MD_FLOAT force = pairForceScalar(rsq, sigma6, epsilon, qi, qj, mol_i, mol_j, cs);
            MD_FLOAT partial_force_x = delx * force;
            MD_FLOAT partial_force_y = dely * force;
            MD_FLOAT partial_force_z = delz * force;

            if (j < Nlocal) { // ghost forces are not folded back
                atomicAdd(&atom_fx(j), -partial_force_x);
                atomicAdd(&atom_fy(j), -partial_force_y);
                atomicAdd(&atom_fz(j), -partial_force_z);
            }

            fix += partial_force_x;
            fiy += partial_force_y;
            fiz += partial_force_z;
        }
    };

    // Chunks of T neighbors; j-slices: block (bx, s) handles the chunks s, s+S, s+2S, ... of
    // its atoms (S = gridDim.y). Without slices the loop keeps unit stride, which nvcc
    // unrolls; the strided loop is not unrolled by itself and up to 46% slower at S = 1
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

    if (T > 1) {
        // the T lanes of an atom are consecutive lanes of one warp
        const unsigned int lane = threadIdx.x & 31;
        const unsigned int group = T >= 32 ? 0xffffffffu : ((1u << (T % 32)) - 1u);
        const unsigned int mask  = group << (lane & ~(unsigned int)(T - 1));
#pragma unroll
        for (int off = T / 2; off > 0; off /= 2) {
            fix += __shfl_down_sync(mask, fix, off, T);
            fiy += __shfl_down_sync(mask, fiy, off, T);
            fiz += __shfl_down_sync(mask, fiz, off, T);
        }
    }
    if (t == 0) {
        atomicAdd(&atom_fx(i), fix);
        atomicAdd(&atom_fy(i), fiy);
        atomicAdd(&atom_fz(i), fiz);
    }
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
static void packAtomsForForceVL(Atom* atom, MD_FLOAT4** xq, MD_FLOAT4** par)
{
    static MD_FLOAT4* d_xq  = NULL;
    static MD_FLOAT4* d_par = NULL;
    static int nalloc       = 0;
    static int packed_build = -1;
    const int n             = atom->Nlocal + atom->Nghost;
    int pack_par            = neighborBuildCountCUDA() != packed_build;
    if (n > nalloc) {
        nalloc = atom->Nmax > n ? atom->Nmax : n;
        d_xq   = (MD_FLOAT4*)reallocateGPU(d_xq, sizeof(MD_FLOAT4) * nalloc);
#if VL_PACK_PAR
        d_par = (MD_FLOAT4*)reallocateGPU(d_par, sizeof(MD_FLOAT4) * nalloc);
#endif
        pack_par = 1;
    }
    if (n > 0) {
        const int threads = 256;
        packAtomsVL<<<(n + threads - 1) / threads, threads>>>(atom->d_atom, n, d_xq, d_par, pack_par);
        cuda_assert("packAtomsVL", cudaPeekAtLastError());
    }
    packed_build = neighborBuildCountCUDA();
    *xq          = d_xq;
    *par         = d_par;
}

/* j-slices (grid.y) of the force kernels, one thread per atom and slice; N_b blocks without
 * slices, n_tot neighbor entries, bps resident blocks per SM:
 *     K > 0:  S = K * sqrt(N_SM * n_tot) / N_b    (as the cluster-pair kernels)
 *     K = 0:  S = W * N_SM * bps / N_b            (W waves of resident blocks)
 * at most SLICE_MAX_WAVES waves. Tuned on GH200 (argon 1k-1M atoms, r 0.9-2.8; mean loss
 * against the best S of each case 0.06-1.6%): the half kernel with SOA gains from slices up to
 * 125k atoms (per-pair j-force atomics), the other variants stop gaining at about 0.3 waves.
 * CSR (contiguous per atom like AOS) uses the AOS values, not tuned. */
#if defined(NBLIST_AOS) || defined(NBLIST_CSR)
#define VL_SLICE_K_HALF 0.0
#define VL_SLICE_W_HALF 0.35
#define VL_SLICE_K_FULL 0.0
#define VL_SLICE_W_FULL 0.23
#else
#define VL_SLICE_K_HALF 0.13
#define VL_SLICE_W_HALF 0.0
#define VL_SLICE_K_FULL 0.0
#define VL_SLICE_W_FULL 0.3
#endif
#ifndef SLICE_MAX_WAVES
#define SLICE_MAX_WAVES 10
#endif

// Resident blocks per SM of a force kernel: one wave is N_SM of these. Queried once per kernel.
static int residentBlocksPerSMVL(const void* kernel, int threads)
{
    static const void* kernels[16];
    static int blocks[16];
    static int nkernels = 0;
    for (int i = 0; i < nkernels; i++) {
        if (kernels[i] == kernel) return blocks[i];
    }
    int n = 0;
    cuda_assert("cudaOccupancyMaxActiveBlocksPerMultiprocessor",
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n, kernel, threads, 0));
    if (n < 1) n = 1;
    if (nkernels < 16) {
        kernels[nkernels] = kernel;
        blocks[nkernels]  = n;
        nkernels++;
    }
    return n;
}

static int chooseJSlicesVL(const Neighbor* neighbor, const void* kernel, int nblocks, int threads)
{
    static int forced = -2; // VL_JSLICES=S fixes S (S >= 1)
    if (forced == -2) {
        const char* e = getenv("VL_JSLICES");
        forced        = e != NULL ? atoi(e) : -1;
    }
    if (forced >= 1) {
        return forced;
    }
    static int nsm = 0;
    if (nsm == 0) {
        int dev = 0;
        cudaGetDevice(&dev);
        cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev);
    }
    if (nblocks <= 0) return 1;
    const int bps       = residentBlocksPerSMVL(kernel, threads);
    const double k      = neighbor->half_neigh ? VL_SLICE_K_HALF : VL_SLICE_K_FULL;
    const double w      = neighbor->half_neigh ? VL_SLICE_W_HALF : VL_SLICE_W_FULL;
    const double s      = k > 0.0 ? k * sqrt((double)nsm) * sqrt((double)numNeighborsTotalCUDA()) /
                                   (double)nblocks
                                  : w * nsm * bps / (double)nblocks;
    int S               = (int)lround(s);
    const long max_blocks = (long)SLICE_MAX_WAVES * nsm * bps;
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
    MD_FLOAT cutforcesq             = param->cutforce * param->cutforce;
    MD_FLOAT sigma6                 = param->sigma6;
    MD_FLOAT epsilon                = param->epsilon;
    const CoulombScalar cs          = coulombScalarInit(param);

    /*
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    size_t free, total;
    for(int i = 0; i < nDevices; ++i) {
        cudaMemGetInfo( &free, &total );
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("DEVICE %d/%d NAME: %s\r\n with %ld MB/%ld MB memory used", i + 1,
    nDevices, prop.name, free / 1024 / 1024, total / 1024 / 1024);
    }
    */

    // HINT: Run with cuda-memcheck ./MDBench-NVCC in case of error
    // memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nlocal * 3);

    GPU_PROFILE_START("force_lj");
    double S = getTimeStamp();
    LIKWID_MARKER_START("force");

    MD_FLOAT4 *xq, *par;
    packAtomsForForceVL(atom, &xq, &par);

    // Lanes per atom, fixed per list type (VL_TPA_LOG2_HALF/FULL in neighbor.h)
    const int half     = neighbor->half_neigh;
    const int lanes    = 1 << (half ? VL_TPA_LOG2_HALF : VL_TPA_LOG2_FULL);
    const void* kernel = half ? (const void*)computeForceLJCudaHalfNeigh
                              : (const void*)computeForceLJCudaFullNeigh;
    if (num_threads_per_block % lanes != 0) { // the lanes of an atom must stay in one block
        fprintf(stderr, "Error: NUM_THREADS must be a multiple of %d!\n", lanes);
        exit(EXIT_FAILURE);
    }
    const int bx = (int)(((long)Nlocal * lanes + num_threads_per_block - 1) /
                         num_threads_per_block);
    const dim3 grid(bx, chooseJSlicesVL(neighbor, kernel, bx, num_threads_per_block), 1);
    if (half || grid.y > 1) { // these add to the forces
#ifdef ATOM_POSITION_AOS
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax * 3);
#else
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fy, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fz, 0, sizeof(MD_FLOAT) * Nmax);
#endif
    }
    if (half) {
        computeForceLJCudaHalfNeigh<<<grid, num_threads_per_block>>>(
            atom->d_atom,
            cutforcesq,
            sigma6,
            epsilon,
            Nlocal,
            neighbor->d_neighbor,
            neighbor->d_neighbor.neighbors,
            neighbor->d_neighbor.numneigh,
            atom->ntypes,
            cs,
            xq,
            par);
    } else {
        computeForceLJCudaFullNeigh<<<grid, num_threads_per_block>>>(
            atom->d_atom,
            cutforcesq,
            sigma6,
            epsilon,
            Nlocal,
            neighbor->d_neighbor,
            neighbor->d_neighbor.neighbors,
            neighbor->d_neighbor.numneigh,
            atom->ntypes,
            cs,
            xq,
            par);
    }

    cuda_assert("computeForceLJCuda", cudaPeekAtLastError());
    cuda_assert("computeForceLJCuda", cudaDeviceSynchronize());
    GPU_PROFILE_STOP();

    LIKWID_MARKER_STOP("force");
    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJCUDA end\n");
    return E - S;
}
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
