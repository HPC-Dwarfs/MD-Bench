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

// par[] is only read with geometric LJ or exclusions (the full rule reads mol directly)
#if LJ_COMB_RULE == LJ_COMB_GEOM || (USE_EXCLUSIONS && LJ_COMB_RULE != LJ_COMB_FULL)
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

// xq = {x, y, z, type (full rule) or charge}, par = {sqrt_epsilon, sigma3, mol, 0}
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
#if LJ_COMB_RULE == LJ_COMB_FULL
    const MD_FLOAT w = int_as_md_float(atom->type[i]);
#elif COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT w = atom->q[i];
#else
    const MD_FLOAT w = (MD_FLOAT)0.0;
#endif
    xq[i] = make_md_float4(atom_x(i), atom_y(i), atom_z(i), w);
#if VL_PACK_PAR
    if (pack_par) {
#if LJ_COMB_RULE == LJ_COMB_GEOM
        const MD_FLOAT p0 = atom->sqrt_epsilon[i];
        const MD_FLOAT p1 = atom->sigma3[i];
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

/* Force kernels: T = 2^LGT lanes per atom (TPA, as in the LAMMPS GPU package), T = 16
 * with full and 8 with half lists. Lane t of atom i handles neighbors t, t+T, t+2T, ...;
 * the list is laid out so that these T lanes read consecutive entries (neighs()). The
 * host keeps the T lanes of an atom consecutive in one warp.
 *
 * Reduction: a tree over log2(T) steps into lane 0 (lane t += lane t+off). Only lane 0's
 * result is used, and it only reads lanes of its own atom.
 *  - NVIDIA: shuffle mask = this atom's T lanes only (atoms i >= Nlocal have returned)
 *  - AMD: DPP row shifts, which stay inside a 16-lane row, hence T <= 16
 */
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

    constexpr int LGT = VL_TPA_LOG2_FULL; // T = 2^LGT lanes per atom, see above
    constexpr int T   = 1 << LGT;
    const int gtid    = blockIdx.x * blockDim.x + threadIdx.x;
    const int i       = gtid >> LGT;
    const int t       = gtid & (T - 1);
#if LJ_COMB_RULE == LJ_COMB_FULL
    extern __shared__ MD_FLOAT4 lj_table[];
    for (int k = threadIdx.x; k < ntypes * ntypes; k += blockDim.x) {
        lj_table[k] = make_md_float4(a.cutforcesq[k], a.sigma6[k], a.epsilon[k], (MD_FLOAT)0.0);
    }
    __syncthreads();
#endif
    if (i >= Nlocal) {
        return; 
    }
    // The host set the same value; as a constant, neighs() indexes with immediate shifts
    neigh.tpa_log2 = LGT;

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
    const int row_i = md_float_as_int(xqi.w) * ntypes; // row of i in lj_table
#endif
#if COULOMB_RULE != COULOMB_NONE && LJ_COMB_RULE == LJ_COMB_FULL
    const MD_FLOAT qi = cs.epsfac * atom->q[i];
#elif COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT qi = cs.epsfac * xqi.w;
#else
    const MD_FLOAT qi = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS && LJ_COMB_RULE == LJ_COMB_FULL
    const int mol_i = atom->mol[i];
#elif USE_EXCLUSIONS
    const int mol_i = md_float_as_int(pi.z);
#else
    const int mol_i = 0;
#endif

    auto pair = [&](int k) {
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
        const MD_FLOAT4 lj        = lj_table[row_i + md_float_as_int(xqj.w)];
        const MD_FLOAT cutforcesq = lj.x;
        const MD_FLOAT sigma6     = lj.y;
        const MD_FLOAT epsilon    = lj.z;
#endif

        if (rsq < cutforcesq) {
#if COULOMB_RULE != COULOMB_NONE && LJ_COMB_RULE == LJ_COMB_FULL
            const MD_FLOAT qj = __ldg(&atom->q[j]);
#elif COULOMB_RULE != COULOMB_NONE
            const MD_FLOAT qj = xqj.w;
#else
            const MD_FLOAT qj = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS && LJ_COMB_RULE == LJ_COMB_FULL
            const int mol_j = __ldg(&atom->mol[j]);
#elif USE_EXCLUSIONS
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
        const VL_LANE_MASK group = T >= VL_WARP_SIZE ? ~(VL_LANE_MASK)0
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
        if (gridDim.y == 1) { 
            atom_fx(i) = fix;
            atom_fy(i) = fiy;
            atom_fz(i) = fiz;
        } else {
            atomicAdd(&atom_fx(i), fix);
            atomicAdd(&atom_fy(i), fiy);
            atomicAdd(&atom_fz(i), fiz);
        }
    }
}


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
#if LJ_COMB_RULE == LJ_COMB_FULL
    extern __shared__ MD_FLOAT4 lj_table[];
    for (int k = threadIdx.x; k < ntypes * ntypes; k += blockDim.x) {
        lj_table[k] = make_md_float4(a.cutforcesq[k], a.sigma6[k], a.epsilon[k], (MD_FLOAT)0.0);
    }
    __syncthreads();
#endif
    if (i >= Nlocal) {
        return; 
    }
    neigh.tpa_log2 = LGT;

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
    const int row_i = md_float_as_int(xqi.w) * ntypes; // row of i in lj_table
#endif
#if COULOMB_RULE != COULOMB_NONE && LJ_COMB_RULE == LJ_COMB_FULL
    const MD_FLOAT qi = cs.epsfac * atom->q[i]; 
#elif COULOMB_RULE != COULOMB_NONE
    const MD_FLOAT qi = cs.epsfac * xqi.w;
#else
    const MD_FLOAT qi = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS && LJ_COMB_RULE == LJ_COMB_FULL
    const int mol_i = atom->mol[i];
#elif USE_EXCLUSIONS
    const int mol_i = md_float_as_int(pi.z);
#else
    const int mol_i = 0;
#endif

    auto pair = [&](int k) {
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
        const MD_FLOAT4 lj        = lj_table[row_i + md_float_as_int(xqj.w)];
        const MD_FLOAT cutforcesq = lj.x;
        const MD_FLOAT sigma6     = lj.y;
        const MD_FLOAT epsilon    = lj.z;
#endif

        if (rsq < cutforcesq) {
#if COULOMB_RULE != COULOMB_NONE && LJ_COMB_RULE == LJ_COMB_FULL
            const MD_FLOAT qj = __ldg(&atom->q[j]);
#elif COULOMB_RULE != COULOMB_NONE
            const MD_FLOAT qj = xqj.w;
#else
            const MD_FLOAT qj = (MD_FLOAT)0.0;
#endif
#if USE_EXCLUSIONS && LJ_COMB_RULE == LJ_COMB_FULL
            const int mol_j = __ldg(&atom->mol[j]);
#elif USE_EXCLUSIONS
            const int mol_j = md_float_as_int(pj.z);
#else
            const int mol_j = 0;
#endif
            MD_FLOAT force = pairForceScalar(rsq, sigma6, epsilon, qi, qj, mol_i, mol_j, cs);
            MD_FLOAT partial_force_x = delx * force;
            MD_FLOAT partial_force_y = dely * force;
            MD_FLOAT partial_force_z = delz * force;

            if (j < Nlocal) { /
                atomicAdd(&atom_fx(j), -partial_force_x);
                atomicAdd(&atom_fy(j), -partial_force_y);
                atomicAdd(&atom_fz(j), -partial_force_z);
            }

            fix += partial_force_x;
            fiy += partial_force_y;
            fiz += partial_force_z;
        }
    };

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
#ifdef __HIP_PLATFORM_AMD__
        static_assert(T <= 16, "the lanes of an atom must lie in one DPP row");
        if (T > 8) VL_ROW_STEP(fix, fiy, fiz, 8);
        if (T > 4) VL_ROW_STEP(fix, fiy, fiz, 4);
        if (T > 2) VL_ROW_STEP(fix, fiy, fiz, 2);
        VL_ROW_STEP(fix, fiy, fiz, 1);
#else
        const unsigned int lane  = threadIdx.x & (VL_WARP_SIZE - 1);
        const VL_LANE_MASK group = T >= VL_WARP_SIZE ? ~(VL_LANE_MASK)0
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

// j-slices (grid.y) of the force kernels: 0 = chosen per force call (auto), > 0 = fixed
#ifndef GPU_J_SLICES
#define GPU_J_SLICES 0
#endif

#if GPU_J_SLICES == 0
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
#endif

// j-slices (grid.y) of this force call, T threads per atom and slice, see GPU_J_SLICES
static int chooseJSlicesVL(
    const Neighbor* neighbor, const void* kernel, int nblocks, int threads, size_t shmem)
{
#if GPU_J_SLICES > 0
    return GPU_J_SLICES;
#else
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
        cudaGetDevice(&dev);
        cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev);
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
#endif
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

    // Lanes per atom T, see the comment above computeForceLJCudaFullNeigh
    const int half     = neighbor->half_neigh;
    const int lanes    = 1 << (half ? VL_TPA_LOG2_HALF : VL_TPA_LOG2_FULL);
    const void* kernel = half ? (const void*)computeForceLJCudaHalfNeigh
                              : (const void*)computeForceLJCudaFullNeigh;
    // the lanes of an atom must stay in one block (and so in one warp) for the reduction
    if (num_threads_per_block % lanes != 0) {
        fprintf(stderr, "Error: NUM_THREADS must be a multiple of %d!\n", lanes);
        exit(EXIT_FAILURE);
    }
    const size_t shmem = LJ_COMB_RULE == LJ_COMB_FULL // the type-pair table
                             ? (size_t)atom->ntypes * atom->ntypes * sizeof(MD_FLOAT4)
                             : 0;
    if (shmem > 48 * 1024) {
        fprintf(stderr, "Error: the table of %d atom types does not fit in shared memory!\n",
            atom->ntypes);
        exit(EXIT_FAILURE);
    }
    // grid.x: T threads per local atom; grid.y: j-slices
    const int bx = (int)(((long)Nlocal * lanes + num_threads_per_block - 1) /
                         num_threads_per_block);
    const dim3 grid(bx, chooseJSlicesVL(neighbor, kernel, bx, num_threads_per_block, shmem), 1);
    if (half || grid.y > 1) { // the kernel accumulates with atomicAdd
#ifdef ATOM_POSITION_AOS
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax * 3);
#else
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fy, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fz, 0, sizeof(MD_FLOAT) * Nmax);
#endif
    }
    if (half) {
        computeForceLJCudaHalfNeigh<<<grid, num_threads_per_block, shmem>>>(
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
        computeForceLJCudaFullNeigh<<<grid, num_threads_per_block, shmem>>>(
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

// Reference force kernel: one thread per atom, no TPA

// Charge and molecule of atom i, 0 if the build does not use them
__device__ __forceinline__ MD_FLOAT chargeVL(const DeviceAtom* atom, int i)
{
#if COULOMB_RULE != COULOMB_NONE
    return atom->q[i];
#else
    return (MD_FLOAT)0.0;
#endif
}
__device__ __forceinline__ int molVL(const DeviceAtom* atom, int i)
{
#if USE_EXCLUSIONS
    return atom->mol[i];
#else
    return 0;
#endif
}

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
    CoulombScalar cs)
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

    for (int k = 0; k < numneighs; k++) {
        // neighs() hides the list layout (AOS, SOA with lanes per atom, CSR)
        int j         = neighs(neigh_neighbors, i, k, Nlocal, &neigh);
        MD_FLOAT delx = xtmp - atom_x(j);
        MD_FLOAT dely = ytmp - atom_y(j);
        MD_FLOAT delz = ztmp - atom_z(j);
        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
        MD_FLOAT sigma6  = atom->sigma3[i] * atom->sigma3[j];
        MD_FLOAT epsilon = atom->sqrt_epsilon[i] * atom->sqrt_epsilon[j];
#elif LJ_COMB_RULE == LJ_COMB_FULL
        int type_ij         = atom->type[i] * ntypes + atom->type[j];
        MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
        MD_FLOAT sigma6     = atom->sigma6[type_ij];
        MD_FLOAT epsilon    = atom->epsilon[type_ij];
#endif

        if (rsq < cutforcesq) {
            MD_FLOAT force = pairForceScalar(rsq,
                sigma6,
                epsilon,
                cs.epsfac * chargeVL(atom, i),
                chargeVL(atom, j),
                molVL(atom, i),
                molVL(atom, j),
                cs);
            MD_FLOAT fx    = delx * force;
            MD_FLOAT fy    = dely * force;
            MD_FLOAT fz    = delz * force;

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
}

extern "C" double computeForceLJCUDARef(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCUDARef begin\n");

    const int num_threads_per_block = get_cuda_num_threads();
    int Nlocal                      = atom->Nlocal;
    int Nmax                        = atom->Nmax;
    MD_FLOAT cutforcesq             = param->cutforce * param->cutforce;
    MD_FLOAT sigma6                 = param->sigma6;
    MD_FLOAT epsilon                = param->epsilon;
    const CoulombScalar cs          = coulombScalarInit(param);

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
    computeForceLJCudaRef_warp<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
        cutforcesq,
        sigma6,
        epsilon,
        Nlocal,
        neighbor->d_neighbor,
        neighbor->d_neighbor.neighbors,
        neighbor->d_neighbor.numneigh,
        atom->ntypes,
        neighbor->half_neigh,
        cs);

    cuda_assert("computeForceLJCudaRef", cudaPeekAtLastError());
    cuda_assert("computeForceLJCudaRef", cudaDeviceSynchronize());
    GPU_PROFILE_STOP();

    LIKWID_MARKER_STOP("force");
    double E = getTimeStamp();
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
