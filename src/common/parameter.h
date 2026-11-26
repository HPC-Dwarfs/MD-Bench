/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#ifndef __PARAMETER_H_
#define __PARAMETER_H_

#include <stdbool.h>
#include <stdint.h>

/* CPU SIMD kernels whose reciprocal is an exact division (SVE, and the scalar SIMD
 * emulation, see SIMD_RECIPROCAL_EXACT in simd_math.h): they are always accurate, and
 * --no-accurate-recip is turned into --accurate-recip with a warning. */
#if (defined(__ISA_SVE__) || defined(USE_SCALAR_KERNEL)) && !defined(CUDA_TARGET)
#define CPU_RECIPROCAL_EXACT
#ifdef __ISA_SVE__
#define RECIPROCAL_EXACT_TARGET "SVE"
#else
#define RECIPROCAL_EXACT_TARGET "the scalar SIMD emulation"
#endif
#define ACCURATE_RECIP_STR(param) "yes (exact division)"
#else
#define ACCURATE_RECIP_STR(param)                                                        \
    ((param)->accurate_reciprocal ? "yes" : "no (hardware estimate)")
#endif

// Portable vector types compatible with CUDA/HIP.
// Skip our definitions whenever CUDA/HIP runtime headers are (or may be)
// pulled in — including for .c files compiled by NVCC's host compiler,
// where __CUDACC__ is not defined but vector_types.h is still included
// transitively via device.h.
#if !defined(__CUDACC__) && !defined(__HIPCC__) && !defined(CUDA_TARGET)
typedef struct {
    float x, y, z;
} float3;
typedef struct {
    float x, y, z, w;
} float4;
typedef struct {
    double x, y, z;
} double3;
typedef struct {
    double x, y, z, w;
} double4;
#endif

#if PRECISION == 1
#define MD_FLOAT  float
#define MD_FLOAT3 float3
#define MD_FLOAT4 float4
#define MD_UINT   unsigned int
#else
#define MD_FLOAT  double
#define MD_FLOAT3 double3
#define MD_FLOAT4 double4
#define MD_UINT   unsigned long long int
#endif

// LJ combination rules (GROMACS terminology), param->lj_comb_rule (--lj-comb). The force
// kernels are compiled once per rule with -DLJ_COMB_RULE=<rule> (kernel_variant.h), all
// other code reads the rule at runtime and must not use LJ_COMB_RULE.
#define LJ_COMB_SINGLE 0 // Single atom type: broadcast global epsilon/sigma
#define LJ_COMB_GEOM   1 // Geometric: sqrt(eps_i*eps_j), sigma3_i*sigma3_j
#define LJ_COMB_MATRIX 2 // Type-pair matrix lookup (no combination formula)
#define LJ_COMB_AUTO   3 // Inferred from the atom types (resolveLJCombRule())

// Verlet-list CPU force kernel (--force-simd) and neighbor-list build (--neighbor-simd)
#define SIMD_MODE_AUTOVEC             0 // Plain C, vectorized by the compiler
#define SIMD_MODE_INTRINSICS          1 // Hand-written SIMD intrinsics
#define SIMD_MODE_INTRINSICS_COMPRESS 2 // Intrinsics, compressed cut-off lanes (force)

// Index variable of the tabulated LJ forces, param->lj_table_index (--lj-table-index)
#define LJ_TABLE_INDEX_R   0 // Uniform grid in r (matches GROMACS), a sqrt per pair
#define LJ_TABLE_INDEX_RSQ 1 // Uniform grid in r^2, no sqrt

// Coulomb treatment compile-time macros
// Use -DCOULOMB_RULE=<value> at compile time
#define COULOMB_NONE        0 // No electrostatics (LJ only)
#define COULOMB_RF          1 // Reaction-field with infinite dielectric, shifted at the cut-off
#define COULOMB_EWALD       2 // Real-space Ewald with analytical force correction
#define COULOMB_EWALD_TABLE 3 // Real-space Ewald with tabulated (cubic spline) correction

#ifndef COULOMB_RULE
#define COULOMB_RULE COULOMB_NONE
#endif

// Both Ewald variants share the splitting coefficient and the energy terms
#define COULOMB_IS_EWALD                                                                 \
    (COULOMB_RULE == COULOMB_EWALD || COULOMB_RULE == COULOMB_EWALD_TABLE)

// Lower limit of the squared distance in the force kernels, as GROMACS'
// c_nbnxmMinDistanceSquared: below it 1/r^12 overflows. The padding atoms of two clusters
// are at the same point, the limit keeps their (zero) force finite.
#if PRECISION == 1
#define MIN_DISTANCE_SQ ((MD_FLOAT)3.82e-07)
#else
#define MIN_DISTANCE_SQ ((MD_FLOAT)1.0e-36)
#endif

#if COULOMB_RULE == COULOMB_NONE
#define COULOMB_RULE_NAME "none"
#elif COULOMB_RULE == COULOMB_RF
#define COULOMB_RULE_NAME "reaction-field"
#elif COULOMB_RULE == COULOMB_EWALD
#define COULOMB_RULE_NAME "ewald"
#else
#define COULOMB_RULE_NAME "ewald (tabulated correction)"
#endif

typedef struct {
    int force_field;
    char* param_file;
    char* input_file;
    char* vtk_file;
    char* xtc_file;
    char* write_atom_file;
    char* types_file;
    MD_FLOAT epsilon;
    MD_FLOAT sigma;
    MD_FLOAT sigma6;
    MD_FLOAT temp;
    MD_FLOAT rho;
    MD_FLOAT mass;
    int ntypes;
    MD_FLOAT* epsilon_per_type;
    MD_FLOAT* sigma_per_type;
    int ntimes;
    int nstat;
    int reneigh_every;
    int displacement_reneigh;
    int resort_every; // --sort-every: 0 off, < 0 at every neighbor-list build
    int prune_every;
    int x_out_every;
    int v_out_every;
    int half_neigh;
    MD_FLOAT dt;
    MD_FLOAT dtforce;
    MD_FLOAT skin;
    MD_FLOAT outer_skin;
    MD_FLOAT cutforce;
    MD_FLOAT cutneigh;
    int nx, ny, nz;
    int pbc_x, pbc_y, pbc_z;
    MD_FLOAT lattice;
    MD_FLOAT xlo, xhi, ylo, yhi, zlo, zhi;
    MD_FLOAT xprd, yprd, zprd;
    double proc_freq;
    int super_clustering;
    int ref_kernel; // GPU force kernels: use the reference (simple) version
    char* eam_file;
    int lj_table_points;
    // MPI implementation
    int balance;
    int method;
    int balance_every;
    int setup;
    int verbose;
    int size;
    bool gmxbenchmark;
    // Electrostatics (only used when COULOMB_RULE != COULOMB_NONE)
    MD_FLOAT epsfac;     // Coulomb prefactor (1/(4 pi eps0) in the unit system used)
    MD_FLOAT ewald_rtol; // Relative tolerance at the cut-off used to derive ewald_beta
    MD_FLOAT ewald_beta; // Ewald splitting coefficient
    // Ewald correction table of COULOMB_EWALD_TABLE (see ewald_table.h)
    MD_FLOAT ewald_table_scale; // points per length unit
    int ewald_table_size;
    MD_FLOAT* ewald_table_f;    // force correction F
    MD_FLOAT* ewald_table_v;    // potential correction V
    MD_FLOAT* ewald_table_fdv0; // F, F[i+1] - F[i], V, 0 per point
    // Kernel-only benchmark: run the force kernel this many times on a fixed list
    int kernel_bench;
    char* dump_forces_file;
    // Energy kernels (-energy): LJ and Coulomb energy of the last force call
    bool compute_energy;
    double energy_lj;
    double energy_coul;
    // Half-LJ i-clusters (-halflj), as GROMACS' half-LJ kernels
    bool half_lj;
    // Atoms per molecule of the system, > 1 turns on the molecular exclusions: atoms of
    // the same molecule only interact through the Coulomb exclusion correction
    int molsize;
    // LJ combination rule (LJ_COMB_*), LJ_COMB_AUTO until resolveLJCombRule()
    int lj_comb_rule;
    // Newton-Raphson refined 1/r^2 and 1/r in the SIMD kernels (--accurate-recip)
    int accurate_reciprocal;
    // j-cluster prefetch of the CPU cluster-pair SIMD kernels (--no-prefetch disables)
    int kernel_prefetch;
    // Verlet-list CPU force kernel and neighbor-list build (SIMD_MODE_*)
    int force_simd;
    int neighbor_simd;
    // Index variable of the tabulated LJ forces (LJ_TABLE_INDEX_*)
    int lj_table_index;
    // j-slices (grid.y) of the GPU force kernels, 0 chooses them per force call (auto)
    int gpu_j_slices;
    // Pinned host memory for the buffers transferred to/from the GPU
    int pinned_memory;
} Parameter;

void initParameter(Parameter*);
void readParameter(Parameter*, const char*);
void readTypesFile(Parameter*);
void printParameter(Parameter*);
int str2ljcomb(const char*);
const char* ljcomb2str(int);
int str2simdmode(const char*);
const char* simdmode2str(int);
int str2ljtableindex(const char*);
const char* ljtableindex2str(int);
int str2jslices(const char*);
void resolveLJCombRule(Parameter*,
    int ntypes,
    const MD_FLOAT* sqrt_epsilon_per_type,
    const MD_FLOAT* sigma3_per_type);
double computeEwaldCoefficient(double cutoff, double rtol);
void computePerTypeLJParameters(int, Parameter*, MD_FLOAT*, MD_FLOAT*);
void computeTypePairLJParameters(int, MD_FLOAT*, MD_FLOAT*, MD_FLOAT*, MD_FLOAT*);

#endif
