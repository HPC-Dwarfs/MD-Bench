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
/*
#ifdef USE_REFERENCE_KERNEL
#define MD_SIMD_FLOAT float
#define MD_SIMD_MASK  uint16_t
#endif
*/
#else
#define MD_FLOAT  double
#define MD_FLOAT3 double3
#define MD_FLOAT4 double4
#define MD_UINT   unsigned long long int
/*
#ifdef USE_REFERENCE_KERNEL
#define MD_SIMD_FLOAT double
#define MD_SIMD_MASK  uint8_t
#endif
*/
#endif

// LJ combination rule compile-time macros (Gromacs terminology)
// Use -DLJ_COMB_RULE=<value> at compile time
#define LJ_COMB_SINGLE 0 // Single atom type: broadcast global epsilon/sigma
#define LJ_COMB_GEOM   1 // Geometric: sqrt(eps_i*eps_j), sigma3_i*sigma3_j
#define LJ_COMB_FULL   2 // Full type-pair matrix lookup (no combination formula)

// Default to geometric combination rule if not specified
#ifndef LJ_COMB_RULE
#define LJ_COMB_RULE LJ_COMB_GEOM
#endif

// String names for printing
#if LJ_COMB_RULE == LJ_COMB_SINGLE
#define LJ_COMB_RULE_NAME "single"
#elif LJ_COMB_RULE == LJ_COMB_GEOM
#define LJ_COMB_RULE_NAME "geometric"
#else
#define LJ_COMB_RULE_NAME "full"
#endif

// Coulomb treatment compile-time macros
// Use -DCOULOMB_RULE=<value> at compile time
#define COULOMB_NONE  0 // No electrostatics (LJ only)
#define COULOMB_RF    1 // Reaction-field with infinite dielectric, shifted at the cut-off
#define COULOMB_EWALD 2 // Real-space Ewald with analytical force correction

#ifndef COULOMB_RULE
#define COULOMB_RULE COULOMB_NONE
#endif

// Molecular exclusions: atoms of the same molecule do not interact via LJ and the
// direct Coulomb term, only through the Coulomb exclusion correction
#ifndef USE_EXCLUSIONS
#define USE_EXCLUSIONS 0
#endif

#if COULOMB_RULE == COULOMB_NONE
#define COULOMB_RULE_NAME "none"
#elif COULOMB_RULE == COULOMB_RF
#define COULOMB_RULE_NAME "reaction-field"
#else
#define COULOMB_RULE_NAME "ewald"
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
    int resort_every;
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
    // Kernel-only benchmark: run the force kernel this many times on a fixed list
    int kernel_bench;
    char* dump_forces_file;
} Parameter;

void initParameter(Parameter*);
void readParameter(Parameter*, const char*);
void readTypesFile(Parameter*);
void printParameter(Parameter*);
double computeEwaldCoefficient(double cutoff, double rtol);
void computePerTypeLJParameters(int, Parameter*, MD_FLOAT*, MD_FLOAT*);
void computeTypePairLJParameters(int, MD_FLOAT*, MD_FLOAT*, MD_FLOAT*, MD_FLOAT*);

#endif
