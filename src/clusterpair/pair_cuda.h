/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Scalar force of one atom pair for the GPU cluster-pair kernels, the same physics as
 * pairForces() in pair_rows.h computes for the CPU kernels: Lennard-Jones,
 * electrostatics (COULOMB_RULE) and molecular exclusions (USE_EXCLUSIONS). */
#ifndef __PAIR_CUDA_H_
#define __PAIR_CUDA_H_

#include <atom.h>
#include <coulomb_scalar.h>
#include <parameter.h>

/* Charge of the cluster slot sca, 0 without electrostatics (cl_q is then not on the GPU) */
static __device__ __forceinline__ MD_FLOAT chargeCuda(const MD_FLOAT* cl_q, int sca)
{
#if COULOMB_RULE != COULOMB_NONE
    return cl_q[sca];
#else
    (void)cl_q;
    (void)sca;
    return (MD_FLOAT)0.0;
#endif
}

/* Molecule of the cluster slot sca, 0 without exclusions (cl_mol is then not on the GPU) */
static __device__ __forceinline__ int moleculeCuda(const int* cl_mol, int sca)
{
#if USE_EXCLUSIONS
    return cl_mol[sca];
#else
    (void)cl_mol;
    (void)sca;
    return 0;
#endif
}

/* Returns f_scalar of a pair inside the cut-off, the force on i is f_scalar * (r_i - r_j).
 * qi is the charge of i multiplied by epsfac. Atoms of the same molecule are excluded:
 * they skip Lennard-Jones and the direct Coulomb term, but still receive the
 * electrostatic correction term, as in GROMACS. */
static __device__ __forceinline__ MD_FLOAT pairForceCuda(MD_FLOAT rsq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    MD_FLOAT qi,
    MD_FLOAT qj,
    int mol_i,
    int mol_j,
    CoulombScalar cs)
{
    // Padding atoms of two clusters can sit at the same point. Like GROMACS, clamp the
    // distance so that their force is zero instead of NaN
    rsq = fmax(rsq, MIN_DISTANCE_SQ);

#if COULOMB_RULE == COULOMB_NONE
    MD_FLOAT sr2 = (MD_FLOAT)(1.0) / rsq;
#else
    MD_FLOAT rinv = rsqrt(rsq);
    MD_FLOAT sr2  = rinv * rinv;
#endif
    MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
    MD_FLOAT force = (MD_FLOAT)(48.0) * sr6 * (sr6 - (MD_FLOAT)(0.5)) * sr2 * epsilon;

#if USE_EXCLUSIONS
    const int interact = mol_i != mol_j;
    if (!interact) {
        force = (MD_FLOAT)0.0;
    }
#endif

#if COULOMB_RULE != COULOMB_NONE
    MD_FLOAT rinv_coul = rinv;
#if USE_EXCLUSIONS
    if (!interact) {
        rinv_coul = (MD_FLOAT)0.0;
    }
#endif
    const MD_FLOAT qq = qi * qj;
#if COULOMB_RULE == COULOMB_RF
    force += qq * (rinv_coul * sr2 - cs.k_rf2);
#else
    MD_FLOAT brsq = cs.beta2 * rsq;
    MD_FLOAT ewc  = cs.beta * pmeForceCorrectionScalar(brsq);
    force += qq * (ewc * brsq + rinv_coul) * sr2;
#endif
#endif

    (void)qi;
    (void)qj;
    (void)mol_i;
    (void)mol_j;
    (void)cs;
    return force;
}

#endif
