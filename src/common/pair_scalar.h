/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Scalar force of one atom pair: Lennard-Jones, electrostatics (COULOMB_RULE) and
 * molecular exclusions (excl: the system has molecules), the same physics as pairForces()
 * in clusterpair/pair_rows.h computes for the cluster-pair SIMD kernels. Used by the GPU
 * cluster-pair kernels and by the verlet-list kernels (CPU and GPU). */
#ifndef __PAIR_SCALAR_H_
#define __PAIR_SCALAR_H_

#include <math.h>

#include <coulomb_scalar.h>
#include <parameter.h>

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
// acc (--accurate-recip): IEEE 1/x and 1/sqrt(x), else the fast hardware versions that
// -use_fast_math would give; HIP has no float rsqrt overload, so use the f versions
#if PRECISION == 1
#define PAIR_RECIP(x, acc) ((acc) ? 1.0f / (x) : __fdividef(1.0f, (x)))
#define PAIR_RSQRT(x, acc) ((acc) ? 1.0f / sqrtf(x) : rsqrtf(x))
#else
#define PAIR_RECIP(x, acc) (1.0 / (x))
#define PAIR_RSQRT(x, acc) ((acc) ? 1.0 / sqrt(x) : rsqrt(x))
#endif
#define PAIR_FMAX(a, b) fmax(a, b)
#else
#define PAIR_RECIP(x, acc) ((MD_FLOAT)1.0 / (x))
#define PAIR_RSQRT(x, acc) ((MD_FLOAT)1.0 / sqrt(x))
#define PAIR_FMAX(a, b)    ((a) > (b) ? (a) : (b))
#endif

/* Returns f_scalar of a pair inside the cut-off, the force on i is f_scalar * (r_i -
 * r_j). qi is the charge of i multiplied by epsfac. Atoms of the same molecule are
 * excluded: they skip Lennard-Jones and the direct Coulomb term, but still receive the
 * electrostatic correction term, as in GROMACS. */
static inline MD_HOST_DEVICE MD_FLOAT pairForceScalarAcc(MD_FLOAT rsq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    MD_FLOAT qi,
    MD_FLOAT qj,
    int mol_i,
    int mol_j,
    const int excl,
    const int accurate,
    CoulombScalar cs)
{
    // Padding atoms can sit at the same point. Like GROMACS, clamp the distance so that
    // their force is zero instead of NaN
    rsq = PAIR_FMAX(rsq, MIN_DISTANCE_SQ);

#if COULOMB_RULE == COULOMB_NONE
    MD_FLOAT sr2 = PAIR_RECIP(rsq, accurate);
#else
    MD_FLOAT rinv = PAIR_RSQRT(rsq, accurate);
    MD_FLOAT sr2  = rinv * rinv;
#endif
    MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
    MD_FLOAT force = (MD_FLOAT)(48.0) * sr6 * (sr6 - (MD_FLOAT)(0.5)) * sr2 * epsilon;

    const int interact = !excl || mol_i != mol_j;
    if (!interact) {
        force = (MD_FLOAT)0.0;
    }

#if COULOMB_RULE != COULOMB_NONE
    MD_FLOAT rinv_coul = rinv;
    if (!interact) {
        rinv_coul = (MD_FLOAT)0.0;
    }
    const MD_FLOAT qq = qi * qj;
#if COULOMB_RULE == COULOMB_RF
    force += qq * (rinv_coul * sr2 - cs.k_rf2);
#elif COULOMB_RULE == COULOMB_EWALD
    MD_FLOAT brsq = cs.beta2 * rsq;
    MD_FLOAT ewc  = cs.beta * pmeForceCorrectionScalar(brsq);
    force += qq * (ewc * brsq + rinv_coul) * sr2;
#else
    const MD_FLOAT r = rsq * rinv;
    force += qq * (rinv_coul - ewaldTableScalar(&cs, r, NULL) * r) * sr2;
#endif
#endif

    (void)qi;
    (void)qj;
    (void)mol_i;
    (void)mol_j;
    (void)interact;
    (void)cs;
    (void)accurate;
    return force;
}

// LJ and Coulomb energy of a pair inside the cut-off (see pairForceScalar())
static inline MD_HOST_DEVICE void pairEnergyScalarAcc(MD_FLOAT rsq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    MD_FLOAT qi,
    MD_FLOAT qj,
    int mol_i,
    int mol_j,
    const int excl,
    const int accurate,
    CoulombScalar cs,
    MD_FLOAT* vlj,
    MD_FLOAT* vc)
{
    rsq = PAIR_FMAX(rsq, MIN_DISTANCE_SQ);
#if COULOMB_RULE == COULOMB_NONE
    MD_FLOAT sr2 = PAIR_RECIP(rsq, accurate);
#else
    MD_FLOAT rinv = PAIR_RSQRT(rsq, accurate);
    MD_FLOAT sr2  = rinv * rinv;
#endif
    MD_FLOAT sr6 = sr2 * sr2 * sr2 * sigma6;
    MD_FLOAT e   = (MD_FLOAT)(4.0) * epsilon * sr6 * (sr6 - (MD_FLOAT)(1.0));

    const int interact = !excl || mol_i != mol_j;
    if (!interact) {
        e = (MD_FLOAT)0.0;
    }
    *vlj = e;

#if COULOMB_RULE != COULOMB_NONE
    MD_FLOAT rinv_coul = rinv;
    if (!interact) {
        rinv_coul = (MD_FLOAT)0.0;
    }
#if COULOMB_RULE == COULOMB_RF
    *vc = qi * qj * (rinv_coul + cs.k_rf * rsq - cs.c_rf);
#elif COULOMB_RULE == COULOMB_EWALD
    *vc = qi * qj * (rinv_coul - cs.beta * pmePotentialCorrectionScalar(cs.beta2 * rsq));
#else
    MD_FLOAT vcorr;
    ewaldTableScalar(&cs, rsq * rinv, &vcorr);
    *vc = qi * qj * (rinv_coul - vcorr);
#endif
#else
    *vc = (MD_FLOAT)0.0;
#endif

    (void)qi;
    (void)qj;
    (void)mol_i;
    (void)mol_j;
    (void)interact;
    (void)cs;
    (void)accurate;
}

// pairForceScalarAcc() with IEEE reciprocals, for the CPU kernels
static inline MD_HOST_DEVICE MD_FLOAT pairForceScalar(MD_FLOAT rsq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    MD_FLOAT qi,
    MD_FLOAT qj,
    int mol_i,
    int mol_j,
    const int excl,
    CoulombScalar cs)
{
    return pairForceScalarAcc(rsq, sigma6, epsilon, qi, qj, mol_i, mol_j, excl, 1, cs);
}

// pairEnergyScalarAcc() with IEEE reciprocals, for the CPU kernels
static inline MD_HOST_DEVICE void pairEnergyScalar(MD_FLOAT rsq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    MD_FLOAT qi,
    MD_FLOAT qj,
    int mol_i,
    int mol_j,
    const int excl,
    CoulombScalar cs,
    MD_FLOAT* vlj,
    MD_FLOAT* vc)
{
    pairEnergyScalarAcc(rsq, sigma6, epsilon, qi, qj, mol_i, mol_j, excl, 1, cs, vlj, vc);
}

#endif
