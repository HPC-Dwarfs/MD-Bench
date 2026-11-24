/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Scalar electrostatics of one atom pair, usable from C and from CUDA/HIP kernels.
 * Same formulation as the SIMD version in clusterpair/coulomb.h: a pair contributes
 *   f_scalar = qq * (interact * 1/r^3 + <correction>)      (f = f_scalar * (r_i - r_j))
 *   reaction-field (eps_rf = inf):  correction = -2 k_rf,   k_rf = 1/(2 rc^3)
 *   Ewald (real space):             correction = beta^3 * pmeForceCorrection(beta^2 r^2)
 * The Ewald force correction is the analytical rational approximation used by GROMACS
 * (src/gromacs/simd/.../scalar_math.h, LGPL).
 */
#ifndef __COULOMB_SCALAR_H_
#define __COULOMB_SCALAR_H_

#include <parameter.h>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MD_HOST_DEVICE __host__ __device__
#else
#define MD_HOST_DEVICE
#endif

typedef struct {
    MD_FLOAT epsfac; // Coulomb prefactor, the charges are multiplied by it
    MD_FLOAT k_rf2;  // 2 k_rf = 1/rc^3
    MD_FLOAT beta;   // Ewald splitting coefficient
    MD_FLOAT beta2;
} CoulombScalar;

static inline CoulombScalar coulombScalarInit(Parameter* param)
{
    CoulombScalar cs;
    cs.epsfac = param->epsfac;
    cs.k_rf2  = 1.0 / (param->cutforce * param->cutforce * param->cutforce);
    cs.beta   = param->ewald_beta;
    cs.beta2  = param->ewald_beta * param->ewald_beta;
    return cs;
}

/* Ewald force correction as a function of z2 = (beta r)^2, scalar copy of
 * pmeForceCorrectionSimd() in clusterpair/coulomb.h */
static inline MD_HOST_DEVICE MD_FLOAT pmeForceCorrectionScalar(MD_FLOAT z2)
{
    const MD_FLOAT z4 = z2 * z2;
#if PRECISION == 1
    MD_FLOAT polyFD0 = (MD_FLOAT)0.0011193462567257629232 * z4 +
                       (MD_FLOAT)0.11583842382862377919;
    MD_FLOAT polyFD1 = (MD_FLOAT)0.014866955030185295499 * z4 +
                       (MD_FLOAT)0.50736591960530292870;
    polyFD0          = polyFD0 * z4 + (MD_FLOAT)1.0;
    polyFD0          = polyFD1 * z2 + polyFD0;

    MD_FLOAT polyFN0 = (MD_FLOAT)-1.7357322914161492954e-8 * z4 +
                       (MD_FLOAT)-0.000053401640219807709149;
    MD_FLOAT polyFN1 = (MD_FLOAT)1.4703624142580877519e-6 * z4 +
                       (MD_FLOAT)0.0010054721316683106153;
    polyFN0          = polyFN0 * z4 + (MD_FLOAT)-0.019278317264888380590;
    polyFN1          = polyFN1 * z4 + (MD_FLOAT)0.069670166153766424023;
    polyFN0          = polyFN0 * z4 + (MD_FLOAT)-0.75225204789749321333;
    polyFN0          = polyFN1 * z2 + polyFN0;
#else
    MD_FLOAT polyFD1 = (MD_FLOAT)0.000016009278224355026701 * z4 +
                       (MD_FLOAT)0.0081803507497974289008;
    polyFD1          = polyFD1 * z4 + (MD_FLOAT)0.41543303143712535988;
    polyFD1          = polyFD1 * z2;
    MD_FLOAT polyFD0 = (MD_FLOAT)0.00051055686934806966046 * z4 +
                       (MD_FLOAT)0.077181146026670287235;
    polyFD0          = polyFD0 * z4 + (MD_FLOAT)1.0;
    polyFD0          = polyFD0 + polyFD1;

    MD_FLOAT polyFN0 = (MD_FLOAT)-8.0072854618360083154e-14 * z4 +
                       (MD_FLOAT)-8.1490406329798423616e-10;
    polyFN0          = polyFN0 * z4 + (MD_FLOAT)-9.9471420832602741006e-7;
    polyFN0          = polyFN0 * z4 + (MD_FLOAT)-0.00031991745139313364005;
    polyFN0          = polyFN0 * z4 + (MD_FLOAT)-0.031750380176100813405;
    polyFN0          = polyFN0 * z4 + (MD_FLOAT)-0.75225277815249618847;
    MD_FLOAT polyFN1 = (MD_FLOAT)1.1859116242260148027e-11 * z4 +
                       (MD_FLOAT)3.4404793543907847655e-8;
    polyFN1          = polyFN1 * z4 + (MD_FLOAT)0.000020740315999115847456;
    polyFN1          = polyFN1 * z4 + (MD_FLOAT)0.0035074449373659008203;
    polyFN1          = polyFN1 * z4 + (MD_FLOAT)0.13884101728898463426;
    polyFN0          = polyFN1 * z2 + polyFN0;
#endif
    return polyFN0 / polyFD0;
}

#endif
