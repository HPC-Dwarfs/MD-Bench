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
 *   Ewald (tabulated):              correction = -Fcorr(r) / r, see ewald_table.h
 * Energy V = qq (interact/r - <potential correction>), plus the self energy.
 * The Ewald corrections are the analytical rational approximations used by GROMACS
 * (src/gromacs/simd/.../scalar_math.h, LGPL).
 */
#ifndef __COULOMB_SCALAR_H_
#define __COULOMB_SCALAR_H_

#include <math.h>

#include <parameter.h>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MD_HOST_DEVICE __host__ __device__
#else
#define MD_HOST_DEVICE
#endif

typedef struct {
    MD_FLOAT epsfac; // Coulomb prefactor, the charges are multiplied by it
    MD_FLOAT k_rf2;  // 2 k_rf = 1/rc^3
    MD_FLOAT k_rf;   // k_rf = 1/(2 rc^3)
    MD_FLOAT c_rf;   // reaction-field potential shift, 1/rc + k_rf rc^2
    MD_FLOAT beta;   // Ewald splitting coefficient
    MD_FLOAT beta2;
    MD_FLOAT tab_scale;    // Ewald table points per length unit
    const MD_FLOAT* tab_f; // Ewald table force correction
    const MD_FLOAT* tab_v; // Ewald table potential correction
} CoulombScalar;

static inline CoulombScalar coulombScalarInit(Parameter* param)
{
    CoulombScalar cs;
    const MD_FLOAT rc = param->cutforce;
    cs.epsfac         = param->epsfac;
    cs.k_rf2          = 1.0 / (rc * rc * rc);
    cs.k_rf           = 0.5 * cs.k_rf2;
    cs.c_rf           = 1.0 / rc + cs.k_rf * rc * rc;
    cs.beta           = param->ewald_beta;
    cs.beta2          = param->ewald_beta * param->ewald_beta;
    cs.tab_scale      = param->ewald_table_scale;
    cs.tab_f          = param->ewald_table_f;
    cs.tab_v          = param->ewald_table_v;
    return cs;
}

// Self energy factor: each atom adds -epsfac q^2 times it
static inline double coulombSelfEnergy(const CoulombScalar* cs)
{
#if COULOMB_RULE == COULOMB_RF
    return 0.5 * cs->c_rf;
#elif COULOMB_RULE == COULOMB_EWALD
    return cs->beta * M_2_SQRTPI * 0.5;
#elif COULOMB_RULE == COULOMB_EWALD_TABLE
    return 0.5 * cs->tab_v[0];
#else
    (void)cs;
    return 0.0;
#endif
}

// Ewald table force correction at r (and potential correction into *v), see ewald_table.h
static inline MD_HOST_DEVICE MD_FLOAT ewaldTableScalar(
    const CoulombScalar* cs, MD_FLOAT r, MD_FLOAT* v)
{
    const MD_FLOAT rs   = r * cs->tab_scale;
    const int ri        = (int)rs;
    const MD_FLOAT frac = rs - (MD_FLOAT)ri;
    const MD_FLOAT f0   = cs->tab_f[ri];
    const MD_FLOAT fc   = f0 + frac * (cs->tab_f[ri + 1] - f0);
    if (v != NULL) {
        *v = cs->tab_v[ri] - ((MD_FLOAT)0.5 / cs->tab_scale) * frac * (f0 + fc);
    }
    return fc;
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

// erf(z)/z of z2 = z^2 (GROMACS' pmePotentialCorrection())
static inline MD_HOST_DEVICE MD_FLOAT pmePotentialCorrectionScalar(MD_FLOAT z2)
{
    const MD_FLOAT z4 = z2 * z2;
#if PRECISION == 1
    MD_FLOAT polyVD1 = (MD_FLOAT)0.0066752224023576045451 * z4 +
                       (MD_FLOAT)0.43336185284710920150;
    MD_FLOAT polyVD0 = (MD_FLOAT)0.078647795836373922256 * z4 + (MD_FLOAT)1.0;
    polyVD0          = polyVD1 * z2 + polyVD0;

    MD_FLOAT polyVN0 = (MD_FLOAT)1.9296833005951166339e-8 * z4 +
                       (MD_FLOAT)0.000041603292906656984871;
    MD_FLOAT polyVN1 = (MD_FLOAT)-1.4213390571557850962e-6 * z4 +
                       (MD_FLOAT)-0.00013134036773265025626;
    polyVN0          = polyVN0 * z4 + (MD_FLOAT)0.038657983986041781264;
    polyVN1          = polyVN1 * z4 + (MD_FLOAT)0.11285044772717598220;
    polyVN0          = polyVN0 * z4 + (MD_FLOAT)1.1283802385263030286;
    polyVN0          = polyVN1 * z2 + polyVN0;
#else
    MD_FLOAT polyVD1 = (MD_FLOAT)0.000021784709867336150342 * z4 +
                       (MD_FLOAT)0.0096311444822588683504;
    MD_FLOAT polyVD0 = (MD_FLOAT)0.00064293662010911388448 * z4 +
                       (MD_FLOAT)0.085608012351550627051;
    polyVD1          = polyVD1 * z4 + (MD_FLOAT)0.43652499166614811084;
    polyVD0          = polyVD0 * z4 + (MD_FLOAT)1.0;
    polyVD0          = polyVD1 * z2 + polyVD0;

    MD_FLOAT polyVN1 = (MD_FLOAT)-9.3723776169321855475e-13 * z4 +
                       (MD_FLOAT)-7.3562157912251309487e-9;
    MD_FLOAT polyVN0 = (MD_FLOAT)1.2280156762674215741e-10 * z4 +
                       (MD_FLOAT)2.6215886208032517509e-7;
    polyVN1          = polyVN1 * z4 + (MD_FLOAT)-4.9532491651265819499e-6;
    polyVN0          = polyVN0 * z4 + (MD_FLOAT)0.00025907400778966060389;
    polyVN1          = polyVN1 * z4 + (MD_FLOAT)0.0010585044856156469792;
    polyVN0          = polyVN0 * z4 + (MD_FLOAT)0.045247661136833092885;
    polyVN1          = polyVN1 * z4 + (MD_FLOAT)0.11643931522926034421;
    polyVN0          = polyVN0 * z4 + (MD_FLOAT)1.1283791671726767970;
    polyVN0          = polyVN1 * z2 + polyVN0;
#endif
    return polyVN0 / polyVD0;
}

#endif
