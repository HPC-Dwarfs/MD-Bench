/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Electrostatics for the cluster-pair force kernels (see force_lj_4xn.c).
 *
 * The formulation follows GROMACS' NBNxM kernels: a pair contributes the force
 *   f_scalar = qq * (interact * 1/r^3 + <correction>)      (f = f_scalar * (r_i - r_j))
 * where "interact" is 0 for excluded pairs, so excluded pairs inside the cut-off still
 * receive the correction term:
 *   reaction-field (eps_rf = inf):  correction = -2 k_rf,   k_rf = 1/(2 rc^3)
 *   Ewald (real space):             correction = beta^3 * pmeForceCorrection(beta^2 r^2)
 * The Ewald force correction is the analytical rational approximation used by GROMACS
 * (src/gromacs/simd/.../scalar_math.h, LGPL).
 */
#ifndef __COULOMB_H_
#define __COULOMB_H_

#include <math.h>

#include <parameter.h>

#include <kernel_math.h>

#ifndef USE_REFERENCE_KERNEL
typedef struct {
    MD_SIMD_FLOAT epsfac;
    MD_SIMD_FLOAT k_rf2;
    MD_SIMD_FLOAT beta;
    MD_SIMD_FLOAT beta2;
} CoulombSimd;

static inline CoulombSimd coulombSimdInit(Parameter* param)
{
    CoulombSimd cs;
    cs.epsfac = simd_real_broadcast(param->epsfac);
    cs.k_rf2  = simd_real_broadcast(
        1.0 / (param->cutforce * param->cutforce * param->cutforce));
    cs.beta  = simd_real_broadcast(param->ewald_beta);
    cs.beta2 = simd_real_broadcast(param->ewald_beta * param->ewald_beta);
    return cs;
}
#endif /* !USE_REFERENCE_KERNEL */

#if COULOMB_RULE != COULOMB_NONE

#ifndef USE_REFERENCE_KERNEL
#if COULOMB_RULE == COULOMB_EWALD
static inline MD_SIMD_FLOAT pmeForceCorrectionSimd(MD_SIMD_FLOAT z2)
{
    const MD_SIMD_FLOAT z4 = simd_real_mul(z2, z2);
#define C(x) simd_real_broadcast((MD_FLOAT)(x))
#if PRECISION == 1
    MD_SIMD_FLOAT polyFD0 = simd_real_fma(C(0.0011193462567257629232),
        z4,
        C(0.11583842382862377919));
    MD_SIMD_FLOAT polyFD1 = simd_real_fma(C(0.014866955030185295499),
        z4,
        C(0.50736591960530292870));
    polyFD0               = simd_real_fma(polyFD0, z4, C(1.0));
    polyFD0               = simd_real_fma(polyFD1, z2, polyFD0);

    MD_SIMD_FLOAT polyFN0 = simd_real_fma(C(-1.7357322914161492954e-8),
        z4,
        C(-0.000053401640219807709149));
    MD_SIMD_FLOAT polyFN1 = simd_real_fma(C(1.4703624142580877519e-6),
        z4,
        C(0.0010054721316683106153));
    polyFN0               = simd_real_fma(polyFN0, z4, C(-0.019278317264888380590));
    polyFN1               = simd_real_fma(polyFN1, z4, C(0.069670166153766424023));
    polyFN0               = simd_real_fma(polyFN0, z4, C(-0.75225204789749321333));
    polyFN0               = simd_real_fma(polyFN1, z2, polyFN0);
#else
    MD_SIMD_FLOAT polyFD1 = simd_real_fma(C(0.000016009278224355026701),
        z4,
        C(0.0081803507497974289008));
    polyFD1               = simd_real_fma(polyFD1, z4, C(0.41543303143712535988));
    polyFD1               = simd_real_mul(polyFD1, z2);
    MD_SIMD_FLOAT polyFD0 = simd_real_fma(C(0.00051055686934806966046),
        z4,
        C(0.077181146026670287235));
    polyFD0               = simd_real_fma(polyFD0, z4, C(1.0));
    polyFD0               = simd_real_add(polyFD0, polyFD1);

    MD_SIMD_FLOAT polyFN0 = simd_real_fma(C(-8.0072854618360083154e-14),
        z4,
        C(-8.1490406329798423616e-10));
    polyFN0               = simd_real_fma(polyFN0, z4, C(-9.9471420832602741006e-7));
    polyFN0               = simd_real_fma(polyFN0, z4, C(-0.00031991745139313364005));
    polyFN0               = simd_real_fma(polyFN0, z4, C(-0.031750380176100813405));
    polyFN0               = simd_real_fma(polyFN0, z4, C(-0.75225277815249618847));
    MD_SIMD_FLOAT polyFN1 = simd_real_fma(C(1.1859116242260148027e-11),
        z4,
        C(3.4404793543907847655e-8));
    polyFN1               = simd_real_fma(polyFN1, z4, C(0.000020740315999115847456));
    polyFN1               = simd_real_fma(polyFN1, z4, C(0.0035074449373659008203));
    polyFN1               = simd_real_fma(polyFN1, z4, C(0.13884101728898463426));
    polyFN0               = simd_real_fma(polyFN1, z2, polyFN0);
#endif
#undef C
    return simd_real_mul(polyFN0, simd_real_reciprocal_accurate(polyFD0));
}
#endif /* COULOMB_EWALD */

#endif /* !USE_REFERENCE_KERNEL */

#endif /* COULOMB_RULE != COULOMB_NONE */

#endif
