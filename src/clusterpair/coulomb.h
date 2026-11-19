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
 *   Ewald (tabulated):              correction = -Fcorr(r) / r, see ewald_table.h
 * The energies are those of coulomb_scalar.h. The Ewald corrections are the analytical
 * rational approximations used by GROMACS (src/gromacs/simd/.../scalar_math.h, LGPL).
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
    MD_SIMD_FLOAT k_rf; // reaction-field energy: k_rf r^2 - c_rf
    MD_SIMD_FLOAT c_rf;
    MD_SIMD_FLOAT beta;
    MD_SIMD_FLOAT beta2;
    MD_SIMD_FLOAT tab_scale;         // Ewald table points per length unit
    MD_SIMD_FLOAT tab_minus_half_dx; // -1 / (2 scale)
    MD_FLOAT* tab_f;                 // Ewald table, F (double precision layout)
    MD_FLOAT* tab_v;                 // V (double precision layout)
    MD_FLOAT* tab_fdv0;              // F, D, V, 0 (single precision layout)
} CoulombSimd;

static inline CoulombSimd coulombSimdInit(Parameter* param)
{
    CoulombSimd cs;
    const MD_FLOAT rc    = param->cutforce;
    const MD_FLOAT k_rf  = 0.5 / (rc * rc * rc);
    cs.epsfac            = simd_real_broadcast(param->epsfac);
    cs.k_rf2             = simd_real_broadcast(1.0 / (rc * rc * rc));
    cs.k_rf              = simd_real_broadcast(k_rf);
    cs.c_rf              = simd_real_broadcast(1.0 / rc + k_rf * rc * rc);
    cs.beta              = simd_real_broadcast(param->ewald_beta);
    cs.beta2             = simd_real_broadcast(param->ewald_beta * param->ewald_beta);
    cs.tab_scale         = simd_real_broadcast(param->ewald_table_scale);
    cs.tab_minus_half_dx = simd_real_broadcast(
        param->ewald_table_scale > 0 ? -0.5 / param->ewald_table_scale : 0.0);
    cs.tab_f    = param->ewald_table_f;
    cs.tab_v    = param->ewald_table_v;
    cs.tab_fdv0 = param->ewald_table_fdv0;
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

#if COULOMB_IS_EWALD
// erf(z)/z of z2 = z^2, SIMD copy of pmePotentialCorrectionScalar() (coulomb_scalar.h)
static inline MD_SIMD_FLOAT pmePotentialCorrectionSimd(MD_SIMD_FLOAT z2)
{
    const MD_SIMD_FLOAT z4 = simd_real_mul(z2, z2);
#define C(x) simd_real_broadcast((MD_FLOAT)(x))
#if PRECISION == 1
    MD_SIMD_FLOAT polyVD1 = simd_real_fma(C(0.0066752224023576045451),
        z4,
        C(0.43336185284710920150));
    MD_SIMD_FLOAT polyVD0 = simd_real_fma(C(0.078647795836373922256), z4, C(1.0));
    polyVD0               = simd_real_fma(polyVD1, z2, polyVD0);

    MD_SIMD_FLOAT polyVN0 = simd_real_fma(C(1.9296833005951166339e-8),
        z4,
        C(0.000041603292906656984871));
    MD_SIMD_FLOAT polyVN1 = simd_real_fma(C(-1.4213390571557850962e-6),
        z4,
        C(-0.00013134036773265025626));
    polyVN0               = simd_real_fma(polyVN0, z4, C(0.038657983986041781264));
    polyVN1               = simd_real_fma(polyVN1, z4, C(0.11285044772717598220));
    polyVN0               = simd_real_fma(polyVN0, z4, C(1.1283802385263030286));
    polyVN0               = simd_real_fma(polyVN1, z2, polyVN0);
#else
    MD_SIMD_FLOAT polyVD1 = simd_real_fma(C(0.000021784709867336150342),
        z4,
        C(0.0096311444822588683504));
    MD_SIMD_FLOAT polyVD0 = simd_real_fma(C(0.00064293662010911388448),
        z4,
        C(0.085608012351550627051));
    polyVD1               = simd_real_fma(polyVD1, z4, C(0.43652499166614811084));
    polyVD0               = simd_real_fma(polyVD0, z4, C(1.0));
    polyVD0               = simd_real_fma(polyVD1, z2, polyVD0);

    MD_SIMD_FLOAT polyVN1 = simd_real_fma(C(-9.3723776169321855475e-13),
        z4,
        C(-7.3562157912251309487e-9));
    MD_SIMD_FLOAT polyVN0 = simd_real_fma(C(1.2280156762674215741e-10),
        z4,
        C(2.6215886208032517509e-7));
    polyVN1               = simd_real_fma(polyVN1, z4, C(-4.9532491651265819499e-6));
    polyVN0               = simd_real_fma(polyVN0, z4, C(0.00025907400778966060389));
    polyVN1               = simd_real_fma(polyVN1, z4, C(0.0010585044856156469792));
    polyVN0               = simd_real_fma(polyVN0, z4, C(0.045247661136833092885));
    polyVN1               = simd_real_fma(polyVN1, z4, C(0.11643931522926034421));
    polyVN0               = simd_real_fma(polyVN0, z4, C(1.1283791671726767970));
    polyVN0               = simd_real_fma(polyVN1, z2, polyVN0);
#endif
#undef C
    return simd_real_mul(polyVN0, simd_real_reciprocal_accurate(polyVD0));
}
#endif /* COULOMB_IS_EWALD */

#if COULOMB_RULE == COULOMB_EWALD_TABLE
// Ewald table force correction at r (and potential correction into *v), see ewald_table.h
static inline MD_SIMD_FLOAT ewaldTableSimd(
    const CoulombSimd* cs, MD_SIMD_FLOAT r, MD_SIMD_FLOAT* v)
{
    const MD_SIMD_FLOAT rs   = simd_real_mul(r, cs->tab_scale);
    const MD_SIMD_INT ti     = simd_i32_from_real(rs);
    const MD_SIMD_FLOAT frac = simd_real_sub(rs, simd_real_from_i32(ti));
    MD_SIMD_FLOAT f0, df;
#if PRECISION == 1
    const MD_SIMD_INT ti2 = simd_i32_add(ti, ti);
    const MD_SIMD_INT ti4 = simd_i32_add(ti2, ti2);
    f0                    = simd_real_gather(ti4, cs->tab_fdv0, sizeof(MD_FLOAT));
    df                    = simd_real_gather(ti4, cs->tab_fdv0 + 1, sizeof(MD_FLOAT));
#else
    f0 = simd_real_gather(ti, cs->tab_f, sizeof(MD_FLOAT));
    df = simd_real_sub(simd_real_gather(ti, cs->tab_f + 1, sizeof(MD_FLOAT)), f0);
#endif
    const MD_SIMD_FLOAT fc = simd_real_fma(frac, df, f0);
    if (v != NULL) {
#if PRECISION == 1
        const MD_SIMD_FLOAT v0 = simd_real_gather(ti4,
            cs->tab_fdv0 + 2,
            sizeof(MD_FLOAT));
#else
        const MD_SIMD_FLOAT v0 = simd_real_gather(ti, cs->tab_v, sizeof(MD_FLOAT));
#endif
        *v = simd_real_fma(simd_real_mul(cs->tab_minus_half_dx, frac),
            simd_real_add(f0, fc),
            v0);
    }
    return fc;
}
#endif /* COULOMB_EWALD_TABLE */

#endif /* !USE_REFERENCE_KERNEL */

#endif /* COULOMB_RULE != COULOMB_NONE */

#endif
