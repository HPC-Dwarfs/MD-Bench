/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Reciprocal and reciprocal square root helpers of the cluster-pair SIMD kernels.
 *
 * The hardware estimates (simd_real_reciprocal and the rsqrt estimate used here) are only
 * 12-14 bits accurate. By default the kernels use them as they are, which is the cheapest
 * option; with ACCURATE_RECIPROCAL they are refined with Newton-Raphson steps to the
 * working precision, which makes the forces comparable with codes that do the same
 * (GROMACS iterates the estimates to full precision).
 */
#ifndef __KERNEL_MATH_H_
#define __KERNEL_MATH_H_

#include <parameter.h>

#ifndef USE_REFERENCE_KERNEL
// force.h defines the cluster dimensions that simd.h would otherwise default
#include <force.h>
#include <simd.h>

static inline MD_SIMD_FLOAT simd_real_reciprocal_accurate(MD_SIMD_FLOAT a)
{
    const MD_SIMD_FLOAT two = simd_real_broadcast(2.0);
    MD_SIMD_FLOAT x         = simd_real_reciprocal(a);
    x                       = simd_real_mul(x, simd_real_sub(two, simd_real_mul(a, x)));
#if PRECISION == 2
    x = simd_real_mul(x, simd_real_sub(two, simd_real_mul(a, x)));
#endif
    return x;
}

/* Estimate of 1/sqrt(a). Targets without a hardware estimate use the exact sqrt and
 * reciprocal. */
#if defined(__ISA_AVX512__) && PRECISION == 1
#define SIMD_RSQRT_ESTIMATE_BITS 14
static inline MD_SIMD_FLOAT simd_real_rsqrt_estimate(MD_SIMD_FLOAT a)
{
    return _mm512_rsqrt14_ps(a);
}
#elif defined(__ISA_AVX512__) && PRECISION == 2
#define SIMD_RSQRT_ESTIMATE_BITS 14
static inline MD_SIMD_FLOAT simd_real_rsqrt_estimate(MD_SIMD_FLOAT a)
{
    return _mm512_rsqrt14_pd(a);
}
#elif (defined(__ISA_AVX2__) || defined(__ISA_AVX__)) && PRECISION == 1
#define SIMD_RSQRT_ESTIMATE_BITS 12
static inline MD_SIMD_FLOAT simd_real_rsqrt_estimate(MD_SIMD_FLOAT a)
{
    return _mm256_rsqrt_ps(a);
}
#elif (defined(__ISA_AVX2__) || defined(__ISA_AVX__)) && PRECISION == 2
#define SIMD_RSQRT_ESTIMATE_BITS 12
static inline MD_SIMD_FLOAT simd_real_rsqrt_estimate(MD_SIMD_FLOAT a)
{
    return _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(a)));
}
#else
#define SIMD_RSQRT_ESTIMATE_BITS 0
#endif

#if SIMD_RSQRT_ESTIMATE_BITS > 0
/* One Newton-Raphson step doubles the number of correct bits: x * (3 - a x^2) / 2 */
static inline MD_SIMD_FLOAT simd_real_rsqrt_step(MD_SIMD_FLOAT a, MD_SIMD_FLOAT x)
{
    const MD_SIMD_FLOAT half       = simd_real_broadcast(0.5);
    const MD_SIMD_FLOAT three      = simd_real_broadcast(3.0);
    const MD_SIMD_FLOAT correction = simd_real_sub(three,
        simd_real_mul(simd_real_mul(a, x), x));
    return simd_real_mul(simd_real_mul(half, x), correction);
}

static inline MD_SIMD_FLOAT simd_real_rsqrt_accurate(MD_SIMD_FLOAT a)
{
    MD_SIMD_FLOAT x = simd_real_rsqrt_estimate(a);
    x               = simd_real_rsqrt_step(a, x);
#if PRECISION == 2
    // 12 or 14 bits -> 24 or 28 -> 48 or 56 bits
    x = simd_real_rsqrt_step(a, x);
#if SIMD_RSQRT_ESTIMATE_BITS < 14
    x = simd_real_rsqrt_step(a, x);
#endif
#endif
    return x;
}
#define SIMD_RSQRT_ESTIMATE simd_real_rsqrt_estimate
#else
static inline MD_SIMD_FLOAT simd_real_rsqrt_accurate(MD_SIMD_FLOAT a)
{
    return simd_real_reciprocal_accurate(simd_real_sqrt(a));
}
static inline MD_SIMD_FLOAT simd_real_rsqrt_noaccurate(MD_SIMD_FLOAT a)
{
    return simd_real_reciprocal(simd_real_sqrt(a));
}
#define SIMD_RSQRT_ESTIMATE simd_real_rsqrt_noaccurate
#endif

#ifdef ACCURATE_RECIPROCAL
#define SIMD_RECIPROCAL simd_real_reciprocal_accurate
#define SIMD_RSQRT      simd_real_rsqrt_accurate
#else
#define SIMD_RECIPROCAL simd_real_reciprocal
#define SIMD_RSQRT      SIMD_RSQRT_ESTIMATE
#endif
#endif /* !USE_REFERENCE_KERNEL */

#endif
