/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <immintrin.h>
#include <stdlib.h>
#include <string.h>

#define SIMD_INTRINSICS "avx2_double"

#define MD_SIMD_FLOAT __m256d
#define MD_SIMD_INT   __m128i
#define MD_SIMD_MASK  __m256d

static inline int simd_test_any(MD_SIMD_MASK a) { return _mm256_movemask_pd(a) != 0; }
static inline MD_SIMD_FLOAT simd_real_broadcast(MD_FLOAT scalar)
{
    return _mm256_set1_pd(scalar);
}
static inline MD_SIMD_FLOAT simd_real_zero(void) { return _mm256_set1_pd(0.0); }
static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_add_pd(a, b);
}
static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_sub_pd(a, b);
}
static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_mul_pd(a, b);
}
static inline MD_SIMD_FLOAT simd_real_load(MD_FLOAT* p) { return _mm256_load_pd(p); }
static inline void simd_real_store(MD_FLOAT* p, MD_SIMD_FLOAT a)
{
    _mm256_store_pd(p, a);
}
static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const MD_FLOAT* m)
{
    __m128d t0 = _mm_loadu_pd(m);
    __m256d t1 = _mm256_castpd128_pd256(t0);
    return _mm256_insertf128_pd(t1, t0, 1);
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const MD_FLOAT* m)
{
    __m128d t0 = _mm_load1_pd(m);
    __m128d t1 = _mm_load1_pd(m + 1);
    __m256d t3 = _mm256_castpd128_pd256(t0);
    return _mm256_insertf128_pd(t3, t1, 1);
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    __m256d t0, t1;
    t0 = _mm256_hadd_pd(v0, v0);
    t1 = _mm256_hadd_pd(v1, v1);
    t0 = _mm256_add_pd(t0, t1);

    __m256d mval = _mm256_load_pd(m);
    mval = _mm256_add_pd(mval, t0);
    _mm256_store_pd(m, mval);

    __m128d sum1 = _mm_add_pd(_mm256_castpd256_pd128(t0), _mm256_extractf128_pd(t0, 1));
    __m128d sum2 = _mm_hadd_pd(sum1, sum1);
    return _mm_cvtsd_f64(sum2);
}

static inline MD_FLOAT simd_real_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    __m256d t0, t1, t2;
    __m128d a0, a1;

    t0 = _mm256_hadd_pd(v0, v1);
    t1 = _mm256_hadd_pd(v2, v3);
    t2 = _mm256_permute2f128_pd(t0, t1, 0x21);
    t0 = _mm256_add_pd(t0, t2);
    t1 = _mm256_add_pd(t1, t2);
    t0 = _mm256_blend_pd(t0, t1, 0xC);
    t1 = _mm256_add_pd(t0, _mm256_load_pd(m));
    _mm256_store_pd(m, t1);

    t0 = _mm256_add_pd(t0, _mm256_permute_pd(t0, 0x5));
    a0 = _mm256_castpd256_pd128(t0);
    a1 = _mm256_extractf128_pd(t0, 0x1);
    a0 = _mm_add_sd(a0, a1);
    return *((MD_FLOAT*)&a0);
}

static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK m)
{
    return _mm256_and_pd(a, m);
}
static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    return _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(a)));
}
// static inline MD_SIMD_FLOAT simd_reciprocal(MD_SIMD_FLOAT a) { return
// _mm256_rcp14_pd(a); }
static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return _mm256_fmadd_pd(a, b, c);
}
static inline MD_SIMD_FLOAT simd_real_masked_add(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m)
{
    return simd_real_add(a, _mm256_and_pd(b, m));
}
static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_cmp_pd(a, b, _CMP_LT_OQ);
}
static inline MD_SIMD_MASK simd_mask_i32_cond_lt(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return _mm256_cvtepi32_pd(_mm_cmplt_epi32(a, b));
}
static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return _mm256_and_pd(a, b);
}
// TODO: Initialize all diagonal cases and just select the proper one (all bits set or
// diagonal) based on cond0
static inline MD_SIMD_MASK simd_mask_from_u32(unsigned int a)
{
    const unsigned long long int all  = 0xFFFFFFFFFFFFFFFF;
    const unsigned long long int none = 0x0;
    return _mm256_castsi256_pd(_mm256_set_epi64x((a & 0x8) ? all : none,
        (a & 0x4) ? all : none,
        (a & 0x2) ? all : none,
        (a & 0x1) ? all : none));
}
// TODO: Implement this, althrough it is just required for debugging
static inline unsigned int simd_mask_to_u32(MD_SIMD_MASK a) { return (unsigned int)_mm256_movemask_pd(a); }
static inline MD_FLOAT simd_real_h_reduce_sum(MD_SIMD_FLOAT a)
{
    __m128d a0, a1;
    // test with shuffle & add as an alternative to hadd later
    a  = _mm256_hadd_pd(a, a);
    a0 = _mm256_castpd256_pd128(a);
    a1 = _mm256_extractf128_pd(a, 0x1);
    a0 = _mm_add_sd(a0, a1);
    return *((MD_FLOAT*)&a0);
}

static inline void simd_h_decr(MD_FLOAT* m, MD_SIMD_FLOAT a)
{
    __m128d t0 = _mm256_castpd256_pd128(a);
    __m128d t1 = _mm256_extractf128_pd(a, 1);
    __m256d dup = _mm256_set_m128d(t1, t0);

    __m256d t = _mm256_load_pd(m);
    t = _mm256_sub_pd(t, dup);
    _mm256_store_pd(m, t);
}

static inline void simd_real_h_decr3(
    MD_FLOAT* m, MD_SIMD_FLOAT a0, MD_SIMD_FLOAT a1, MD_SIMD_FLOAT a2)
{
    simd_h_decr(m, a0);
    simd_h_decr(m + CLUSTER_N, a1);
    simd_h_decr(m + CLUSTER_N * 2, a2);
}

static inline MD_SIMD_INT simd_i32_broadcast(int scalar)
{
    return _mm_set1_epi32(scalar);
}
static inline MD_SIMD_INT simd_i32_zero(void) { return _mm_setzero_si128(); }
static inline MD_SIMD_INT simd_i32_seq(void) { return _mm_set_epi32(3, 2, 1, 0); }
static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return _mm_add_epi32(a, b);
}
static inline MD_SIMD_INT simd_i32_mul(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return _mm_mul_epi32(a, b);
}
static inline MD_SIMD_INT simd_i32_load(const int* m)
{
    return _mm_load_si128((__m128i const*)m);
}
static inline void simd_i32_store(int* m, MD_SIMD_INT a)
{
    _mm_store_si128((__m128i*)m, a);
}
static inline MD_SIMD_INT simd_i32_mask_load(const int* m, MD_SIMD_MASK k)
{
    __m128i imask = _mm256_cvtpd_epi32(k);
    return _mm_maskload_epi32(m, imask);
}

static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    return _mm_set_epi32(m[1], m[0], m[1], m[0]);
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    int i1 = m[0] * scale;
    int i2 = m[1] * scale;
    return _mm_set_epi32(i2, i2, i1, i1);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    if (scale == 1) {
        return _mm256_i32gather_pd(base, vidx, 1);
    } else if (scale == 2) {
        return _mm256_i32gather_pd(base, vidx, 2);
    } else if (scale == 4) {
        return _mm256_i32gather_pd(base, vidx, 4);
    } else {
        return _mm256_i32gather_pd(base, vidx, 8);
    }
}

static inline MD_SIMD_INT simd_i32_gather(
    MD_SIMD_INT vidx, int* base, const int scale)
{
    if (scale == 1) {
        return _mm256_i32gather_epi32(base, vidx, 1);
    } else if (scale == 2) {
        return _mm256_i32gather_epi32(base, vidx, 2);
    } else if (scale == 4) {
        return _mm256_i32gather_epi32(base, vidx, 4);
    } else {
        return _mm256_i32gather_epi32(base, vidx, 8);
    }
}
