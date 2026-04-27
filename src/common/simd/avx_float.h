/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <immintrin.h>

#define SIMD_INTRINSICS "avx_float"

#define MD_SIMD_FLOAT __m256
#define MD_SIMD_MASK  __m256
#define MD_SIMD_INT   __m256i

static inline int simd_test_any(MD_SIMD_MASK a) { return _mm256_movemask_ps(a) != 0; }
static inline MD_SIMD_FLOAT simd_real_broadcast(MD_FLOAT scalar)
{
    return _mm256_set1_ps(scalar);
}
static inline MD_SIMD_FLOAT simd_real_zero(void) { return _mm256_set1_ps(0.0); }
static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_add_ps(a, b);
}
static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_sub_ps(a, b);
}
static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_mul_ps(a, b);
}
static inline MD_SIMD_FLOAT simd_real_load(MD_FLOAT* p) { return _mm256_load_ps(p); }
static inline void simd_real_store(MD_FLOAT* p, MD_SIMD_FLOAT a)
{
    _mm256_store_ps(p, a);
}
static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK m)
{
    return _mm256_and_ps(a, m);
}
static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    return _mm256_rcp_ps(a);
}

#ifdef __ISA_AVX_FMA__
static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return _mm256_fmadd_ps(a, b, c);
}
#else
static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return simd_real_add(simd_real_mul(a, b), c);
}
#endif

static inline MD_SIMD_FLOAT simd_real_masked_add(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m)
{
    return _mm256_add_ps(a, _mm256_and_ps(b, m));
}
static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm256_cmp_ps(a, b, _CMP_LT_OQ);
}
static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return _mm256_and_ps(a, b);
}

static inline MD_SIMD_MASK simd_mask_from_u32(unsigned int a)
{
    const unsigned long int all  = 0xFFFFFFFF;
    const unsigned long int none = 0x0;
    return _mm256_castsi256_ps(_mm256_set_epi32((a & 0x80) ? all : none,
        (a & 0x40) ? all : none,
        (a & 0x20) ? all : none,
        (a & 0x10) ? all : none,
        (a & 0x8) ? all : none,
        (a & 0x4) ? all : none,
        (a & 0x2) ? all : none,
        (a & 0x1) ? all : none));
}

static inline unsigned int simd_mask_to_u32(MD_SIMD_MASK a)
{
    return _mm256_movemask_ps(a);
}
static inline MD_FLOAT simd_real_h_reduce_sum(MD_SIMD_FLOAT a)
{
    __m128 t0;
    t0 = _mm_add_ps(_mm256_castps256_ps128(a), _mm256_extractf128_ps(a, 0x1));
    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *((MD_FLOAT*)&t0);
}

static inline MD_FLOAT simd_real_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    __m128 t0, t2;
    v0 = _mm256_hadd_ps(v0, v1);
    v2 = _mm256_hadd_ps(v2, v3);
    v0 = _mm256_hadd_ps(v0, v2);
    t0 = _mm_add_ps(_mm256_castps256_ps128(v0), _mm256_extractf128_ps(v0, 0x1));
    t2 = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t2);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *((MD_FLOAT*)&t0);
}

static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const MD_FLOAT* m)
{
    return _mm256_broadcast_ps((const __m128*)(m));
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const MD_FLOAT* m)
{
    __m128 t0, t1;
    t0 = _mm_broadcast_ss(m);
    t1 = _mm_broadcast_ss(m + 1);
    return _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t1, 0x1);
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    __m128 t0, t1;
    v0 = _mm256_hadd_ps(v0, v1);
    t0 = _mm256_extractf128_ps(v0, 0x1);
    t0 = _mm_hadd_ps(_mm256_castps256_ps128(v0), t0);
    t0 = _mm_permute_ps(t0, _MM_SHUFFLE(3, 1, 2, 0));
    t1 = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t1);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *((MD_FLOAT*)&t0);
}

inline void simd_h_decr(MD_FLOAT* m, MD_SIMD_FLOAT a)
{
    __m128 asum = _mm_add_ps(_mm256_castps256_ps128(a), _mm256_extractf128_ps(a, 0x1));
    _mm_store_ps(m, _mm_sub_ps(_mm_load_ps(m), asum));
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
    return _mm256_set1_epi32(scalar);
}
static inline MD_SIMD_INT simd_i32_load(const int* m)
{
    return _mm256_load_si256((MD_SIMD_INT*)m);
}
static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    __m128i val = _mm_load_si128((const __m128i*)(m));
    __m256i res = _mm256_castsi128_si256(val);
    return _mm256_insertf128_si256(res, val, 1);
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    __m128i t0     = _mm_set1_epi32(m[0] * scale);
    __m128i t1     = _mm_set1_epi32(m[1] * scale);
    __m256i result = _mm256_castsi128_si256(t0);
    return _mm256_insertf128_si256(result, t1, 0x1);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    __m128i vidx_lo   = _mm256_castsi256_si128(vidx);
    __m128i vidx_hi   = _mm256_extractf128_si256(vidx, 1);
    __m128i scaled_lo = vidx_lo; // _mm_mullo_epi32(vidx_lo, _mm_set1_epi32(scale));
    __m128i scaled_hi = vidx_hi; // _mm_mullo_epi32(vidx_hi, _mm_set1_epi32(scale));

    int i0 = _mm_extract_epi32(scaled_lo, 0);
    int i1 = _mm_extract_epi32(scaled_lo, 1);
    int i2 = _mm_extract_epi32(scaled_lo, 2);
    int i3 = _mm_extract_epi32(scaled_lo, 3);
    int i4 = _mm_extract_epi32(scaled_hi, 0);
    int i5 = _mm_extract_epi32(scaled_hi, 1);
    int i6 = _mm_extract_epi32(scaled_hi, 2);
    int i7 = _mm_extract_epi32(scaled_hi, 3);

    __m128 gath_lo = _mm_set_ps(base[i3], base[i2], base[i1], base[i0]);
    __m128 gath_hi = _mm_set_ps(base[i7], base[i6], base[i5], base[i4]);
    return _mm256_insertf128_ps(_mm256_castps128_ps256(gath_lo), gath_hi, 1);
}

static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    __m128i low_add  = _mm_add_epi32(_mm256_extractf128_si256(a, 0),
        _mm256_extractf128_si256(b, 0));
    __m128i high_add = _mm_add_epi32(_mm256_extractf128_si256(a, 1),
        _mm256_extractf128_si256(b, 1));

    return _mm256_set_m128i(high_add, low_add);
}

static inline MD_SIMD_INT simd_i32_zero(void) { return _mm256_setzero_si256(); }
static inline MD_SIMD_INT simd_i32_seq(void)
{
    return _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
}
static inline MD_SIMD_INT simd_i32_mul(MD_SIMD_INT a, MD_SIMD_INT b)
{
    __m128i low_mul  = _mm_mullo_epi32(_mm256_extractf128_si256(a, 0),
        _mm256_extractf128_si256(b, 0));
    __m128i high_mul = _mm_mullo_epi32(_mm256_extractf128_si256(a, 1),
        _mm256_extractf128_si256(b, 1));
    return _mm256_set_m128i(high_mul, low_mul);
}
static inline MD_SIMD_INT simd_i32_mask_load(const int* m, MD_SIMD_MASK k)
{
    // AVX doesn't have masked load; emulate with scalar fallback
    unsigned int u32_mask                      = simd_mask_to_u32(k);
    int result[8] __attribute__((aligned(32))) = { 0 };
    for (int i = 0; i < 8; i++) {
        if ((u32_mask >> i) & 1) {
            result[i] = m[i];
        }
    }
    return _mm256_load_si256((const __m256i*)result);
}
static inline MD_SIMD_MASK simd_mask_i32_cond_lt(MD_SIMD_INT a, MD_SIMD_INT b)
{
    __m128i low_cmp  = _mm_cmplt_epi32(_mm256_extractf128_si256(a, 0),
        _mm256_extractf128_si256(b, 0));
    __m128i high_cmp = _mm_cmplt_epi32(_mm256_extractf128_si256(a, 1),
        _mm256_extractf128_si256(b, 1));
    __m256i imask    = _mm256_set_m128i(high_cmp, low_cmp);
    return _mm256_castsi256_ps(imask);
}
