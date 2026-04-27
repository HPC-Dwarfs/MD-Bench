/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <immintrin.h>
#include <stdlib.h>
#include <string.h>
#ifndef NO_ZMM_INTRIN
#include <zmmintrin.h>
#endif

#define SIMD_INTRINSICS "avx512_float"

#define MD_SIMD_FLOAT   __m512
#define MD_SIMD_MASK    __mmask16
#define MD_SIMD_IBOOL   __mmask16
#define MD_SIMD_INT     __m512i
#define MD_SIMD_BITMASK __m512i

static inline int simd_test_any(MD_SIMD_MASK a) { return _mm512_kortestz(a, a) == 0; }
static inline MD_SIMD_BITMASK simd_load_bitmask(const int* m)
{
    return _mm512_load_si512(m);
}

static inline MD_SIMD_INT simd_i32_broadcast(int a) { return _mm512_set1_epi32(a); }

static inline MD_SIMD_IBOOL simd_test_bits(MD_SIMD_FLOAT a)
{
    return _mm512_test_epi32_mask(_mm512_castps_si512(a), _mm512_castps_si512(a));
}

static inline MD_SIMD_MASK cvtIB2B(MD_SIMD_IBOOL a) { return a; }
static inline MD_SIMD_FLOAT simd_real_broadcast(float scalar)
{
    return _mm512_set1_ps(scalar);
}
static inline MD_SIMD_FLOAT simd_real_zero(void) { return _mm512_set1_ps(0.0f); }
static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm512_add_ps(a, b);
}
static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm512_sub_ps(a, b);
}
static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm512_mul_ps(a, b);
}
static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return _mm512_fmadd_ps(a, b, c);
}
static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    return _mm512_rcp14_ps(a);
}
static inline MD_SIMD_FLOAT simd_real_masked_add(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m)
{
    return _mm512_mask_add_ps(a, m, a, b);
}
static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return _kand_mask16(a, b);
}
static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return _mm512_cmp_ps_mask(a, b, _CMP_LT_OQ);
}
static inline MD_SIMD_MASK simd_mask_from_u32(unsigned int a)
{
    return _cvtu32_mask16(a);
}
static inline unsigned int simd_mask_to_u32(MD_SIMD_MASK a) { return _cvtmask16_u32(a); }
static inline MD_SIMD_FLOAT simd_real_load(MD_FLOAT* p) { return _mm512_load_ps(p); }
static inline void simd_real_store(MD_FLOAT* p, MD_SIMD_FLOAT a)
{
    _mm512_store_ps(p, a);
}
static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK m)
{
    return _mm512_mask_mov_ps(_mm512_setzero_ps(), m, a);
}
static inline MD_FLOAT simd_real_h_reduce_sum(MD_SIMD_FLOAT a)
{
    // This would only be called in a Mx16 configuration, which is not valid in GROMACS
    fprintf(stderr,
        "simd_h_reduce_sum(): Called with AVX512 intrinsics and single-precision which "
        "is not valid!\n");
    exit(-1);
    return 0.0;
}

static inline MD_FLOAT simd_real_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    // This would only be called in a Mx16 configuration, which is not valid in GROMACS
    fprintf(stderr,
        "simd_h_reduce_sum(): Called with AVX512 intrinsics and single-precision which "
        "is not valid!\n");
    exit(-1);
    return 0.0;
}

static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const float* m)
{
    return _mm512_castpd_ps(_mm512_broadcast_f64x4(_mm256_load_pd((const double*)(m))));
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const float* m)
{
    return _mm512_shuffle_f32x4(_mm512_broadcastss_ps(_mm_load_ss(m)),
        _mm512_broadcastss_ps(_mm_load_ss(m + 1)),
        0x44);
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    float* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    __m512 t0, t1;
    __m128 t2, t3;

    t0 = _mm512_shuffle_f32x4(v0, v1, 0x88);
    t1 = _mm512_shuffle_f32x4(v0, v1, 0xdd);
    t0 = _mm512_add_ps(t0, t1);
    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0x4e));
    t0 = _mm512_add_ps(t0, _mm512_permute_ps(t0, 0xb1));
    t0 = _mm512_maskz_compress_ps(simd_mask_from_u32(0x1111ul), t0);
    t3 = _mm512_castps512_ps128(t0);
    t2 = _mm_load_ps(m);
    t2 = _mm_add_ps(t2, t3);
    _mm_store_ps(m, t2);

    t3 = _mm_add_ps(t3, _mm_permute_ps(t3, 0x4e));
    t3 = _mm_add_ps(t3, _mm_permute_ps(t3, 0xb1));
    return _mm_cvtss_f32(t3);
}

static inline void simd_h_decr(MD_FLOAT* m, MD_SIMD_FLOAT a)
{
    __m256 t;
    a = _mm512_add_ps(a, _mm512_shuffle_f32x4(a, a, 0xee));
    t = _mm256_load_ps(m);
    t = _mm256_sub_ps(t, _mm512_castps512_ps256(a));
    _mm256_store_ps(m, t);
}

static inline void simd_real_h_decr3(
    MD_FLOAT* m, MD_SIMD_FLOAT a0, MD_SIMD_FLOAT a1, MD_SIMD_FLOAT a2)
{
    simd_h_decr(m, a0);
    simd_h_decr(m + CLUSTER_N, a1);
    simd_h_decr(m + CLUSTER_N * 2, a2);
}

static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return _mm512_add_epi32(a, b);
}

static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    return _mm512_broadcast_i32x8(_mm256_load_si256((const __m256i*)m));
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    return _mm512_inserti32x8(_mm512_set1_epi32(m[0] * scale),
        _mm256_set1_epi32(m[1] * scale),
        1);
}

static inline MD_SIMD_INT simd_i32_load(const int* m)
{
    return _mm512_load_si512((const MD_SIMD_INT*)m);
}

static inline void simd_i32_store(int* m, MD_SIMD_INT a)
{
    _mm512_store_si512((MD_SIMD_INT*)m, a);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    if (scale == 1) {
        return _mm512_i32gather_ps(vidx, base, 1);
    } else if (scale == 2) {
        return _mm512_i32gather_ps(vidx, base, 2);
    } else if (scale == 4) {
        return _mm512_i32gather_ps(vidx, base, 4);
    } else {
        return _mm512_i32gather_ps(vidx, base, 8);
    }
}

static inline MD_SIMD_INT simd_i32_gather(MD_SIMD_INT vidx, int* base, const int scale)
{
    if (scale == 1) {
        return _mm512_i32gather_epi32(vidx, base, 1);
    } else if (scale == 2) {
        return _mm512_i32gather_epi32(vidx, base, 2);
    } else if (scale == 4) {
        return _mm512_i32gather_epi32(vidx, base, 4);
    } else {
        return _mm512_i32gather_epi32(vidx, base, 8);
    }
}
static inline void simd_real_masked_scatter_sub(
    MD_FLOAT* base, MD_SIMD_INT vidx, MD_SIMD_FLOAT v, MD_SIMD_MASK mask)
{
    MD_FLOAT vals[16] __attribute__((aligned(64)));
    int idx[16] __attribute__((aligned(64)));
    simd_real_store(vals, v);
    simd_i32_store(idx, vidx);
    for (int i = 0; i < 16; i++) {
        if ((mask >> i) & 1) {
            _Pragma("omp atomic") base[idx[i]] -= vals[i];
        }
    }
}
