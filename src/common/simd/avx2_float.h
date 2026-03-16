/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <immintrin.h>

#define SIMD_INTRINSICS "avx2_float"

#define MD_SIMD_FLOAT __m256
#define MD_SIMD_MASK  __m256
#define MD_SIMD_INT   __m256i

static inline int simd_test_any(MD_SIMD_MASK a) { return _mm256_movemask_ps(a) != 0; }
static inline MD_SIMD_FLOAT simd_real_broadcast(MD_FLOAT scalar)
{
    return _mm256_set1_ps(scalar);
}
static inline MD_SIMD_FLOAT simd_real_zero() { return _mm256_set1_ps(0.0); }
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
static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return _mm256_fmadd_ps(a, b, c);
}
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
    __m256i broadcast_mask = _mm256_set1_epi32(a);
    __m256i index          = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i shift          = _mm256_and_si256(broadcast_mask,
        _mm256_sllv_epi32(_mm256_set1_epi32(1), index));
    __m256i result         = _mm256_cmpgt_epi32(shift, _mm256_setzero_si256());
    return _mm256_castsi256_ps(result);
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
static inline void simd_i32_store(int* m, MD_SIMD_INT a)
{
    _mm256_store_si256((MD_SIMD_INT*)m, a);
}
static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    return _mm256_broadcastsi128_si256(_mm_load_si128((const __m128i*)(m)));
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    __m128i t0, t1;
    t0 = _mm_set1_epi32(m[0] * scale);
    t1 = _mm_set1_epi32(m[1] * scale);
    return _mm256_inserti128_si256(_mm256_castsi128_si256(t0), t1, 0x1);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    if (scale == 1) {
        return _mm256_i32gather_ps(base, vidx, 1);
    } else if (scale == 2) {
        return _mm256_i32gather_ps(base, vidx, 2);
    } else if (scale == 4) {
        return _mm256_i32gather_ps(base, vidx, 4);
    } else {
        return _mm256_i32gather_ps(base, vidx, 8);
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
// AVX2 has no hardware scatter; implement as scalar fallback
static inline void simd_real_masked_scatter_sub(
    MD_FLOAT* base, MD_SIMD_INT vidx, MD_SIMD_FLOAT v, MD_SIMD_MASK mask)
{
    unsigned int m = simd_mask_to_u32(mask);
    MD_FLOAT vals[8] __attribute__((aligned(32)));
    int idx[8] __attribute__((aligned(32)));
    simd_real_store(vals, v);
    simd_i32_store(idx, vidx);
    for (int i = 0; i < 8; i++) {
        if ((m >> i) & 1) {
            _Pragma("omp atomic")
            base[idx[i]] -= vals[i];
        }
    }
}

static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return _mm256_add_epi32(a, b);
}
