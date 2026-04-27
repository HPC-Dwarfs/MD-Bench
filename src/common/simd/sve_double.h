/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <arm_acle.h>
#include <arm_sve.h>
#include <stdlib.h>

#define SIMD_INTRINSICS "sve_double"

#define MD_SIMD_FLOAT svfloat64_t
#define MD_SIMD_MASK  svbool_t
#define MD_SIMD_INT   svint64_t

static inline int simd_test_any(MD_SIMD_MASK a) { return svptest_any(svptrue_b64(), a); }
static inline MD_SIMD_FLOAT simd_real_broadcast(MD_FLOAT value)
{
    return svdup_f64(value);
}
static inline MD_SIMD_FLOAT simd_real_zero(void) { return svdup_f64(0.0); }
static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svsub_f64_z(svptrue_b64(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_load(const MD_FLOAT* ptr)
{
    return svld1_f64(svptrue_b64(), ptr);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    svint64_t offsets = svmul_n_s64_x(svptrue_b64(), vidx, scale);
    return svld1_gather_s64offset_f64(svptrue_b64(), base, offsets);
    // return svld1_gather_s64index_f64(svptrue_b64(), base, vidx);
}

static inline void simd_real_store(MD_FLOAT* ptr, MD_SIMD_FLOAT vec)
{
    svst1_f64(svptrue_b64(), ptr, vec);
}

static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svadd_f64_z(svptrue_b64(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svmul_f64_z(svptrue_b64(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return svmad_f64_z(svptrue_b64(), a, b, c);
}

static inline MD_SIMD_MASK simd_mask_from_u32(uint32_t a)
{
    return svdupq_n_b64(a & 0x1 ? 1 : 0, a & 0x2 ? 1 : 0);
}

static inline uint32_t simd_mask_to_u32(MD_SIMD_MASK mask)
{
    svuint64_t seq    = svindex_u64(0, 1);
    uint32_t result   = 0;
    MD_SIMD_MASK next = svpnext_b64(svptrue_b64(), mask);

    while (svptest_any(svptrue_b64(), next)) {
        result |= 1 << (uint32_t)svaddv_u64(next, seq);
        mask = svand_b_z(svptrue_b64(), mask, svnot_b_z(svptrue_b64(), next));
    }

    return result;
}

static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return svand_b_z(svptrue_b64(), a, b);
}

static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svcmplt_f64(svptrue_b64(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    MD_SIMD_FLOAT reciprocal = svrecpe_f64(a);
    reciprocal = svmul_f64_z(svptrue_b64(), reciprocal, svrecps_f64(reciprocal, a));
    return reciprocal;
}

static inline MD_FLOAT simd_real_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    svbool_t pg = svptrue_b64();
    double sum[4];
    sum[0] = svadda_f64(pg, 0.0, v0);
    sum[1] = svadda_f64(pg, 0.0, v1);
    sum[2] = svadda_f64(pg, 0.0, v2);
    sum[3] = svadda_f64(pg, 0.0, v3);
#if VECTOR_WIDTH >= 4
    pg             = SVE_DOUBLE4_MASK;
    svfloat64_t _m = svld1_f64(pg, m);
    svfloat64_t _s = svld1_f64(pg, sum);
    svst1_f64(pg, m, svadd_f64_x(pg, _m, _s));
    return svadda_f64(pg, 0.0, _s);
#else
    double res = 0;
    for (int i = 0; i < 4; i++) {
        m[i] += sum[i];
        res += sum[i];
    }
    return res;
#endif

    /*
    MD_SIMD_FLOAT sum0 = svaddp_f64_m(svptrue_b64(), v0, v1);
    MD_SIMD_FLOAT sum1 = svaddp_f64_m(svptrue_b64(), v2, v3);
    MD_SIMD_FLOAT odd  = svuzp2_f64(sum0, sum1);
    MD_SIMD_FLOAT even = svuzp1_f64(sum0, sum1);
    MD_SIMD_FLOAT sum  = svaddp_f64_m(svptrue_b64(), even, odd);

    MD_SIMD_FLOAT mem = svld1_f64(svptrue_b64(), m);
    sum               = svadd_f64_m(svptrue_b64(), sum, mem);

    svst1_f64(svptrue_b64(), m, sum);
    return svaddv_f64(svptrue_b64(), sum);
    */
}

static inline MD_SIMD_FLOAT simd_real_masked_add(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m)
{
    return svadd_f64_m(m, a, b);
}

static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK mask)
{
    return svsel_f64(mask, a, svdup_f64(0.0));
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_dual(): Not implemented for SVE with double precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_duplicate(): Not implemented for SVE with double precision!");
    exit(-1);
    return ret;
}

static inline void simd_real_h_decr3(
    MD_FLOAT* m, MD_SIMD_FLOAT a0, MD_SIMD_FLOAT a1, MD_SIMD_FLOAT a2)
{
    fprintf(stderr,
        "simd_real_h_decr3(): Not implemented for SVE with double precision!");
    exit(-1);
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    fprintf(stderr,
        "simd_real_h_dual_incr_reduced_sum(): Not implemented for SVE with double "
        "precision!");
    exit(-1);
    return 0.0f;
}

static inline MD_SIMD_INT simd_i32_broadcast(int a) { return svdup_s64((int64_t)a); }

static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svadd_s64_x(svptrue_b64(), a, b);
}

// Create sequence [0, 1, 2, ...] for SVE
static inline MD_SIMD_INT simd_i32_seq(void) { return svindex_s64(0, 1); }

// Integer multiply (for type indices)
static inline MD_SIMD_INT simd_i32_mul(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svmul_s64_x(svptrue_b64(), a, b);
}

// Integer comparison to mask
static inline MD_SIMD_MASK simd_mask_i32_cond_lt(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svcmplt_s64(svptrue_b64(), a, b);
}

// Masked integer load
static inline MD_SIMD_INT simd_i32_mask_load(const int* ptr, MD_SIMD_MASK mask)
{
    // Load 32-bit integers and extend to 64-bit
    svbool_t pg32        = svwhilelt_b32(0, VECTOR_WIDTH);
    svbool_t pg32_masked = svand_b_z(svptrue_b32(), pg32, svunpklo_b(mask));
    svint32_t loaded     = svld1_s32(pg32_masked, ptr);
    return svunpklo_s64(loaded);
}

// Gather integers (for multi-atom-type)
static inline MD_SIMD_INT simd_i32_gather(MD_SIMD_INT vidx, int* base, const int scale)
{
    // SVE gather for 32-bit integers with 64-bit indices
    svint64_t offsets  = svmul_n_s64_x(svptrue_b64(), vidx, sizeof(int));
    svint32_t gathered = svld1_gather_s64offset_s32(svptrue_b64(), base, offsets);
    // Sign-extend 32-bit to 64-bit
    return svunpklo_s64(gathered);
}

// Horizontal sum reduction
static inline MD_FLOAT simd_real_h_reduce_sum(MD_SIMD_FLOAT a)
{
    return svaddv_f64(svptrue_b64(), a);
}

// Masked scatter-subtract (for half-neighbor lists)
static inline void simd_real_masked_scatter_sub(
    MD_FLOAT* base, MD_SIMD_INT vidx, MD_SIMD_FLOAT v, MD_SIMD_MASK mask)
{
    // SVE doesn't have atomic scatter, use scalar fallback
    svbool_t pg = svptrue_b64();

    // Extract values and indices
    MD_FLOAT vals[VECTOR_WIDTH] __attribute__((aligned(64)));
    int64_t idx[VECTOR_WIDTH] __attribute__((aligned(64)));
    svst1_f64(pg, vals, v);
    svst1_s64(pg, idx, vidx);

    // Scalar atomic updates
    for (int i = 0; i < VECTOR_WIDTH; i++) {
        if (svptest_any(pg, mask)) {
#pragma omp atomic
            base[idx[i]] -= vals[i];
        }
        mask = svpnext_b64(pg, mask);
    }
}

static inline MD_SIMD_INT simd_i32_load(const int* m)
{
    svbool_t pg = svwhilelt_b32(0, VECTOR_WIDTH);
    return svunpklo_s64(svld1_s32(pg, m));
}

static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_duplicate(): Not implemented for SVE with double precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_dual_scaled(): Not implemented for SVE with double precision!");
    exit(-1);
    return ret;
}
