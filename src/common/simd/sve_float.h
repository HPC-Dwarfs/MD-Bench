/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <arm_acle.h>
#include <arm_sve.h>
#include <stdlib.h>

#define SIMD_INTRINSICS "sve_float"

#define MD_SIMD_FLOAT svfloat32_t
#define MD_SIMD_MASK  svbool_t
#define MD_SIMD_INT   svint32_t

static inline int simd_test_any(MD_SIMD_MASK a) { return svptest_any(svptrue_b32(), a); }
static inline MD_SIMD_FLOAT simd_real_broadcast(float value) { return svdup_f32(value); }
static inline MD_SIMD_FLOAT simd_real_zero(void) { return svdup_f32(0.0f); }
static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svsub_f32_z(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_load(const float* ptr)
{
    return svld1_f32(svptrue_b32(), ptr);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    svint32_t offsets = svmul_n_s32_x(svptrue_b32(), vidx, scale);
    return svld1_gather_s32offset_f32(svptrue_b32(), base, offsets);
    //return svld1_gather_s32index_f32(svptrue_b32(), base, vidx);
}

static inline void simd_real_store(MD_FLOAT* ptr, MD_SIMD_FLOAT vec)
{
    svst1_f32(svptrue_b32(), ptr, vec);
}

static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svadd_f32_z(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svmul_f32_z(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return svmad_f32_z(svptrue_b32(), a, b, c);
}

static inline MD_SIMD_MASK simd_mask_from_u32(uint32_t a)
{
    svbool_t predicate = svdupq_n_b32(a & 0x1 ? 1 : 0,
        a & 0x2 ? 1 : 0,
        a & 0x4 ? 1 : 0,
        a & 0x8 ? 1 : 0);

    return predicate;
}

static inline uint32_t simd_mask_to_u32(MD_SIMD_MASK mask)
{
    svuint32_t seq    = svindex_u32(0, 1);
    uint32_t result   = 0;
    MD_SIMD_MASK next = svpnext_b32(svptrue_b32(), mask);

    while (svptest_any(svptrue_b32(), next)) {
        result |= 1 << svaddv_u32(next, seq);
        mask = svand_b_z(svptrue_b32(), mask, svnot_b_z(svptrue_b32(), next));
    }

    return result;
}

static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return svand_b_z(svptrue_b32(), a, b);
}

static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svcmplt_f32(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    MD_SIMD_FLOAT reciprocal = svrecpe_f32(a);
    reciprocal = svmul_f32_z(svptrue_b32(), reciprocal, svrecps_f32(reciprocal, a));
    return reciprocal;
}

static inline float simd_real_incr_reduced_sum(
    float* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    /*
    svbool_t    pg = svptrue_b32();
    svfloat32_t _m, _s;
    float32_t   sum[4];
    sum[0] = svadda_f32(pg, 0.0f, v0);
    sum[1] = svadda_f32(pg, 0.0f, v1);
    sum[2] = svadda_f32(pg, 0.0f, v2);
    sum[3] = svadda_f32(pg, 0.0f, v3);
    pg     = svptrue_pat_b32(SV_VL4);
    _m     = svld1_f32(pg, m);
    _s     = svld1_f32(pg, sum);
    svst1_f32(pg, m, svadd_f32_x(pg, _m, _s));
    return svadda_f32(pg, 0.0f, _s);
    */

    MD_SIMD_FLOAT sum0 = svaddp_f32_m(svptrue_b32(), v0, v1);
    MD_SIMD_FLOAT sum1 = svaddp_f32_m(svptrue_b32(), v2, v3);
    MD_SIMD_FLOAT odd  = svuzp2_f32(sum0, sum1);
    MD_SIMD_FLOAT even = svuzp1_f32(sum0, sum1);
    MD_SIMD_FLOAT sum  = svaddp_f32_m(svptrue_b32(), even, odd);

    MD_SIMD_FLOAT mem = svld1_f32(svptrue_b32(), m);
    sum               = svadd_f32_m(svptrue_b32(), sum, mem);

    svst1_f32(svptrue_b32(), m, sum);
    return svaddv_f32(svptrue_b32(), sum);
}

static inline MD_SIMD_FLOAT simd_real_masked_add(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m)
{
    return svadd_f32_m(m, a, b);
}

static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK mask)
{
    return svsel_f32(mask, a, svdup_f32(0.0f));
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_dual(): Not implemented for SVE with single precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_duplicate(): Not implemented for SVE with single precision!");
    exit(-1);
    return ret;
}

static inline void simd_real_h_decr3(
    MD_FLOAT* m, MD_SIMD_FLOAT a0, MD_SIMD_FLOAT a1, MD_SIMD_FLOAT a2)
{
    fprintf(stderr,
        "simd_real_h_decr3(): Not implemented for SVE with single precision!");
    exit(-1);
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    fprintf(stderr,
        "simd_real_h_dual_incr_reduced_sum(): Not implemented for SVE with single "
        "precision!");
    exit(-1);
    return 0.0f;
}

static inline MD_SIMD_INT simd_i32_broadcast(int a) { return svdup_s32(a); }

static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svadd_s32_x(svptrue_b32(), a, b);
}

static inline MD_SIMD_INT simd_i32_load(const int* m)
{
    return svld1_s32(svptrue_b32(), m);
}

static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_duplicate(): Not implemented for SVE with single precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_dual_scaled(): Not implemented for SVE with single precision!");
    exit(-1);
    return ret;
}
