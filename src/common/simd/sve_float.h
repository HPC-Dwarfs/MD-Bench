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
    return svsub_f32_x(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_load(const float* ptr)
{
    return svld1_f32(svptrue_b32(), ptr);
}

// Every caller passes scale == sizeof(MD_FLOAT), so the index-scaled gather
// form applies directly (see sve_double.h's simd_real_gather for the same
// fix): it multiplies vidx by the destination element size (4 bytes for f32)
// as part of the gather's addressing mode, avoiding a separate
// svmul_n_s32_x to pre-compute byte offsets.
static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    (void)scale;
    return svld1_gather_s32index_f32(svptrue_b32(), base, vidx);
}


static inline void simd_real_store(MD_FLOAT* ptr, MD_SIMD_FLOAT vec)
{
    svst1_f32(svptrue_b32(), ptr, vec);
}

static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svadd_f32_x(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svmul_f32_x(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return svmad_f32_x(svptrue_b32(), a, b, c);
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
    svuint32_t seq  = svindex_u32(0, 1);
    uint32_t result = 0;
    svbool_t next   = svpfalse_b();

    while (svptest_any(mask, next = svpnext_b32(mask, next))) {
        result |= 1u << svaddv_u32(next, seq);
    }

    return result;
}

static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return svand_b_z(svptrue_b32(), a, b);
}

static inline MD_SIMD_MASK simd_mask_not(MD_SIMD_MASK a)
{
    return svnot_b_z(svptrue_b32(), a);
}

static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return svcmplt_f32(svptrue_b32(), a, b);
}

static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    // svrecpe_f32(a) + one svrecps_f32 Newton-Raphson step was measured slower here (gather-bound kernel, div hides in memory stalls better than the dependency chain):
    // MD_SIMD_FLOAT reciprocal = svrecpe_f32(a);
    // reciprocal = svmul_f32_x(svptrue_b32(), reciprocal, svrecps_f32(reciprocal, a));
    // return reciprocal;
    return svdiv_f32_x(svptrue_b32(), svdup_f32(1.0f), a);
}

static inline float simd_real_incr_reduced_sum(
    float* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    svbool_t pg = svptrue_b32();
    float sum[4];
    sum[0] = svadda_f32(pg, 0.0f, v0);
    sum[1] = svadda_f32(pg, 0.0f, v1);
    sum[2] = svadda_f32(pg, 0.0f, v2);
    sum[3] = svadda_f32(pg, 0.0f, v3);
#if VECTOR_WIDTH >= 4
    pg                 = svptrue_pat_b32(SV_VL4);
    svfloat32_t _m     = svld1_f32(pg, m);
    svfloat32_t _s     = svld1_f32(pg, sum);
    svst1_f32(pg, m, svadd_f32_x(pg, _m, _s));
    return svadda_f32(pg, 0.0f, _s);
#else
    float res = 0;
    for (int i = 0; i < 4; i++) {
        m[i] += sum[i];
        res += sum[i];
    }
    return res;
#endif
}

static inline float simd_real_incr_reduced_sum_j2(
     float* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    svbool_t pg = svptrue_b32();
    float sum[2];
    sum[0] = svadda_f32(pg, 0.0f, v0);
    sum[1] = svadda_f32(pg, 0.0f, v1);
#if VECTOR_WIDTH >= 2
    pg                 = svptrue_pat_b32(SV_VL2);
    svfloat32_t _m     = svld1_f32(pg, m);
    svfloat32_t _s     = svld1_f32(pg, sum);
    svst1_f32(pg, m, svadd_f32_x(pg, _m, _s));
    return svadda_f32(pg, 0.0f, _s);
#else
    float res = 0;
    for (int i = 0; i < 2; i++) {
        m[i] += sum[i];
        res += sum[i];
    }
    return res;
#endif
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

static inline MD_SIMD_INT simd_i32_mul(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svmul_s32_x(svptrue_b32(), a, b);
}

static inline MD_SIMD_INT simd_i32_seq(void) { return svindex_s32(0, 1); }

// Single whilelt instruction -- the idiomatic SVE way to predicate a loop's remaining trip count, used by the VLA kernels.
static inline MD_SIMD_MASK simd_mask_from_remaining(int remaining) { return svwhilelt_b32(0, remaining); }

static inline MD_SIMD_MASK simd_mask_i32_cond_lt(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svcmplt_s32(svptrue_b32(), a, b);
}

static inline MD_SIMD_MASK simd_mask_i32_cond_eq(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svcmpeq_s32(svptrue_b32(), a, b);
}

// SVE natively supports predicated loads: inactive lanes are zeroed, matching
// the semantics buildNeighborCPU()'s SIMD path relies on for its tail mask.
static inline MD_SIMD_INT simd_i32_mask_load(const int* ptr, MD_SIMD_MASK mask)
{
    return svld1_s32(mask, ptr);
}

static inline MD_SIMD_INT simd_i32_gather(MD_SIMD_INT vidx, int* base, const int scale)
{
    svint32_t offsets = svmul_n_s32_x(svptrue_b32(), vidx, scale);
    return svld1_gather_s32offset_s32(svptrue_b32(), base, offsets);
}

// Horizontal sum reduction
static inline MD_FLOAT simd_real_h_reduce_sum(MD_SIMD_FLOAT a)
{
    return svaddv_f32(svptrue_b32(), a);
}

// Masked scatter-subtract (for half-neighbor lists). SVE has no atomic
// scatter, so fall back to a scalar loop
static inline void simd_real_masked_scatter_sub(
    MD_FLOAT* base, MD_SIMD_INT vidx, MD_SIMD_FLOAT v, MD_SIMD_MASK mask)
{
    unsigned int m = simd_mask_to_u32(mask);
    MD_FLOAT vals[VECTOR_WIDTH] __attribute__((aligned(64)));
    int32_t idx[VECTOR_WIDTH] __attribute__((aligned(64)));
    svst1_f32(svptrue_b32(), vals, v);
    svst1_s32(svptrue_b32(), idx, vidx);

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        if ((m >> i) & 1) {
#pragma omp atomic
            base[idx[i]] -= vals[i];
        }
    }
}

static inline MD_SIMD_INT simd_i32_load(const int* m)
{
    return svld1_s32(svptrue_b32(), m);
}

// SVE's svld1 has no alignment requirement, so unaligned is the same load.
static inline MD_SIMD_INT simd_i32_loadu(const int* m)
{
    return svld1_s32(svptrue_b32(), m);
}

static inline void simd_i32_store(int* m, MD_SIMD_INT a)
{
    svst1_s32(svptrue_b32(), m, a);
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
static inline MD_SIMD_FLOAT simd_real_sqrt(MD_SIMD_FLOAT v)
{
    return svsqrt_f32_x(svptrue_b32(), v);
}
static inline MD_SIMD_INT simd_i32_from_real(MD_SIMD_FLOAT v)
{
    return svcvt_s32_f32_x(svptrue_b32(), v);
}
static inline MD_SIMD_FLOAT simd_real_from_i32(MD_SIMD_INT v)
{
    return svcvt_f32_s32_x(svptrue_b32(), v);
}
static inline MD_SIMD_INT simd_i32_min(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return svmin_s32_x(svptrue_b32(), a, b);
}
