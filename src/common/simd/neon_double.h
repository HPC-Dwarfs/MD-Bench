#include <arm_acle.h>
#include <arm_neon.h>
#include <stdlib.h>

#define SIMD_INTRINSICS "neon_double"

#define MD_SIMD_FLOAT float64x2_t
#define MD_SIMD_MASK  uint64x2_t
#define MD_SIMD_INT   int64x2_t

static inline int simd_test_any(MD_SIMD_MASK a)
{
    return vgetq_lane_u64(a, 0) != 0 || vgetq_lane_u64(a, 1) != 0;
}
static inline MD_SIMD_FLOAT simd_real_broadcast(MD_FLOAT value)
{
    return vdupq_n_f64(value);
}
static inline MD_SIMD_FLOAT simd_real_zero(void) { return vdupq_n_f64(0.0f); }
static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vsubq_f64(a, b);
}
static inline MD_SIMD_FLOAT simd_real_load(const MD_FLOAT* ptr) { return vld1q_f64(ptr); }
static inline void simd_real_store(MD_FLOAT* ptr, MD_SIMD_FLOAT vec)
{
    vst1q_f64(ptr, vec);
}
static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vaddq_f64(a, b);
}
static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vmulq_f64(a, b);
}
static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return vfmaq_f64(c, a, b);
}

static inline MD_SIMD_MASK simd_mask_from_u32(uint32_t a)
{
    const uint64_t all  = 0xFFFFFFFFFFFFFFFFULL;
    const uint64_t none = 0x0;
    MD_SIMD_MASK result;

    result = vsetq_lane_u64((a & 0x1) ? all : none, result, 0);
    result = vsetq_lane_u64((a & 0x2) ? all : none, result, 1);
    return result;
}

static inline uint32_t simd_mask_to_u32(MD_SIMD_MASK mask) { return 0; }

static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return vandq_u64(a, b);
}

static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vcltq_f64(a, b);
}

static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    MD_SIMD_FLOAT reciprocal = vrecpeq_f64(a);
    reciprocal               = vmulq_f64(reciprocal, vrecpsq_f64(reciprocal, a));
    return reciprocal;
}

static inline MD_FLOAT simd_real_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    float64x2_t t1, t2, t3, t4;

    t1 = vuzp1q_f64(v0, v1);
    t2 = vuzp2q_f64(v0, v1);
    t3 = vuzp1q_f64(v2, v3);
    t4 = vuzp2q_f64(v2, v3);

    t1 = vaddq_f64(t1, t2);
    t3 = vaddq_f64(t3, t4);

    t2 = vaddq_f64(t1, vld1q_f64(m));
    t4 = vaddq_f64(t3, vld1q_f64(m + 2));
    vst1q_f64(m, t2);
    vst1q_f64(m + 2, t4);

    t1 = vaddq_f64(t1, t3);
    t2 = vpaddq_f64(t1, t1);

    return vgetq_lane_f64(t2, 0);

    /*
    float64x2_t sum0 = vpaddq_f64(v0, v1);
    float64x2_t sum1 = vpaddq_f64(v2, v3);
    float64x2_t sum  = vpaddq_f64(sum0, sum1);

    float64x2_t mem = vld1q_f64(m);
    sum             = vaddq_f64(sum, mem);
    vst1q_f64(m, sum);
    return vget_lane_f64(vget_low_f64(sum) + vget_high_f64(sum), 0);
    */
}

static inline MD_FLOAT simd_real_incr_reduced_sum_j2(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    float64x2_t t1, t2;
 
    t1 = vuzp1q_f64(v0, v1);
    t2 = vuzp2q_f64(v0, v1);
     
    t1 = vaddq_f64(t1, t2);
 
    t2 = vaddq_f64(t1, vld1q_f64(m));
    vst1q_f64(m, t2);
 
    t2 = vpaddq_f64(t1, t1);
 
    return vgetq_lane_f64(t2, 0);
}

static inline MD_SIMD_FLOAT simd_real_masked_add(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m)
{
    MD_SIMD_FLOAT masked_b = vreinterpretq_f64_u64(
        vandq_u64(vreinterpretq_u64_f64(b), m));
    return vaddq_f64(a, masked_b);
}

static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK mask)
{
    return vbslq_f64(mask, a, vdupq_n_f64(0.0f));
}

static inline MD_SIMD_INT simd_i32_load(const int* ptr)
{
    return vmovl_s32(vld1_s32(ptr));
}
static inline MD_SIMD_INT simd_i32_broadcast(int value) { return vdupq_n_s64(value); }
static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return vaddq_s64(a, b);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
    MD_SIMD_FLOAT result = vdupq_n_f64(0);
    result               = vld1q_lane_f64(&base[vgetq_lane_s64(vidx, 0)], result, 0);
    result               = vld1q_lane_f64(&base[vgetq_lane_s64(vidx, 1)], result, 1);
    return result;
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_dual(): Not implemented for NEON with double precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_duplicate(): Not implemented for NEON with double precision!");
    exit(-1);
    return ret;
}

static inline void simd_real_h_decr3(
    MD_FLOAT* m, MD_SIMD_FLOAT a0, MD_SIMD_FLOAT a1, MD_SIMD_FLOAT a2)
{
    fprintf(stderr,
        "simd_real_h_decr3(): Not implemented for NEON with double precision!");
    exit(-1);
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    fprintf(stderr,
        "simd_real_h_dual_incr_reduced_sum(): Not implemented for NEON with double "
        "precision!");
    exit(-1);
    return 0.0f;
}

static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_duplicate(): Not implemented for NEON with double precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_dual_scaled(): Not implemented for NEON with double precision!");
    exit(-1);
    return ret;
}
