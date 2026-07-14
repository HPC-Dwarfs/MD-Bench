#include <arm_acle.h>
#include <arm_neon.h>
#include <stdlib.h>

#define SIMD_INTRINSICS "neon_float"

#define MD_SIMD_FLOAT float32x4_t
#define MD_SIMD_MASK  uint32x4_t
#define MD_SIMD_INT   int32x4_t

static inline int simd_test_any(MD_SIMD_MASK a) { return vmaxvq_u32(a) != 0; }
static inline MD_SIMD_FLOAT simd_real_broadcast(float value)
{
    return vdupq_n_f32(value);
}
static inline MD_SIMD_FLOAT simd_real_zero(void) { return vdupq_n_f32(0.0f); }
static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vsubq_f32(a, b);
}
static inline MD_SIMD_FLOAT simd_real_load(const float* ptr) { return vld1q_f32(ptr); }
static inline void simd_real_store(MD_FLOAT* ptr, MD_SIMD_FLOAT vec)
{
    vst1q_f32(ptr, vec);
}
static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vaddq_f32(a, b);
}
static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vmulq_f32(a, b);
}
static inline MD_SIMD_FLOAT simd_real_fma(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c)
{
    return vfmaq_f32(c, a, b);
}

static inline MD_SIMD_MASK simd_mask_from_u32(uint32_t a)
{
    const uint32_t all  = 0xFFFFFFFF;
    const uint32_t none = 0x0;
    return vsetq_lane_u32((a & 0x8) ? all : none,
        vsetq_lane_u32((a & 0x4) ? all : none,
            vsetq_lane_u32((a & 0x2) ? all : none,
                vsetq_lane_u32((a & 0x1) ? all : none, vdupq_n_u32(0), 0),
                1),
            2),
        3);
}

static inline uint32_t simd_mask_to_u32(MD_SIMD_MASK mask) { return 0; }

static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b)
{
    return vandq_u32(a, b);
}

static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b)
{
    return vcltq_f32(a, b);
}

static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a)
{
    MD_SIMD_FLOAT reciprocal = vrecpeq_f32(a);
    reciprocal               = vmulq_f32(reciprocal, vrecpsq_f32(reciprocal, a));
    return reciprocal;
}

static inline float simd_real_incr_reduced_sum(
    float* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3)
{
    float32x4_t sum0 = vpaddq_f32(v0, v1);
    float32x4_t sum1 = vpaddq_f32(v2, v3);
    float32x4_t sum  = vpaddq_f32(sum0, sum1);

    float32x4_t mem = vld1q_f32(m);
    sum             = vaddq_f32(sum, mem);
    vst1q_f32(m, sum);

    float32x2_t sum_low  = vget_low_f32(sum);
    float32x2_t sum_high = vget_high_f32(sum);
    sum_low              = vadd_f32(sum_low, sum_high);
    sum_low              = vpadd_f32(sum_low, sum_low);
    return vget_lane_f32(sum_low, 0);
}

static inline float simd_real_incr_reduced_sum_j2(
     float* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
 {
    float32x4_t sum0 = vpaddq_f32(v0, v1); 
    float32x2_t sum = vpadd_f32(vget_low_f32(sum0), vget_high_f32(sum0)); 
 
     float32x2_t mem = vld1_f32(m); 
     sum             = vadd_f32(sum, mem); 
     vst1_f32(m, sum);  
     
     return vget_lane_f32(vpadd_f32(sum, sum), 0); 
 }

static inline MD_SIMD_FLOAT simd_real_masked_add(
    MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m)
{
    MD_SIMD_FLOAT masked_b = vreinterpretq_f32_u32(
        vandq_u32(vreinterpretq_u32_f32(b), m));
    return vaddq_f32(a, masked_b);
}

static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK mask)
{
    return vbslq_f32(mask, a, vdupq_n_f32(0.0f));
}

static inline MD_SIMD_INT simd_i32_load(const int* ptr) { return vld1q_s32(ptr); }
static inline MD_SIMD_INT simd_i32_broadcast(int value) { return vdupq_n_s32(value); }
static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b)
{
    return vaddq_s32(a, b);
}

static inline MD_SIMD_FLOAT simd_real_gather(
    MD_SIMD_INT vidx, MD_FLOAT* base, const int scale)
{
#if defined(__ARM_FEATURE_MVE)
    return vldrwq_gather_offset_f32(base, vreinterpretq_u32_s32(vidx));
#else
    MD_SIMD_FLOAT result = vdupq_n_f32(0);

    result = vld1q_lane_f32(&base[vgetq_lane_s32(vidx, 0)], result, 0);
    result = vld1q_lane_f32(&base[vgetq_lane_s32(vidx, 1)], result, 1);
    result = vld1q_lane_f32(&base[vgetq_lane_s32(vidx, 2)], result, 2);
    result = vld1q_lane_f32(&base[vgetq_lane_s32(vidx, 3)], result, 3);
    return result;
#endif
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_dual(): Not implemented for NEON with single precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const MD_FLOAT* m)
{
    MD_SIMD_FLOAT ret;
    fprintf(stderr,
        "simd_real_load_h_duplicate(): Not implemented for NEON with single precision!");
    exit(-1);
    return ret;
}

static inline void simd_real_h_decr3(
    MD_FLOAT* m, MD_SIMD_FLOAT a0, MD_SIMD_FLOAT a1, MD_SIMD_FLOAT a2)
{
    fprintf(stderr,
        "simd_real_h_decr3(): Not implemented for NEON with single precision!");
    exit(-1);
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1)
{
    fprintf(stderr,
        "simd_real_h_dual_incr_reduced_sum(): Not implemented for NEON with single "
        "precision!");
    exit(-1);
    return 0.0f;
}

static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int* m)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_duplicate(): Not implemented for NEON with single precision!");
    exit(-1);
    return ret;
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int* m, int scale)
{
    MD_SIMD_INT ret;
    fprintf(stderr,
        "simd_i32_load_h_dual_scaled(): Not implemented for NEON with single precision!");
    exit(-1);
    return ret;
}
