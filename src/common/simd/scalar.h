/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdbool.h>
#include <string.h>

#define SIMD_INTRINSICS "scalar"

// At least for x86, this works
typedef MD_FLOAT SimdRealRegT
    __attribute__((__vector_size__(VECTOR_WIDTH * sizeof(MD_FLOAT)), __may_alias__));
typedef int SimdIntRegT
    __attribute__((__vector_size__(VECTOR_WIDTH * sizeof(MD_FLOAT)), __may_alias__));

typedef union {
    MD_FLOAT val[VECTOR_WIDTH];
    SimdRealRegT reg;
} MD_SIMD_FLOAT;

typedef union {
    int val[VECTOR_WIDTH];
    SimdIntRegT reg;
} MD_SIMD_INT;

typedef unsigned short MD_SIMD_MASK;

static inline int simd_test_any(MD_SIMD_MASK a) { return a != 0; }

static inline MD_SIMD_FLOAT simd_real_broadcast(MD_FLOAT scalar) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = scalar;
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_zero(void) { return simd_real_broadcast(0.0); }
static inline MD_SIMD_FLOAT simd_real_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = a.val[i] + b.val[i];
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_sub(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = a.val[i] - b.val[i];
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_mul(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = a.val[i] * b.val[i];
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_fma(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_FLOAT c) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = a.val[i] * b.val[i] + c.val[i];
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_reciprocal(MD_SIMD_FLOAT a) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = 1.0 / a.val[i];
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_masked_add(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b, MD_SIMD_MASK m) {
    MD_SIMD_FLOAT result = a;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        if (m & (1 << i)) {
            result.val[i] = a.val[i] + b.val[i];
        }
    }

    return result;
}

static inline MD_SIMD_MASK simd_mask_and(MD_SIMD_MASK a, MD_SIMD_MASK b) { return a & b; }
static inline MD_SIMD_MASK simd_mask_cond_lt(MD_SIMD_FLOAT a, MD_SIMD_FLOAT b) {
    MD_SIMD_MASK result = 0;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        if (a.val[i] < b.val[i]) {
            result |= (1 << i);
        }
    }

    return result;
}

static inline MD_SIMD_MASK simd_mask_from_u32(unsigned int a) { return (MD_SIMD_MASK)a; }
static inline unsigned int simd_mask_to_u32(MD_SIMD_MASK a) { return (unsigned int)a; }

static inline MD_SIMD_FLOAT simd_real_load(MD_FLOAT* p) {
    MD_SIMD_FLOAT result;
    memcpy(result.val, p, sizeof(result.val));
    return result;
}

static inline void simd_real_store(MD_FLOAT* p, MD_SIMD_FLOAT a) { memcpy(p, a.val, sizeof(a.val)); }
static inline MD_SIMD_FLOAT simd_real_select_by_mask(MD_SIMD_FLOAT a, MD_SIMD_MASK m) {
    MD_SIMD_FLOAT result = simd_real_zero();

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        if (m & (1 << i)) {
            result.val[i] = a.val[i];
        }
    }

    return result;
}

static inline MD_FLOAT simd_real_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1, MD_SIMD_FLOAT v2, MD_SIMD_FLOAT v3) {
    MD_FLOAT sum = 0;
    MD_FLOAT partial_sums[4];

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        partial_sums[0] += v0.val[i];
        partial_sums[1] += v1.val[i];
        partial_sums[2] += v2.val[i];
        partial_sums[3] += v3.val[i];
    }

    for (int i = 0; i < 4; i++) {
        m[i] += partial_sums[i];
        sum += partial_sums[i];
    }

    return sum;
}

static inline MD_SIMD_FLOAT simd_real_load_h_duplicate(const MD_FLOAT* m) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH / 2; i++) {
        result.val[i]                      = m[i];
        result.val[i + (VECTOR_WIDTH / 2)] = result.val[i];
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_load_h_dual(const MD_FLOAT* m) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH / 2; i++) {
        result.val[i]                      = m[0];
        result.val[i + (VECTOR_WIDTH / 2)] = m[1];
    }

    return result;
}

static inline MD_FLOAT simd_real_h_dual_incr_reduced_sum(
    MD_FLOAT* m, MD_SIMD_FLOAT v0, MD_SIMD_FLOAT v1) {
    MD_FLOAT sum = (MD_FLOAT)(0.0);
    MD_FLOAT partial_sums[4] = {0};

    for (int i = 0; i < VECTOR_WIDTH / 2; i++) {
        partial_sums[0] += v0.val[i];
        partial_sums[1] += v0.val[i + (VECTOR_WIDTH / 2)];
        partial_sums[2] += v1.val[i];
        partial_sums[3] += v1.val[i + (VECTOR_WIDTH / 2)];
    }

    for (int i = 0; i < 4; i++) {
        m[i] += partial_sums[i];
        sum += partial_sums[i];
    }

    return sum;
}

static inline void simd_real_h_decr(MD_FLOAT* m, MD_SIMD_FLOAT a) {
    for (int i = 0; i < VECTOR_WIDTH / 2; i++) {
        m[i] -= a.val[i] + a.val[(VECTOR_WIDTH / 2) + i];
    }
}

static inline void simd_real_h_decr3(
    MD_FLOAT* m, MD_SIMD_FLOAT a0, MD_SIMD_FLOAT a1, MD_SIMD_FLOAT a2) {
    simd_real_h_decr(m, a0);
    simd_real_h_decr(m + CLUSTER_N, a1);
    simd_real_h_decr(m + CLUSTER_N * 2, a2);
}

static inline MD_SIMD_INT simd_i32_broadcast(int scalar) {
    MD_SIMD_INT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = scalar;
    }

    return result;
}

static inline MD_SIMD_INT simd_i32_zero(void) { return simd_i32_broadcast(0); }

static inline MD_SIMD_INT simd_i32_seq(void) {
    MD_SIMD_INT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = i;
    }

    return result;
}

static inline MD_SIMD_INT simd_i32_load(const int* m) {
    MD_SIMD_INT result;
    memcpy(result.val, m, sizeof(result.val));
    return result;
}

static inline MD_SIMD_INT simd_i32_add(MD_SIMD_INT a, MD_SIMD_INT b) {
    MD_SIMD_INT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = a.val[i] + b.val[i];
    }

    return result;
}

static inline MD_SIMD_INT simd_i32_mul(MD_SIMD_INT a, MD_SIMD_INT b) {
    MD_SIMD_INT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = a.val[i] * b.val[i];
    }

    return result;
}

static inline MD_SIMD_INT simd_i32_mask_load(const int* m, MD_SIMD_MASK k) {
    MD_SIMD_INT result = simd_i32_zero();

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        if (k & (1 << i)) {
            result.val[i] = m[i];
        }
    }

    return result;
}

static inline MD_SIMD_MASK simd_mask_int_cond_lt(MD_SIMD_INT a, MD_SIMD_INT b) {
    MD_SIMD_MASK result = 0;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        if (a.val[i] < b.val[i]) {
            result |= (1 << i);
        }
    }

    return result;
}

static inline MD_SIMD_INT simd_i32_load_h_duplicate(const int *m) {
    MD_SIMD_INT result;

    for (int i = 0; i < VECTOR_WIDTH / 2; i++) {
        result.val[i] = m[i];
    }

    for (int i = VECTOR_WIDTH / 2; i < VECTOR_WIDTH; i++) {
        result.val[i] = m[i - (VECTOR_WIDTH / 2)];
    }

    return result;
}

static inline MD_SIMD_INT simd_i32_load_h_dual_scaled(const int *m, int scale) {
    MD_SIMD_INT result;

    for (int i = 0; i < VECTOR_WIDTH / 2; i++) {
        result.val[i] = m[0] * scale;
    }

    for (int i = VECTOR_WIDTH / 2; i < VECTOR_WIDTH; i++) {
        result.val[i] = m[1] * scale;
    }

    return result;
}

static inline MD_SIMD_FLOAT simd_real_gather(MD_SIMD_INT vidx, MD_FLOAT *base, const int scale) {
    MD_SIMD_FLOAT result;

    for (int i = 0; i < VECTOR_WIDTH; i++) {
        result.val[i] = base[vidx.val[i]];
    }

    return result;
}
