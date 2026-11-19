/*
 * Minimal test runner helpers for MD-Bench unit tests.
 */
#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline void tr_log(const char* msg) { fprintf(stdout, "%s\n", msg); }

#define ASSERT_TRUE(cond, msg)                                                           \
    do {                                                                                 \
        if (!(cond)) {                                                                   \
            fprintf(stderr,                                                              \
                "ASSERT_TRUE failed: %s (at %s:%d)\n",                                   \
                msg,                                                                     \
                __FILE__,                                                                \
                __LINE__);                                                               \
            return 1;                                                                    \
        }                                                                                \
    } while (0)

#define ASSERT_INT_EQ(got, exp, msg)                                                     \
    do {                                                                                 \
        if ((got) != (exp)) {                                                            \
            fprintf(stderr,                                                              \
                "ASSERT_INT_EQ failed: %s (got=%d, exp=%d at %s:%d)\n",                  \
                msg,                                                                     \
                (int)(got),                                                              \
                (int)(exp),                                                              \
                __FILE__,                                                                \
                __LINE__);                                                               \
            return 1;                                                                    \
        }                                                                                \
    } while (0)

/* Values stored as MD_FLOAT cannot be compared tighter than their precision: in single
 * precision the tolerance is at least TEST_SP_RTOL relative to the expected value */
#if PRECISION == 1
#define TEST_SP_RTOL 1e-5
#else
#define TEST_SP_RTOL 0.0
#endif

#define ASSERT_NEAR(got, exp, tol, msg)                                                  \
    do {                                                                                 \
        double _g = (double)(got);                                                       \
        double _e = (double)(exp);                                                       \
        double _d = fabs(_g - _e);                                                       \
        double _t = fmax((double)(tol), TEST_SP_RTOL * fmax(1.0, fabs(_e)));             \
        if (_d > _t) {                                                                   \
            fprintf(stderr,                                                              \
                "ASSERT_NEAR failed: %s (got=%g, exp=%g, diff=%g, tol=%g at %s:%d)\n",   \
                msg,                                                                     \
                _g,                                                                      \
                _e,                                                                      \
                _d,                                                                      \
                _t,                                                                      \
                __FILE__,                                                                \
                __LINE__);                                                               \
            return 1;                                                                    \
        }                                                                                \
    } while (0)
