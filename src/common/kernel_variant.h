/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Variants of the CPU kernels (force kernels and neighbor-list builds), selected at
 * runtime without a branch in the kernel. A kernel file is a template on the macros
 * below; the Makefile compiles it once per variant (VARIANTS_<file>),
 * RULE_VARIANT_NAME/VARIANT_NAME give the functions of each variant their own name and
 * the SELECT_* and CALL_RULE_VARIANT macros pick one. LJ_COMB_RULE:        LJ combination
 * rule (--lj-comb), every copy has one ACCURATE_RECIPROCAL: Newton-Raphson refined 1/r
 * and 1/r^2 (--accurate-recip) KERNEL_NO_PREFETCH:  no j-cluster prefetch
 * (--no-prefetch), cluster-pair only KERNEL_EXCLUSIONS:   molecular exclusions (systems
 * with molecules) A runtime flag inside the kernel instead costs 1.5-3%, as GCC then lays
 * out the kernel code differently. */
#ifndef __KERNEL_VARIANT_H_
#define __KERNEL_VARIANT_H_

#include <parameter.h>

#ifdef LJ_COMB_RULE
#if LJ_COMB_RULE == LJ_COMB_SINGLE
#define KV_RULE Single
#elif LJ_COMB_RULE == LJ_COMB_GEOM
#define KV_RULE Geom
#elif LJ_COMB_RULE == LJ_COMB_MATRIX
#define KV_RULE Matrix
#else
#error "LJ_COMB_RULE must be LJ_COMB_SINGLE, LJ_COMB_GEOM or LJ_COMB_MATRIX"
#endif
#endif

#ifdef ACCURATE_RECIPROCAL
#define KV_ACCURATE Accurate
#else
#define KV_ACCURATE
#endif
#ifdef KERNEL_NO_PREFETCH
#define KV_NO_PREFETCH NoPrefetch
#else
#define KV_NO_PREFETCH
#endif
#ifdef KERNEL_EXCLUSIONS
#define KV_EXCL     Excl
#define KERNEL_EXCL 1
#else
#define KV_EXCL
#define KERNEL_EXCL 0
#endif

#define KV_PASTE(name, r, a, p, x)  name##r##a##p##x
#define KV_EXPAND(name, r, a, p, x) KV_PASTE(name, r, a, p, x)
// Functions that only depend on the combination rule
#define RULE_VARIANT_NAME(name) KV_EXPAND(name, KV_RULE, , , )
// Kernel entry points: the rule and the options of the copy
#define VARIANT_NAME(name) KV_EXPAND(name, KV_RULE, KV_ACCURATE, KV_NO_PREFETCH, KV_EXCL)

// Selection of a copy by the runtime options in param (name: the base name)
#define KV_SEL_RULE(param, name, SEL)                                                    \
    ((param)->lj_comb_rule == LJ_COMB_SINGLE                                             \
            ? SEL(param, name##Single)                                                   \
            : ((param)->lj_comb_rule == LJ_COMB_GEOM ? SEL(param, name##Geom)            \
                                                     : SEL(param, name##Matrix)))
#define KV_SEL_NONE(param, name) name
#define KV_SEL_EXCL(param, name) ((param)->molsize > 1 ? name##Excl : name)
#define KV_SEL_ACCURATE(param, name)                                                     \
    ((param)->accurate_reciprocal ? name##Accurate : name)
#define KV_SEL_OPTIONS(param, name)                                                      \
    ((param)->accurate_reciprocal                                                        \
            ? ((param)->kernel_prefetch ? KV_SEL_EXCL(param, name##Accurate)             \
                                        : KV_SEL_EXCL(param, name##AccurateNoPrefetch))  \
            : ((param)->kernel_prefetch ? KV_SEL_EXCL(param, name)                       \
                                        : KV_SEL_EXCL(param, name##NoPrefetch)))

// Rule only (reference, table and EAM kernels)
#define SELECT_RULE_VARIANT(param, name) KV_SEL_RULE(param, name, KV_SEL_NONE)
// Rule and exclusions (verlet-list C kernels)
#define SELECT_RULE_EXCL_VARIANT(param, name) KV_SEL_RULE(param, name, KV_SEL_EXCL)
// Rule and --accurate-recip (verlet-list SIMD kernels)
#define SELECT_RULE_ACCURATE_VARIANT(param, name)                                        \
    KV_SEL_RULE(param, name, KV_SEL_ACCURATE)
// Rule, --accurate-recip, --no-prefetch and exclusions (cluster-pair SIMD kernels)
#define SELECT_SIMD_VARIANT(param, name) KV_SEL_RULE(param, name, KV_SEL_OPTIONS)

// Declarations of the copies
#define KV_DECL(name) extern double name(Parameter*, Atom*, Neighbor*, Stats*);
#define KV_DECL_RULES(name, D)                                                           \
    D(name##Single)                                                                      \
    D(name##Geom)                                                                        \
    D(name##Matrix)
#define KV_DECL_EXCL(name)     KV_DECL(name) KV_DECL(name##Excl)
#define KV_DECL_ACCURATE(name) KV_DECL(name) KV_DECL(name##Accurate)
#define KV_DECL_OPTIONS(name)                                                            \
    KV_DECL_EXCL(name)                                                                   \
    KV_DECL_EXCL(name##Accurate)                                                         \
    KV_DECL_EXCL(name##NoPrefetch)                                                       \
    KV_DECL_EXCL(name##AccurateNoPrefetch)
#define DECLARE_RULE_VARIANTS(name)          KV_DECL_RULES(name, KV_DECL)
#define DECLARE_RULE_EXCL_VARIANTS(name)     KV_DECL_RULES(name, KV_DECL_EXCL)
#define DECLARE_RULE_ACCURATE_VARIANTS(name) KV_DECL_RULES(name, KV_DECL_ACCURATE)
#define DECLARE_SIMD_VARIANTS(name)          KV_DECL_RULES(name, KV_DECL_OPTIONS)

// Functions with other signatures that only depend on the rule (neighbor-list builds)
#define DECLARE_RULE_VARIANT_FUNCTIONS(ret, name, ...)                                   \
    ret name##Single(__VA_ARGS__);                                                       \
    ret name##Geom(__VA_ARGS__);                                                         \
    ret name##Matrix(__VA_ARGS__);
#define CALL_RULE_VARIANT(rule, name, ...)                                               \
    ((rule) == LJ_COMB_SINGLE ? name##Single(__VA_ARGS__)                                \
                              : ((rule) == LJ_COMB_GEOM ? name##Geom(__VA_ARGS__)        \
                                                        : name##Matrix(__VA_ARGS__)))

#endif
