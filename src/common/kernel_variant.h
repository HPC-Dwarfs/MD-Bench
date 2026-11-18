/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Compiled copies of the CPU SIMD force kernels, selected at runtime without a branch in
 * the kernel. The Makefile compiles each kernel file once more for every combination of
 * the macros below (KERNEL_COPIES), KERNEL_FN gives the entry points of each copy their
 * own name and initForce() picks one with SELECT_KERNEL.
 *   ACCURATE_RECIPROCAL: Newton-Raphson refined 1/r and 1/r^2 (--accurate-recip)
 *   KERNEL_NO_PREFETCH:  no j-cluster prefetch (--no-prefetch), cluster-pair only
 *   KERNEL_EXCLUSIONS:   molecular exclusions (systems with molecules), cluster-pair only
 * A runtime flag inside the kernel instead costs 1.5-3%, as GCC then lays out the kernel
 * code differently. */
#ifndef __KERNEL_VARIANT_H_
#define __KERNEL_VARIANT_H_

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

#define KV_PASTE(name, a, p, x)  name##a##p##x
#define KV_EXPAND(name, a, p, x) KV_PASTE(name, a, p, x)
#define KERNEL_FN(name)          KV_EXPAND(name, KV_ACCURATE, KV_NO_PREFETCH, KV_EXCL)

// The copy of a kernel for the runtime options in param
#define KV_SELECT_EXCL(param, name) ((param)->molsize > 1 ? name##Excl : name)
#define SELECT_KERNEL(param, name)                                                       \
    ((param)->accurate_reciprocal                                                        \
            ? ((param)->kernel_prefetch                                                  \
                      ? KV_SELECT_EXCL(param, name##Accurate)                            \
                      : KV_SELECT_EXCL(param, name##AccurateNoPrefetch))                 \
            : ((param)->kernel_prefetch ? KV_SELECT_EXCL(param, name)                    \
                                        : KV_SELECT_EXCL(param, name##NoPrefetch)))

#define KV_DECLARE(name) extern double name(Parameter*, Atom*, Neighbor*, Stats*)
#define KV_DECLARE_X(name)                                                               \
    KV_DECLARE(name);                                                                    \
    KV_DECLARE(name##Excl)
#define DECLARE_KERNEL_VARIANTS(name)                                                    \
    KV_DECLARE_X(name);                                                                  \
    KV_DECLARE_X(name##Accurate);                                                        \
    KV_DECLARE_X(name##NoPrefetch);                                                      \
    KV_DECLARE_X(name##AccurateNoPrefetch)

#endif
