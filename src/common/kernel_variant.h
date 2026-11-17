/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Compiled copies of the CPU SIMD force kernels, selected by runtime options without a
 * branch in the kernel: wrapper files (force_lj_*_accurate.c, *_noprefetch.c, ...)
 * define the macros below and include the kernel file, KERNEL_FN gives the entry points
 * of each copy its own name and initForce() picks one.
 *   ACCURATE_RECIPROCAL: Newton-Raphson refined 1/r and 1/r^2 (--accurate-recip)
 *   KERNEL_NO_PREFETCH:  no j-cluster prefetch, cluster-pair only (--no-prefetch)
 */
#ifndef __KERNEL_VARIANT_H_
#define __KERNEL_VARIANT_H_

#if defined(ACCURATE_RECIPROCAL) && defined(KERNEL_NO_PREFETCH)
#define KERNEL_FN(name) name##AccurateNoPrefetch
#elif defined(ACCURATE_RECIPROCAL)
#define KERNEL_FN(name) name##Accurate
#elif defined(KERNEL_NO_PREFETCH)
#define KERNEL_FN(name) name##NoPrefetch
#else
#define KERNEL_FN(name) name
#endif

// The copy of a kernel for the runtime options in param
#define SELECT_KERNEL(param, name)                                                       \
    ((param)->accurate_reciprocal                                                        \
            ? ((param)->kernel_prefetch ? name##Accurate : name##AccurateNoPrefetch)     \
            : ((param)->kernel_prefetch ? name : name##NoPrefetch))

#define DECLARE_KERNEL_VARIANTS(name)                                                    \
    extern double name(Parameter*, Atom*, Neighbor*, Stats*);                            \
    extern double name##Accurate(Parameter*, Atom*, Neighbor*, Stats*);                  \
    extern double name##NoPrefetch(Parameter*, Atom*, Neighbor*, Stats*);                \
    extern double name##AccurateNoPrefetch(Parameter*, Atom*, Neighbor*, Stats*)

#endif
