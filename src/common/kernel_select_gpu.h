/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// Selection of the GPU force kernel instantiation for the runtime options: the kernels
// are templates on the LJ combination rule (RULE), the energies (ENERGY), the molecular
// exclusions (EXCL) and the IEEE reciprocals (ACCURATE, --accurate-recip), so none of
// these options costs a branch inside the kernel.
#ifndef __KERNEL_SELECT_GPU_H_
#define __KERNEL_SELECT_GPU_H_

#include <parameter.h>

// K<R> for the combination rule rule (LJ_COMB_SINGLE/GEOM/FULL)
#define SELECT_GPU_KERNEL_RULE(K, rule)                                                  \
    ((rule) == LJ_COMB_SINGLE    ? K<LJ_COMB_SINGLE>                                     \
        : (rule) == LJ_COMB_GEOM ? K<LJ_COMB_GEOM>                                       \
                                 : K<LJ_COMB_FULL>)

// K<R, E, excl, accurate> for a fixed rule R and energies E
#define SELECT_GPU_KERNEL_XA(K, R, E, excl, accurate)                                    \
    ((excl) ? ((accurate) ? K<R, E, true, true> : K<R, E, true, false>)                  \
            : ((accurate) ? K<R, E, false, true> : K<R, E, false, false>))

// K<R, energy, excl, accurate> for a fixed rule R
#define SELECT_GPU_KERNEL_EXA(K, R, energy, excl, accurate)                              \
    ((energy) ? SELECT_GPU_KERNEL_XA(K, R, true, excl, accurate)                         \
              : SELECT_GPU_KERNEL_XA(K, R, false, excl, accurate))

// K<rule, energy, excl, accurate> for the combination rule, energies, exclusions and
// reciprocals
#define SELECT_GPU_KERNEL(K, rule, energy, excl, accurate)                               \
    ((rule) == LJ_COMB_SINGLE                                                            \
            ? SELECT_GPU_KERNEL_EXA(K, LJ_COMB_SINGLE, energy, excl, accurate)           \
        : (rule) == LJ_COMB_GEOM                                                         \
            ? SELECT_GPU_KERNEL_EXA(K, LJ_COMB_GEOM, energy, excl, accurate)             \
            : SELECT_GPU_KERNEL_EXA(K, LJ_COMB_FULL, energy, excl, accurate))

#endif
