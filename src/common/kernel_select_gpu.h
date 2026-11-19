/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// Selection of the GPU force kernel instantiation for the runtime options: the kernels
// are templates on the LJ combination rule (RULE), the energies (ENERGY) and the
// molecular exclusions (EXCL), so none of these options costs a branch inside the kernel.
#ifndef __KERNEL_SELECT_GPU_H_
#define __KERNEL_SELECT_GPU_H_

#include <parameter.h>

// K<R> for the combination rule rule (LJ_COMB_SINGLE/GEOM/FULL)
#define SELECT_GPU_KERNEL_RULE(K, rule)                                                  \
    ((rule) == LJ_COMB_SINGLE    ? K<LJ_COMB_SINGLE>                                     \
        : (rule) == LJ_COMB_GEOM ? K<LJ_COMB_GEOM>                                       \
                                 : K<LJ_COMB_FULL>)

// K<R, energy, excl> for a fixed rule R
#define SELECT_GPU_KERNEL_EE(K, R, energy, excl)                                         \
    ((energy) ? ((excl) ? K<R, true, true> : K<R, true, false>)                          \
              : ((excl) ? K<R, false, true> : K<R, false, false>))

// K<rule, energy, excl> for the combination rule, energies and exclusions
#define SELECT_GPU_KERNEL(K, rule, energy, excl)                                         \
    ((rule) == LJ_COMB_SINGLE    ? SELECT_GPU_KERNEL_EE(K, LJ_COMB_SINGLE, energy, excl) \
        : (rule) == LJ_COMB_GEOM ? SELECT_GPU_KERNEL_EE(K, LJ_COMB_GEOM, energy, excl)   \
                                 : SELECT_GPU_KERNEL_EE(K, LJ_COMB_FULL, energy, excl))

#endif
