/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <eam.h>
#include <kernel_variant.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>

#ifndef __FORCE_H_
#define __FORCE_H_

typedef double (*ComputeForceFunction)(Parameter*, Atom*, Neighbor*, Stats*);
extern ComputeForceFunction computeForce;

enum forcetype { FF_LJ = 0, FF_EAM, FF_LJ_TABLE };

extern void initForce(Parameter*);
// Copies per LJ combination rule (and exclusions / --accurate-recip), see
// kernel_variant.h
DECLARE_RULE_EXCL_VARIANTS(computeForceLJHalfNeigh)
DECLARE_RULE_EXCL_VARIANTS(computeForceLJFullNeigh)
DECLARE_RULE_VARIANTS(computeForceLJTableHalfNeigh)
DECLARE_RULE_VARIANTS(computeForceLJTableFullNeigh)
DECLARE_RULE_VARIANTS(computeForceEam)

// The intrinsics kernels (--force-simd) need a SIMD header (src/common/simd/), none for
// SSE and scalar builds; SVE/SVE2 use the single-loop VLA kernel instead of main+tail
#if !defined(CUDA_TARGET) &&                                                             \
    (defined(__ISA_AVX__) || defined(__ISA_AVX2__) || defined(__ISA_AVX512__) ||         \
        defined(__ISA_NEON__) || defined(__ISA_SVE__))
#define HAVE_SIMD_FORCE
#if defined(__ISA_SVE__)
#define SIMD_FORCE_VLA
#endif
#endif

#ifdef HAVE_SIMD_FORCE
DECLARE_RULE_ACCURATE_VARIANTS(computeForceLJHalfNeigh_simd_compress)
DECLARE_RULE_ACCURATE_VARIANTS(computeForceLJFullNeigh_simd_compress)
#ifdef SIMD_FORCE_VLA
DECLARE_RULE_ACCURATE_VARIANTS(computeForceLJHalfNeigh_simd_vla)
DECLARE_RULE_ACCURATE_VARIANTS(computeForceLJFullNeigh_simd_vla)
#else
DECLARE_RULE_ACCURATE_VARIANTS(computeForceLJHalfNeigh_simd)
DECLARE_RULE_ACCURATE_VARIANTS(computeForceLJFullNeigh_simd)
#endif
#endif

#ifdef CUDA_TARGET
extern double computeForceLJCUDA(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJCUDARef(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJTableCUDA(Parameter*, Atom*, Neighbor*, Stats*);
#define KERNEL_NAME "CUDA"
#endif

#endif // __FORCE_H_
