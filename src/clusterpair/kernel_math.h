/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Reciprocal and reciprocal square root helpers of the cluster-pair SIMD kernels, see
 * simd_math.h (shared with the verlet-list SIMD kernels).
 */
#ifndef __KERNEL_MATH_H_
#define __KERNEL_MATH_H_

#include <parameter.h>

#ifndef USE_REFERENCE_KERNEL
// force.h defines the cluster dimensions that simd.h would otherwise default
#include <force.h>
#include <simd.h>
#include <simd_math.h>
#endif /* !USE_REFERENCE_KERNEL */

#endif
