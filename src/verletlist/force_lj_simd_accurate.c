/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// The SIMD kernels for --accurate-recip, see kernel_variant.h
#ifdef __SIMD_KERNEL__
#define ACCURATE_RECIPROCAL
#include "force_lj_simd.c"
#endif
