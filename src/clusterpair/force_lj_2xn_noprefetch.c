/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// The 2xN kernels for --no-prefetch, see kernel_variant.h
#define KERNEL_NO_PREFETCH
#include "force_lj_2xn.c"
