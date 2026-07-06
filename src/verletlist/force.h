/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <eam.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>

#ifndef __FORCE_H_
#define __FORCE_H_

typedef double (*ComputeForceFunction)(Parameter*, Atom*, Neighbor*, Stats*);
extern ComputeForceFunction computeForce;

enum forcetype { FF_LJ = 0, FF_EAM, FF_LJ_TABLE };

extern void initForce(Parameter*);
extern double computeForceLJHalfNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJFullNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJTableHalfNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJTableFullNeigh(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceEam(Parameter*, Atom*, Neighbor*, Stats*);

#ifdef __SIMD_KERNEL__
extern double computeForceLJHalfNeigh_simd(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJFullNeigh_simd(Parameter*, Atom*, Neighbor*, Stats*);
#ifdef __SIMD_COMPRESS__
extern double computeForceLJHalfNeigh_simd_compress(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJFullNeigh_simd_compress(Parameter*, Atom*, Neighbor*, Stats*);
#endif
#endif

#ifdef CUDA_TARGET
extern double computeForceLJCUDA(Parameter*, Atom*, Neighbor*, Stats*);
extern double computeForceLJTableCUDA(Parameter*, Atom*, Neighbor*, Stats*);
#define KERNEL_NAME "CUDA"
#else
#ifdef __SIMD_KERNEL__
#define KERNEL_NAME "SIMD"
#else
#define KERNEL_NAME "SCALAR"
#endif
#endif

#endif // __FORCE_H_
