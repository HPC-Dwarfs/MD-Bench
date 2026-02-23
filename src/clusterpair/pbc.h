/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdbool.h>

#include <atom.h>
#include <parameter.h>

#ifndef __PBC_H_
#define __PBC_H_

typedef void (*UpdatePbcFunction)(Atom*, Parameter*, bool);
extern UpdatePbcFunction updatePbc;
extern UpdatePbcFunction updateAtomsPbc;

extern void initPbc(Atom*);
extern void updatePbcCPU(Atom*, Parameter*, bool);
extern void updateAtomsPbcCPU(Atom*, Parameter*, bool);
extern void setupPbc(Atom*, Parameter*);

#ifdef CUDA_TARGET
extern void updatePbcCUDA(Atom*, Parameter*, bool);
// extern void updateAtomsPbcCUDA(Atom*, Parameter*, bool);
extern void setupPbcGPU(Atom*, Parameter*);
#endif
#endif
