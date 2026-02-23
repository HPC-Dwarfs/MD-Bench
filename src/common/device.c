/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
//---
#include <device.h>

#ifndef CUDA_TARGET
void GPUfree(void* any) {}
void initDevice(Parameter* param, Atom* atom, Neighbor* neighbor) {}
void* allocateGPU(size_t bytesize) { return NULL; }
void* reallocateGPU(void* ptr, size_t new_bytesize) { return NULL; }
void memcpyToGPU(void* d_ptr, void* h_ptr, size_t bytesize) {}
void memcpyFromGPU(void* h_ptr, void* d_ptr, size_t bytesize) {}
void memsetGPU(void* d_ptr, int value, size_t bytesize) {}
#endif
