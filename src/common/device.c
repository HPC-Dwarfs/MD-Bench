/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
//---
#include <allocate.h>
#include <device.h>

#ifndef CUDA_TARGET
void GPUfree(void* any) {}
void initDeviceContext(void) {}
void initDevice(Parameter* param, Atom* atom, Neighbor* neighbor) {}
void* allocateGPU(size_t bytesize) { return NULL; }
void* reallocateGPU(void* ptr, size_t new_bytesize) { return NULL; }
void* reallocateGPUKeep(void* ptr, size_t new_bytesize, size_t old_bytesize)
{
    return NULL;
}
void memcpyToGPU(void* d_ptr, void* h_ptr, size_t bytesize) {}
void memcpyFromGPU(void* h_ptr, void* d_ptr, size_t bytesize) {}
void memcpyOnGPU(void* d_dst, void* d_src, size_t bytesize) {}
void memsetGPU(void* d_ptr, int value, size_t bytesize) {}

// No GPU to transfer to/from on this build, so these are plain host
// allocations; callers (e.g. growClusters()) use them unconditionally
// regardless of target.
void* allocateHostPinned(size_t bytesize) { return allocate(ALIGNMENT, bytesize); }
void* reallocateHostPinned(void* ptr, size_t new_bytesize, size_t old_bytesize)
{
    return reallocate(ptr, ALIGNMENT, new_bytesize, old_bytesize);
}
void freeHostPinned(void* ptr) { free(ptr); }
#endif
