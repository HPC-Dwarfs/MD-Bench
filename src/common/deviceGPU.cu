/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//---
#include <device.h>

void cuda_assert(const char* label, error_t err)
{
    if (err != GPU_SUCCESS) {
        printf("[GPU Error]: %s: %s\r\n", label, GPU_ERROR_STR(err));
        exit(-1);
    }
}

void GPUfree(void* any) { cuda_assert("GPUfree", GPU_FREE(any)); }

void* allocateGPU(size_t bytesize)
{
    void* ptr;
#ifdef CUDA_HOST_MEMORY
    cuda_assert("allocateGPU", GPU_MALLOC_HOST((void**)&ptr, bytesize));
#else
    cuda_assert("allocateGPU", GPU_MALLOC((void**)&ptr, bytesize));
#endif
    return ptr;
}

void* reallocateGPU(void* ptr, size_t new_bytesize)
{
    if (ptr != NULL) {
#ifdef CUDA_HOST_MEMORY
        (void)GPU_FREE_HOST(ptr);
#else
        (void)GPU_FREE(ptr);
#endif
    }
    return allocateGPU(new_bytesize);
}

// Growing realloc that preserves the old contents. reallocateGPU() above is
// free+malloc and loses the data, which is only safe for buffers that are
// refilled right after; growAtom() must not lose device state (e.g. velocities
// are never re-uploaded in non-MPI runs).
void* reallocateGPUKeep(void* ptr, size_t new_bytesize, size_t old_bytesize)
{
    void* new_ptr = allocateGPU(new_bytesize);
    if (ptr != NULL) {
        if (old_bytesize > 0) {
#ifdef CUDA_HOST_MEMORY
            memcpy(new_ptr, ptr, old_bytesize);
#else
            cuda_assert("reallocateGPUKeep",
                GPU_MEMCPY(new_ptr, ptr, old_bytesize, GPU_D2D));
#endif
        }
#ifdef CUDA_HOST_MEMORY
        (void)GPU_FREE_HOST(ptr);
#else
        (void)GPU_FREE(ptr);
#endif
    }
    return new_ptr;
}

void memcpyToGPU(void* d_ptr, void* h_ptr, size_t bytesize)
{
#ifndef CUDA_HOST_MEMORY
    cuda_assert("memcpyToGPU", GPU_MEMCPY(d_ptr, h_ptr, bytesize, GPU_H2D));
#endif
}

void memcpyFromGPU(void* h_ptr, void* d_ptr, size_t bytesize)
{
#ifndef CUDA_HOST_MEMORY
    cuda_assert("memcpyFromGPU", GPU_MEMCPY(h_ptr, d_ptr, bytesize, GPU_D2H));
#endif
}

void memcpyOnGPU(void* d_dst, void* d_src, size_t bytesize)
{
    cuda_assert("memcpyOnGPU", GPU_MEMCPY(d_dst, d_src, bytesize, GPU_D2D));
}

void memsetGPU(void* d_ptr, int value, size_t bytesize)
{
    cuda_assert("memsetGPU", GPU_MEMSET(d_ptr, value, bytesize));
}
