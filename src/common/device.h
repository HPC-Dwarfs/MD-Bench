/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stddef.h>
//---
#include <atom.h>
#include <neighbor.h>

#ifndef __DEVICE_H_
#define __DEVICE_H_

#ifdef CUDA_TARGET
#if CUDA_TARGET == 0
#include <cuda_runtime.h>
#define error_t                cudaError_t
#define GPU_SUCCESS            cudaSuccess
#define GPU_ERROR_STR          cudaGetErrorString
#define GPU_MALLOC(p, s)       cudaMalloc(p, s)
#define GPU_MALLOC_HOST(p, s)  cudaMallocHost(p, s)
#define GPU_FREE               cudaFree
#define GPU_FREE_HOST          cudaFreeHost
#define GPU_MEMCPY(d, s, n, k) cudaMemcpy(d, s, n, k)
#define GPU_MEMSET(d, v, n)    cudaMemset(d, v, n)
#define GPU_H2D                cudaMemcpyHostToDevice
#define GPU_D2H                cudaMemcpyDeviceToHost
#define GPU_D2D                cudaMemcpyDeviceToDevice
#elif CUDA_TARGET == 1
#define __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>
#define error_t                hipError_t
#define GPU_SUCCESS            hipSuccess
#define GPU_ERROR_STR          hipGetErrorString
#define GPU_MALLOC(p, s)       hipMalloc(p, s)
#define GPU_MALLOC_HOST(p, s)  hipHostMalloc(p, s, 0)
#define GPU_FREE               hipFree
#define GPU_FREE_HOST          hipHostFree
#define GPU_MEMCPY(d, s, n, k) hipMemcpy(d, s, n, k)
#define GPU_MEMSET(d, v, n)    hipMemset(d, v, n)
#define GPU_H2D                hipMemcpyHostToDevice
#define GPU_D2H                hipMemcpyDeviceToHost
#define GPU_D2D                hipMemcpyDeviceToDevice
#endif
#ifdef __cplusplus
extern "C" {
#endif
extern void cuda_assert(const char* msg, error_t err);
#endif

extern void GPUfree(void*);
extern void initDevice(Parameter*, Atom*, Neighbor*);
extern void* allocateGPU(size_t bytesize);
extern void* reallocateGPU(void* ptr, size_t new_bytesize);
extern void memcpyToGPU(void* d_ptr, void* h_ptr, size_t bytesize);
extern void memcpyFromGPU(void* h_ptr, void* d_ptr, size_t bytesize);
extern void memcpyOnGPU(void* d_dst, void* d_src, size_t bytesize);
extern void memsetGPU(void* d_ptr, int value, size_t bytesize);
#ifdef __cplusplus
}
#endif
#endif
