/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <neighbor.h>
#include <parameter.h>

#if defined(MEM_TRACER) || defined(INDEX_TRACER)
#include <stdio.h>
#include <stdlib.h>
#endif

#ifndef VECTOR_WIDTH
#define VECTOR_WIDTH 8
#endif

#ifndef TRACER_CONDITION
#define TRACER_CONDITION (!(timestep % param->every))
#endif

#ifdef MEM_TRACER
#define MEM_TRACER_INIT                                                                  \
    FILE* mem_tracer_fp;                                                                 \
    if (TRACER_CONDITION) {                                                              \
        char mem_tracer_fn[128];                                                         \
        snprintf(mem_tracer_fn, sizeof mem_tracer_fn, "mem_tracer_%d.out", timestep);    \
        mem_tracer_fp = fopen(mem_tracer_fn, "w");
}

#define MEM_TRACER_END                                                                   \
    if (TRACER_CONDITION) {                                                              \
        fclose(mem_tracer_fp);                                                           \
    }
#define MEM_TRACE(addr, op)                                                              \
    if (TRACER_CONDITION) {                                                              \
        fprintf(mem_tracer_fp, "%c: %p\n", op, (void*)(&(addr)));                        \
    }
#else
#define MEM_TRACER_INIT
#define MEM_TRACER_END
#define MEM_TRACE(addr, op)
#endif

#ifdef INDEX_TRACER
#define INDEX_TRACER_INIT                                                                \
    FILE* index_tracer_fp;                                                               \
    if (TRACER_CONDITION) {                                                              \
        char index_tracer_fn[128];                                                       \
        snprintf(index_tracer_fn,                                                        \
            sizeof index_tracer_fn,                                                      \
            "index_tracer_%d.out",                                                       \
            timestep);                                                                   \
        index_tracer_fp = fopen(index_tracer_fn, "w");                                   \
    }

#define INDEX_TRACER_END                                                                 \
    if (TRACER_CONDITION) {                                                              \
        fclose(index_tracer_fp);                                                         \
    }
#define INDEX_TRACE_NATOMS(nl, ng, mn)                                                   \
    if (TRACER_CONDITION) {                                                              \
        fprintf(index_tracer_fp, "N: %d %d %d\n", nl, ng, mn);                           \
    }
#define INDEX_TRACE_ATOM(a)                                                              \
    if (TRACER_CONDITION) {                                                              \
        fprintf(index_tracer_fp, "A: %d\n", a);                                          \
    }
#define INDEX_TRACE(l, i, e, M, N)                                                       \
    if (TRACER_CONDITION) {                                                              \
        for (int __i = 0; __i < (e); __i += VECTOR_WIDTH) {                              \
            int __e = (((e)-__i) < VECTOR_WIDTH) ? ((e)-__i) : VECTOR_WIDTH;             \
            fprintf(index_tracer_fp, "I: ");                                             \
            for (int __j = 0; __j < __e; ++__j) {                                        \
                fprintf(index_tracer_fp, "%d ", neighs(l, i, __i + __j, M, N));          \
            }                                                                            \
            fprintf(index_tracer_fp, "\n");                                              \
        }                                                                                \
    }

#define DIST_TRACE_SORT(l, i, e, M, N)                                                   \
    if (TRACER_CONDITION) {                                                              \
        for (int __i = 0; __i < (e); __i += VECTOR_WIDTH) {                              \
            int __e = (((e)-__i) < VECTOR_WIDTH) ? ((e)-__i) : VECTOR_WIDTH;             \
            if (__e > 1) {                                                               \
                for (int __j = __i; __j < __i + __e - 1; ++__j) {                        \
                    for (int __k = __i; __k < __i + __e - (__j - __i) - 1; ++__k) {      \
                        if (neighs(l, i, __k, M, N) > neighs(l, i, __k + 1, M, N)) {     \
                            int __t    = neighs(l, i, __k, M, N);                        \
                            neighs(l, i, __k, M, N)     = neighs(l, i, __k + 1, M, N);   \
                            neighs(l, i, __k + 1, M, N) = __t;                           \
                        }                                                                \
                    }                                                                    \
                }                                                                        \
            }                                                                            \
        }                                                                                \
    }

#define DIST_TRACE(l, i, e, M, N)                                                        \
    if (TRACER_CONDITION) {                                                              \
        for (int __i = 0; __i < (e); __i += VECTOR_WIDTH) {                              \
            int __e = (((e)-__i) < VECTOR_WIDTH) ? ((e)-__i) : VECTOR_WIDTH;             \
            if (__e > 1) {                                                               \
                fprintf(index_tracer_fp, "D: ");                                         \
                for (int __j = 0; __j < __e - 1; ++__j) {                                \
                    int __dist = abs(neighs(l, i, __i + __j + 1, M, N) -                 \
                                     neighs(l, i, __i + __j, M, N));                     \
                    fprintf(index_tracer_fp, "%d ", __dist);                             \
                }                                                                        \
                fprintf(index_tracer_fp, "\n");                                          \
            }                                                                            \
        }                                                                                \
    }
#else
#define INDEX_TRACER_INIT
#define INDEX_TRACER_END
#define INDEX_TRACE_NATOMS(nl, ng, mn)
#define INDEX_TRACE_ATOM(a)
#define INDEX_TRACE(l, i, e, M, N)
#define DIST_TRACE_SORT(l, i, e, M, N)
#define DIST_TRACE(l, i, e, M, N)
#endif

extern void traceAddresses(
    Parameter* param, Atom* atom, Neighbor* neighbor, int timestep);
