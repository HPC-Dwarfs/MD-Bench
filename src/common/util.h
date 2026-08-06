/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#ifndef __UTIL_H_
#define __UTIL_H_

#include <stdio.h>
#include <string.h>
#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef ABS
#define ABS(a) ((a) >= 0 ? (a) : -(a))
#endif

#define DEBUG_MESSAGE debug_printf

#ifndef MAXLINE
#define MAXLINE 4096
#endif

#if PRECISION == 1
#define PRECISION_STRING "single"
#else
#define PRECISION_STRING "double"
#endif

// Whole string matching (exact match)
#define STR_EQ(s, lit)                                                                   \
    (strncmp((s), "" lit "", sizeof(lit) - 1) == 0 && (s)[sizeof(lit) - 1] == '\0')

enum { fullShell = 0, halfShell, eightShell, halfStencil };

#define BigOrEqual(a, b) (fabs((a) - (b)) < 1e-9 || (a) > (b))
#define Equal(a, b)      (fabs((a) - (b)) < 1e-6)

extern double myrandom(int*);
extern void random_reset(int* seed, int ibase, double* coord);
extern int str2ff(const char* string);
extern const char* ff2str(int ff);
extern void readline(char* line, FILE* fp);
extern void debug_printf(const char* format, ...);
extern int get_cuda_num_threads(void);
extern void fprintf_once(int me, FILE* stream, const char* format, ...);
#endif
