/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <force.h>
#include <util.h>

/* Park/Miller RNG w/out MASKING, so as to be like f90s version */
#define IA   16807
#define IM   2147483647
#define AM   (1.0 / IM)
#define IQ   127773
#define IR   2836
#define MASK 123459876

double myrandom(int* seed) {
    int k = (*seed) / IQ;
    double ans;

    *seed = IA * (*seed - k * IQ) - IR * k;
    if (*seed < 0) *seed += IM;
    ans = AM * (*seed);
    return ans;
}

void random_reset(int* seed, int ibase, double* coord) {
    int i;
    char* str         = (char*)&ibase;
    int n             = sizeof(int);
    unsigned int hash = 0;

    for (i = 0; i < n; i++) {
        hash += str[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }

    str = (char*)coord;
    n   = 3 * sizeof(double);
    for (i = 0; i < n; i++) {
        hash += str[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }

    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

    // keep 31 bits of unsigned int as new seed
    // do not allow seed = 0, since will cause hang in gaussian()

    *seed = hash & 0x7ffffff;
    if (!(*seed)) *seed = 1;

    // warm up the RNG

    for (i = 0; i < 5; i++)
        myrandom(seed);
    // save = 0;
}

int str2ff(const char* string) {
    if (strncmp(string, "lj", 2) == 0) return FF_LJ;
    if (strncmp(string, "eam", 3) == 0) return FF_EAM;
    return -1;
}

const char* ff2str(int ff) {
    if (ff == FF_LJ) {
        return "lj";
    }

    if (ff == FF_EAM) {
        return "eam";
    }

    return "invalid";
}

int get_cuda_num_threads(void) {
    const char* num_threads_env = getenv("NUM_THREADS");
    return (num_threads_env == NULL) ? 128 : atoi(num_threads_env);
}

void readline(char* line, FILE* fp) {
    if (fgets(line, MAXLINE, fp) == NULL) {
        if (errno != 0) {
            perror("readline()");
            exit(-1);
        }
    }
}

void debug_printf(const char* format, ...) {
#ifdef DEBUG
    va_list arg;
    int ret;

    va_start(arg, format);
    if ((vfprintf(stdout, format, arg)) < 0) {
        perror("debug_printf()");
    }
    va_end(arg);
#endif
}

void fprintf_once(int me, FILE* stream, const char* format, ...)
{
    if(me == 0) {
        va_list arg;
        int ret;

        va_start(arg, format);
        if ((vfprintf(stream, format, arg)) < 0) {
            perror("debug_printf()");
        }
        va_end(arg);
    }
}
