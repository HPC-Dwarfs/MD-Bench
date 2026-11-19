/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// State and helpers of neighbor.c shared with the list-build template neighbor_build.c
#ifndef __NEIGHBOR_INTERNAL_H_
#define __NEIGHBOR_INTERNAL_H_

#include <atom.h>
#include <kernel_variant.h>
#include <neighbor.h>
#include <parameter.h>
#ifdef _OPENMP
#include <omp.h>
#endif

extern int nmax;
extern MD_FLOAT binsizex, binsizey;
extern MD_FLOAT cutneigh;
extern MD_FLOAT cutneighsq;
extern int nstencil;
extern int* stencil;
extern int* bin_clusters;
extern int* bin_nclusters;
extern int clusters_per_bin;
// Reusable scratch for pruneNeighbor*, one maxneighs-sized slice per thread
extern int* is_inner_buf;
extern int me;
extern int method;

void neighborGhost(Atom*, Neighbor*);
DECLARE_RULE_VARIANT_FUNCTIONS(void, buildNeighborCPU, Atom*, Neighbor*)

static inline int max_threads(void)
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

#endif
