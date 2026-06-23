/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <parameter.h>

#ifndef __LJTABLE_H_
#define __LJTABLE_H_

// Cubic-spline tabulated Lennard-Jones forces (GROMACS-style).
//
// Instead of evaluating the LJ force-over-distance analytically, two
// *universal* force-shape functions are tabulated once and combined per
// type-pair through the prefactors (epsilon*sigma^12) and (epsilon*sigma^6):
//
//   fpair(r) = (eps*sigma^12) * Hrep(x) + (eps*sigma^6) * Gdisp(x)
//
// with the (force-over-distance) shape functions
//   Hrep  =  48 / r^14   (repulsion)
//   Gdisp = -24 / r^8    (dispersion)
//
// fx += delx * fpair, exactly as in the analytic kernels. Because the shape
// functions do not depend on the atom types, a single table serves every
// type-pair and every LJ_COMB_RULE mode (the prefactors carry the per-pair
// epsilon/sigma6 the kernels already compute).
//
// The table is indexed by x = r (default, matches GROMACS) or, when built
// with -DLJ_TABLE_RSQ, by x = r^2 (avoids the sqrt). The grid is uniform in
// x with spacing h; within interval m, eps = (x - x_m)*inv_h in [0,1) and
//   shape = c0 + eps*(c1 + eps*(c2 + eps*c3))
// (cubic Hermite, A0..A3 stored per shape). Each knot stores the 8 coeffs of
// both shapes interleaved: [rep0,rep1,rep2,rep3, disp0,disp1,disp2,disp3].

#define LJ_TABLE_STRIDE 8 // cubic coeffs per knot (4 repulsion + 4 dispersion)

#ifdef LJ_TABLE_RSQ
#define LJ_TABLE_COORD(rsq) (rsq)
#define LJ_TABLE_INDEX_NAME "rsq"
#else
#define LJ_TABLE_COORD(rsq) (sqrt(rsq))
#define LJ_TABLE_INDEX_NAME "r"
#endif

typedef struct {
    int n;             // number of intervals (knots = n + 1)
    MD_FLOAT cut;      // tabulated up to cutforce
    MD_FLOAT h;        // grid spacing in the index variable (r or r^2)
    MD_FLOAT inv_h;    // 1 / h
    MD_FLOAT* coeff;   // host coefficients, (n + 1) * LJ_TABLE_STRIDE
    MD_FLOAT* d_coeff; // device copy (GPU builds), NULL until uploaded
} LJTable;

extern LJTable ljtable;

void initLJTable(Parameter* param);

#ifdef CUDA_TARGET
void freeLJTableGPU(void);
#endif
#endif // __LJTABLE_H_
