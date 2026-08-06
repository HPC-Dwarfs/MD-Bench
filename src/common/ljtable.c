/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <allocate.h>
#include <ljtable.h>
#include <parameter.h>

LJTable ljtable;

// Evaluate the two force-shape functions and their first derivatives at the
// index coordinate x. For r-indexed tables x is the distance r, for rsq tables
// x is r^2. Returns force-over-distance shapes so that
//   fpair = (eps*sigma^12)*rep + (eps*sigma^6)*disp.
static void ljShape(
    MD_FLOAT x, MD_FLOAT* rep, MD_FLOAT* drep, MD_FLOAT* disp, MD_FLOAT* ddisp)
{
#ifdef LJ_TABLE_RSQ
    // x = r^2 :  Hrep = 48 x^-7,  Gdisp = -24 x^-4
    MD_FLOAT xi  = 1.0 / x;
    MD_FLOAT xi4 = xi * xi * xi * xi;
    MD_FLOAT xi7 = xi4 * xi * xi * xi;
    *rep         = 48.0 * xi7;
    *drep        = -7.0 * 48.0 * xi7 * xi; // d/dx (48 x^-7) = -336 x^-8
    *disp        = -24.0 * xi4;
    *ddisp       = 4.0 * 24.0 * xi4 * xi; // d/dx (-24 x^-4) = 96 x^-5
#else
    // x = r :  Hrep = 48 r^-14,  Gdisp = -24 r^-8
    MD_FLOAT ri   = 1.0 / x;
    MD_FLOAT ri2  = ri * ri;
    MD_FLOAT ri8  = ri2 * ri2 * ri2 * ri2;
    MD_FLOAT ri14 = ri8 * ri2 * ri2 * ri2;
    *rep          = 48.0 * ri14;
    *drep         = -14.0 * 48.0 * ri14 * ri; // d/dr (48 r^-14) = -672 r^-15
    *disp         = -24.0 * ri8;
    *ddisp        = 8.0 * 24.0 * ri8 * ri; // d/dr (-24 r^-8) = 192 r^-9
#endif
}

void initLJTable(Parameter* param)
{
    int n = param->lj_table_points;
    if (n < 8) {
        n = 1000;
    }

    MD_FLOAT cut = param->cutforce;
#ifdef LJ_TABLE_RSQ
    MD_FLOAT xmax = cut * cut;
#else
    MD_FLOAT xmax = cut;
#endif
    MD_FLOAT h = xmax / (MD_FLOAT)n;

    ljtable.n       = n;
    ljtable.cut     = cut;
    ljtable.h       = h;
    ljtable.inv_h   = 1.0 / h;
    ljtable.d_coeff = NULL;

    int nknots    = n + 1;
    ljtable.coeff = (MD_FLOAT*)allocate(ALIGNMENT,
        nknots * LJ_TABLE_STRIDE * sizeof(MD_FLOAT));
    MD_FLOAT* c   = ljtable.coeff;

    // Floor the index coordinate away from the r->0 singularity. Atoms never
    // approach this closely, so these inner knots are never reached in practice
    // (the force kernels still gate on rsq < cutforcesq); the floor only keeps
    // the stored coefficients finite.
    MD_FLOAT xfloor = 1.0e-2 * h;
    if (xfloor < 1.0e-6) {
        xfloor = 1.0e-6;
    }

    // Per-knot value and eps-scaled derivative for both shapes.
    MD_FLOAT* rv = (MD_FLOAT*)allocate(ALIGNMENT, nknots * sizeof(MD_FLOAT));
    MD_FLOAT* rd = (MD_FLOAT*)allocate(ALIGNMENT, nknots * sizeof(MD_FLOAT));
    MD_FLOAT* dv = (MD_FLOAT*)allocate(ALIGNMENT, nknots * sizeof(MD_FLOAT));
    MD_FLOAT* dd = (MD_FLOAT*)allocate(ALIGNMENT, nknots * sizeof(MD_FLOAT));

    for (int i = 0; i < nknots; i++) {
        MD_FLOAT x = i * h;
        if (x < xfloor) {
            x = xfloor;
        }
        MD_FLOAT rep, drep, disp, ddisp;
        ljShape(x, &rep, &drep, &disp, &ddisp);
        rv[i] = rep;
        dv[i] = disp;
        // derivative wrt the local coordinate eps = (x - x_m)/h is h * d/dx
        rd[i] = h * drep;
        dd[i] = h * ddisp;
    }

    // Build cubic Hermite coefficients per interval m in [0, n-1]:
    //   A0 = f0, A1 = d0, A2 = 3(f1-f0)-2d0-d1, A3 = 2(f0-f1)+d0+d1
    for (int m = 0; m < n; m++) {
        MD_FLOAT* ci = &c[m * LJ_TABLE_STRIDE];
        MD_FLOAT rf0 = rv[m], rf1 = rv[m + 1], rd0 = rd[m], rd1 = rd[m + 1];
        ci[0]        = rf0;
        ci[1]        = rd0;
        ci[2]        = 3.0 * (rf1 - rf0) - 2.0 * rd0 - rd1;
        ci[3]        = 2.0 * (rf0 - rf1) + rd0 + rd1;
        MD_FLOAT df0 = dv[m], df1 = dv[m + 1], dd0 = dd[m], dd1 = dd[m + 1];
        ci[4] = df0;
        ci[5] = dd0;
        ci[6] = 3.0 * (df1 - df0) - 2.0 * dd0 - dd1;
        ci[7] = 2.0 * (df0 - df1) + dd0 + dd1;
    }

    // Defensive copy for the final (unused) knot slot; the kernels clamp the
    // interval index to n-1 so this is never read, but keep it well-defined.
    for (int k = 0; k < LJ_TABLE_STRIDE; k++) {
        c[n * LJ_TABLE_STRIDE + k] = c[(n - 1) * LJ_TABLE_STRIDE + k];
    }

    free(rv);
    free(rd);
    free(dv);
    free(dd);

    printf("Initialized tabulated LJ forces: %d intervals over [0, %.4f] "
           "(index=%s, spacing=%.3e)\n",
        n,
        (double)cut,
        LJ_TABLE_INDEX_NAME,
        (double)h);
}
