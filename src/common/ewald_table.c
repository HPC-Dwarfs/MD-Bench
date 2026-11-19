/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// GROMACS' generateEwaldCorrectionTables() (tables/forcetable.cpp, LGPL)
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <allocate.h>
#include <ewald_table.h>

#if PRECISION == 1
#define TABLE_REAL_EPS FLT_EPSILON
#define TABLE_REAL_MAX FLT_MAX
#else
#define TABLE_REAL_EPS DBL_EPSILON
#define TABLE_REAL_MAX DBL_MAX
#endif

// erf(beta r) / r, its limit at r = 0
static double ewaldLongRange(double beta, double r)
{
    return r == 0 ? beta * 2.0 / sqrt(M_PI) : erf(beta * r) / r;
}

// Table points per length unit (GROMACS' spline3_table_scale())
static double splineTableScale(double third_deriv_max, double x_scale, double func_tol)
{
    const double deriv_tol = FLT_EPSILON;
    const double sc_deriv  = sqrt(third_deriv_max / (6 * 4 * deriv_tol * x_scale)) *
                            x_scale;
    func_tol             = func_tol > TABLE_REAL_EPS ? func_tol : TABLE_REAL_EPS;
    const double sc_func = cbrt(third_deriv_max / (6 * 12 * sqrt(3.0) * func_tol)) *
                           x_scale;
    return sc_deriv > sc_func ? sc_deriv : sc_func;
}

void initEwaldTable(Parameter* param)
{
    const MD_FLOAT beta   = param->ewald_beta;
    const MD_FLOAT cutoff = param->cutforce;
    // Energy tolerance 0.1 times the cut-off jump, max of (erf(x)/x)''' is 1.0522
    const double etol    = 0.1 * erfc((double)beta * cutoff);
    const MD_FLOAT scale = (MD_FLOAT)splineTableScale(1.0522, beta, etol);
    const int n          = (int)(cutoff * scale) + 2;
    const double dx      = 1.0 / (double)scale;
    const double tab_max = TABLE_REAL_MAX * 0.0001;

    double* tf = (double*)malloc(n * sizeof(double));
    double* tv = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        tf[i] = 0.0;
    }

    // Quadratic splines from the end of the table to r = 0, F averaged over two intervals
    int out_of_range = 0;
    int i_inrange    = n;
    double v_inrange = 0.0;
    double dc        = 0.0;
    for (int i = n - 1; i >= 0; i--) {
        const double x_r0 = i * dx;
        const double v_r0 = ewaldLongRange(beta, x_r0);
        double vi;
        if (!out_of_range) {
            i_inrange = i;
            v_inrange = v_r0;
            vi        = v_r0;
        } else {
            vi = v_inrange - dc * (i - i_inrange) * dx;
        }
        tv[i] = vi;
        if (i == 0) {
            continue;
        }

        const double v_r1 = ewaldLongRange(beta, (i - 1) * dx);
        if (v_r1 != v_r1 || v_r1 < -tab_max || v_r1 > tab_max) {
            out_of_range = 1;
        }
        if (!out_of_range) {
            const double a2dx = (v_r0 + v_r1 -
                                    2 * ewaldLongRange(beta, x_r0 - 0.5 * dx)) /
                                (0.25 * dx);
            dc = (v_r0 - v_r1) / dx + 0.5 * a2dx;
        }
        if (i == n - 1) {
            tf[i] = -dc;
        } else {
            tf[i] += -0.5 * dc;
        }
        if (!out_of_range) {
            const double a0     = v_r0;
            const double a1     = dc;
            const double a2dx   = (a1 * dx + v_r1 - a0) * 2 / dx;
            const double dc_new = a1 - a2dx;
            if (dc_new != dc_new || dc_new < -tab_max || dc_new > tab_max) {
                out_of_range = 1;
            } else {
                dc = dc_new;
            }
        }
        tf[i - 1] = -0.5 * dc;
    }
    tf[0] *= 2;

    MD_FLOAT* f    = (MD_FLOAT*)allocate(ALIGNMENT, (n + 1) * sizeof(MD_FLOAT));
    MD_FLOAT* v    = (MD_FLOAT*)allocate(ALIGNMENT, (n + 1) * sizeof(MD_FLOAT));
    MD_FLOAT* fdv0 = (MD_FLOAT*)allocate(ALIGNMENT, 4 * n * sizeof(MD_FLOAT));
    for (int i = 0; i < n; i++) {
        f[i] = (MD_FLOAT)tf[i];
        v[i] = (MD_FLOAT)tv[i];
    }
    // One extra point so that the F[i+1] of the last point can be read
    f[n] = f[n - 1];
    v[n] = v[n - 1];
    for (int i = 0; i < n - 1; i++) {
        fdv0[4 * i]     = f[i];
        fdv0[4 * i + 1] = f[i + 1] - f[i];
        fdv0[4 * i + 2] = v[i];
        fdv0[4 * i + 3] = 0.0;
    }
    fdv0[4 * (n - 1)]     = f[n - 1];
    fdv0[4 * (n - 1) + 1] = -f[n - 1];
    fdv0[4 * (n - 1) + 2] = v[n - 1];
    fdv0[4 * (n - 1) + 3] = 0.0;
    free(tf);
    free(tv);

    param->ewald_table_scale = scale;
    param->ewald_table_size  = n;
    param->ewald_table_f     = f;
    param->ewald_table_v     = v;
    param->ewald_table_fdv0  = fdv0;
}
