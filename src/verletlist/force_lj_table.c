/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>

#include <atom.h>
#include <likwid-marker.h>
#include <ljtable.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>

// Tabulated/spline-interpolated Lennard-Jones forces (GROMACS-style).
// The analytic sr2/sr6/force arithmetic of computeForceLJ*Neigh is replaced by
// a cubic-spline lookup of two universal force-shape functions combined with
// the per-pair prefactors (epsilon*sigma^12) and (epsilon*sigma^6). See
// src/common/ljtable.{h,c} for the table layout.

double computeForceLJTableFullNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJTableFullNeigh begin\n");

    int nlocal = atom->Nlocal;
    int* neighs;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif
    const int nm1         = ljtable.n - 1;
    const MD_FLOAT inv_h  = ljtable.inv_h;
    const MD_FLOAT* coeff = ljtable.coeff;

    for (int i = 0; i < nlocal; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }
    double timeStart = getTimeStamp();

#pragma omp parallel
    {
        // restrict pointers so GCC can vectorize this loop at full width even
        // when outlined for OpenMP (see atom_access_x/y/z in atom.h)
#ifdef ATOM_POSITION_AOS
        MD_FLOAT* restrict x_ptr  = atom->x;
        MD_FLOAT* y_ptr           = x_ptr;
        MD_FLOAT* z_ptr           = x_ptr;
        MD_FLOAT* restrict fx_ptr = atom->fx;
        MD_FLOAT* fy_ptr          = fx_ptr;
        MD_FLOAT* fz_ptr          = fx_ptr;
#else
        MD_FLOAT* restrict x_ptr  = atom->x;
        MD_FLOAT* restrict y_ptr  = atom->y;
        MD_FLOAT* restrict z_ptr  = atom->z;
        MD_FLOAT* restrict fx_ptr = atom->fx;
        MD_FLOAT* restrict fy_ptr = atom->fy;
        MD_FLOAT* restrict fz_ptr = atom->fz;
#endif
        MD_FLOAT* restrict sqrt_epsilon_ptr = atom->sqrt_epsilon;
        MD_FLOAT* restrict sigma3_ptr       = atom->sigma3;
        int* restrict neighbors_ptr         = neighbor->neighbors;
        int* restrict numneigh_inner_ptr    = neighbor->numneigh_inner;

        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < nlocal; i++) {
            int numneighs = numneigh_inner_ptr[i];
            MD_FLOAT xtmp = atom_access_x(x_ptr, i);
            MD_FLOAT ytmp = atom_access_y(y_ptr, i);
            MD_FLOAT ztmp = atom_access_z(z_ptr, i);
            MD_FLOAT fix  = 0;
            MD_FLOAT fiy  = 0;
            MD_FLOAT fiz  = 0;

#if LJ_COMB_RULE == LJ_COMB_GEOM
            const MD_FLOAT sqrt_eps_i = sqrt_epsilon_ptr[i];
            const MD_FLOAT sigma3_i   = sigma3_ptr[i];
#elif LJ_COMB_RULE == LJ_COMB_FULL
            const int type_i = atom->type[i];
#endif

#pragma omp simd reduction(+ : fix, fiy, fiz)
            for (int k = 0; k < numneighs; k++) {
                int j         = neighs(neighbors_ptr, i, k, nlocal, neighbor);
                MD_FLOAT delx = xtmp - atom_access_x(x_ptr, j);
                MD_FLOAT dely = ytmp - atom_access_y(y_ptr, j);
                MD_FLOAT delz = ztmp - atom_access_z(z_ptr, j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
                const MD_FLOAT sigma6  = sigma3_i * sigma3_ptr[j];
                const MD_FLOAT epsilon = sqrt_eps_i * sqrt_epsilon_ptr[j];
#elif LJ_COMB_RULE == LJ_COMB_FULL
                const int type_j          = atom->type[j];
                const int type_ij         = type_i * atom->ntypes + type_j;
                const MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
                const MD_FLOAT sigma6     = atom->sigma6[type_ij];
                const MD_FLOAT epsilon    = atom->epsilon[type_ij];
#endif

                if (rsq < cutforcesq) {
                    MD_FLOAT u = LJ_TABLE_COORD(rsq) * inv_h;
                    int m      = (int)u;
                    if (m > nm1) m = nm1;
                    MD_FLOAT eps       = u - (MD_FLOAT)m;
                    const MD_FLOAT* cc = &coeff[m * LJ_TABLE_STRIDE];
                    MD_FLOAT hrep  = cc[0] + eps * (cc[1] + eps * (cc[2] + eps * cc[3]));
                    MD_FLOAT gdisp = cc[4] + eps * (cc[5] + eps * (cc[6] + eps * cc[7]));
                    MD_FLOAT force = epsilon * sigma6 * (sigma6 * hrep + gdisp);
                    fix += delx * force;
                    fiy += dely * force;
                    fiz += delz * force;
                }
            }

            atom_access_x(fx_ptr, i) += fix;
            atom_access_y(fy_ptr, i) += fiy;
            atom_access_z(fz_ptr, i) += fiz;

            addStat(stats->total_force_neighs, numneighs);
            addStat(stats->total_force_iters,
                (numneighs + VECTOR_WIDTH - 1) / VECTOR_WIDTH);
        }

        LIKWID_MARKER_STOP("force");
    }

    double timeStop = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJTableFullNeigh end\n");
    return timeStop - timeStart;
}

double computeForceLJTableHalfNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJTableHalfNeigh begin\n");

    int nlocal = atom->Nlocal;
    int nghost = atom->Nghost;
    int* neighs;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif
    const int nm1         = ljtable.n - 1;
    const MD_FLOAT inv_h  = ljtable.inv_h;
    const MD_FLOAT* coeff = ljtable.coeff;

    for (int i = 0; i < nlocal + nghost; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    double timeStart = getTimeStamp();

#pragma omp parallel
    {
        // See the matching comment in computeForceLJTableFullNeigh
#ifdef ATOM_POSITION_AOS
        MD_FLOAT* restrict x_ptr  = atom->x;
        MD_FLOAT* y_ptr           = x_ptr;
        MD_FLOAT* z_ptr           = x_ptr;
        MD_FLOAT* restrict fx_ptr = atom->fx;
        MD_FLOAT* fy_ptr          = fx_ptr;
        MD_FLOAT* fz_ptr          = fx_ptr;
#else
        MD_FLOAT* restrict x_ptr  = atom->x;
        MD_FLOAT* restrict y_ptr  = atom->y;
        MD_FLOAT* restrict z_ptr  = atom->z;
        MD_FLOAT* restrict fx_ptr = atom->fx;
        MD_FLOAT* restrict fy_ptr = atom->fy;
        MD_FLOAT* restrict fz_ptr = atom->fz;
#endif
        MD_FLOAT* restrict sqrt_epsilon_ptr = atom->sqrt_epsilon;
        MD_FLOAT* restrict sigma3_ptr       = atom->sigma3;
        int* restrict neighbors_ptr         = neighbor->neighbors;
        int* restrict numneigh_inner_ptr    = neighbor->numneigh_inner;

        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < nlocal; i++) {
            int numneighs = numneigh_inner_ptr[i];
            MD_FLOAT xtmp = atom_access_x(x_ptr, i);
            MD_FLOAT ytmp = atom_access_y(y_ptr, i);
            MD_FLOAT ztmp = atom_access_z(z_ptr, i);
            MD_FLOAT fix  = 0;
            MD_FLOAT fiy  = 0;
            MD_FLOAT fiz  = 0;

#if LJ_COMB_RULE == LJ_COMB_GEOM
            const MD_FLOAT sqrt_eps_i = sqrt_epsilon_ptr[i];
            const MD_FLOAT sigma3_i   = sigma3_ptr[i];
#elif LJ_COMB_RULE == LJ_COMB_FULL
            const int type_i = atom->type[i];
#endif

#pragma omp simd reduction(+ : fix, fiy, fiz)
            for (int k = 0; k < numneighs; k++) {
                int j         = neighs(neighbors_ptr, i, k, nlocal, neighbor);
                MD_FLOAT delx = xtmp - atom_access_x(x_ptr, j);
                MD_FLOAT dely = ytmp - atom_access_y(y_ptr, j);
                MD_FLOAT delz = ztmp - atom_access_z(z_ptr, j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
                const MD_FLOAT sigma6  = sigma3_i * sigma3_ptr[j];
                const MD_FLOAT epsilon = sqrt_eps_i * sqrt_epsilon_ptr[j];
#elif LJ_COMB_RULE == LJ_COMB_FULL
                const int type_j          = atom->type[j];
                const int type_ij         = type_i * atom->ntypes + type_j;
                const MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
                const MD_FLOAT sigma6     = atom->sigma6[type_ij];
                const MD_FLOAT epsilon    = atom->epsilon[type_ij];
#endif

                if (rsq < cutforcesq) {
                    MD_FLOAT u = LJ_TABLE_COORD(rsq) * inv_h;
                    int m      = (int)u;
                    if (m > nm1) m = nm1;
                    MD_FLOAT eps       = u - (MD_FLOAT)m;
                    const MD_FLOAT* cc = &coeff[m * LJ_TABLE_STRIDE];
                    MD_FLOAT hrep  = cc[0] + eps * (cc[1] + eps * (cc[2] + eps * cc[3]));
                    MD_FLOAT gdisp = cc[4] + eps * (cc[5] + eps * (cc[6] + eps * cc[7]));
                    MD_FLOAT force = epsilon * sigma6 * (sigma6 * hrep + gdisp);
                    fix += delx * force;
                    fiy += dely * force;
                    fiz += delz * force;

                    if ((param->half_neigh && j < nlocal) || param->method) {
                        atom_access_x(fx_ptr, j) -= delx * force;
                        atom_access_y(fy_ptr, j) -= dely * force;
                        atom_access_z(fz_ptr, j) -= delz * force;
                    }
                }
            }

            atom_access_x(fx_ptr, i) += fix;
            atom_access_y(fy_ptr, i) += fiy;
            atom_access_z(fz_ptr, i) += fiz;

            addStat(stats->total_force_neighs, numneighs);
            addStat(stats->total_force_iters,
                (numneighs + VECTOR_WIDTH - 1) / VECTOR_WIDTH);
        }

        LIKWID_MARKER_STOP("force");
    }

    double timeStop = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJTableHalfNeigh end\n");
    return timeStop - timeStart;
}
