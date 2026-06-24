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
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < nlocal; i++) {
            int numneighs = neighbor->numneigh_inner[i];
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            MD_FLOAT fix  = 0;
            MD_FLOAT fiy  = 0;
            MD_FLOAT fiz  = 0;

#if LJ_COMB_RULE == LJ_COMB_GEOM
            const MD_FLOAT sqrt_eps_i = atom->sqrt_epsilon[i];
            const MD_FLOAT sigma3_i   = atom->sigma3[i];
#elif LJ_COMB_RULE == LJ_COMB_NONE
            const int type_i = atom->type[i];
#endif

            for (int k = 0; k < numneighs; k++) {
                int j = neighs(neighbor->neighbors, i, k, nlocal, neighbor->maxneighs);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
                const MD_FLOAT sigma6  = sigma3_i * atom->sigma3[j];
                const MD_FLOAT epsilon = sqrt_eps_i * atom->sqrt_epsilon[j];
#elif LJ_COMB_RULE == LJ_COMB_NONE
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

            atom_fx(i) += fix;
            atom_fy(i) += fiy;
            atom_fz(i) += fiz;

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
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < nlocal; i++) {
            int numneighs = neighbor->numneigh_inner[i];
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            MD_FLOAT fix  = 0;
            MD_FLOAT fiy  = 0;
            MD_FLOAT fiz  = 0;

#if LJ_COMB_RULE == LJ_COMB_GEOM
            const MD_FLOAT sqrt_eps_i = atom->sqrt_epsilon[i];
            const MD_FLOAT sigma3_i   = atom->sigma3[i];
#elif LJ_COMB_RULE == LJ_COMB_NONE
            const int type_i = atom->type[i];
#endif

            for (int k = 0; k < numneighs; k++) {
                int j = neighs(neighbor->neighbors, i, k, nlocal, neighbor->maxneighs);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
                const MD_FLOAT sigma6  = sigma3_i * atom->sigma3[j];
                const MD_FLOAT epsilon = sqrt_eps_i * atom->sqrt_epsilon[j];
#elif LJ_COMB_RULE == LJ_COMB_NONE
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
                        atom_fx(j) -= delx * force;
                        atom_fy(j) -= dely * force;
                        atom_fz(j) -= delz * force;
                    }
                }
            }

            atom_fx(i) += fix;
            atom_fy(i) += fiy;
            atom_fz(i) += fiz;

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
