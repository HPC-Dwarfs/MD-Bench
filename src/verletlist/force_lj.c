/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <likwid-marker.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>

void computeForceGhostShell(Parameter*, Atom*, Neighbor*);

double computeForceLJFullNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJFullNeigh begin\n");

    int nlocal = atom->Nlocal;
    int* neighs;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif
    const MD_FLOAT num1  = 1.0;
    const MD_FLOAT num48 = 48.0;
    const MD_FLOAT num05 = 0.5;

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
#if LJ_COMB_RULE == LJ_COMB_GEOM
        MD_FLOAT* restrict sqrt_epsilon_ptr = atom->sqrt_epsilon;
        MD_FLOAT* restrict sigma3_ptr       = atom->sigma3;
#endif
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
                    MD_FLOAT sr2   = num1 / rsq;
                    MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
                    MD_FLOAT force = num48 * sr6 * (sr6 - num05) * sr2 * epsilon;
                    fix += delx * force;
                    fiy += dely * force;
                    fiz += delz * force;
#ifdef USE_REFERENCE_KERNEL
                    addStat(stats->atoms_within_cutoff, 1);
                } else {
                    addStat(stats->atoms_outside_cutoff, 1);
#endif
                }
            }

            atom_access_x(fx_ptr, i) += fix;
            atom_access_y(fy_ptr, i) += fiy;
            atom_access_z(fz_ptr, i) += fiz;

#ifdef USE_REFERENCE_KERNEL
            if (numneighs % VECTOR_WIDTH > 0) {
                addStat(stats->atoms_outside_cutoff,
                    VECTOR_WIDTH - (numneighs % VECTOR_WIDTH));
            }
#endif

            addStat(stats->total_force_neighs, numneighs);
            addStat(stats->total_force_iters,
                (numneighs + VECTOR_WIDTH - 1) / VECTOR_WIDTH);
        }

        LIKWID_MARKER_STOP("force");
    }

    double timeStop = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJFullNeigh end\n");
    return timeStop - timeStart;
}

double computeForceLJHalfNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJHalfNeigh begin\n");

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
    const MD_FLOAT num1  = 1.0;
    const MD_FLOAT num48 = 48.0;
    const MD_FLOAT num05 = 0.5;

    for (int i = 0; i < nlocal + nghost; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    double timeStart = getTimeStamp();

#pragma omp parallel
    {
        // See the matching comment in computeForceLJFullNeigh
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

// Pragma required to vectorize the inner loop
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
                    MD_FLOAT sr2   = num1 / rsq;
                    MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
                    MD_FLOAT force = num48 * sr6 * (sr6 - num05) * sr2 * epsilon;
                    fix += delx * force;
                    fiy += dely * force;
                    fiz += delz * force;

                    // Skip ghost atoms unless shell_method or half stencil needs them
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
        if (param->method == eightShell) computeForceGhostShell(param, atom, neighbor);
        LIKWID_MARKER_STOP("force");
    }

    double timeStop = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJHalfNeigh end\n");
    return timeStop - timeStart;
}

void computeForceGhostShell(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    DEBUG_MESSAGE("computeForceGhostShell begin\n");
    int Nshell = neighbor->Nshell;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
    MD_FLOAT sigma6     = param->sigma6;
    MD_FLOAT epsilon    = param->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif
    const MD_FLOAT num1  = 1.0;
    const MD_FLOAT num48 = 48.0;
    const MD_FLOAT num05 = 0.5;

    for (int i = 0; i < Nshell; i++) {
        int numneigh  = neighbor->numNeighShell[i];
        int iatom     = neighbor->listshell[i];
        MD_FLOAT xtmp = atom_x(iatom);
        MD_FLOAT ytmp = atom_y(iatom);
        MD_FLOAT ztmp = atom_z(iatom);
        MD_FLOAT fix  = 0;
        MD_FLOAT fiy  = 0;
        MD_FLOAT fiz  = 0;

#if LJ_COMB_RULE == LJ_COMB_GEOM
        const MD_FLOAT sqrt_eps_i = atom->sqrt_epsilon[iatom];
        const MD_FLOAT sigma3_i   = atom->sigma3[iatom];
#elif LJ_COMB_RULE == LJ_COMB_FULL
        const int type_i = atom->type[iatom];
#endif

        for (int k = 0; k < numneigh; k++) {
            int jatom     = neighshell(neighbor->neighshell, i, k, neighbor->maxneighs);
            MD_FLOAT delx = xtmp - atom_x(jatom);
            MD_FLOAT dely = ytmp - atom_y(jatom);
            MD_FLOAT delz = ztmp - atom_z(jatom);
            MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
            const MD_FLOAT sigma6  = sigma3_i * atom->sigma3[jatom];
            const MD_FLOAT epsilon = sqrt_eps_i * atom->sqrt_epsilon[jatom];
#elif LJ_COMB_RULE == LJ_COMB_FULL
            const int type_j          = atom->type[jatom];
            const int type_ij         = type_i * atom->ntypes + type_j;
            const MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
            const MD_FLOAT sigma6     = atom->sigma6[type_ij];
            const MD_FLOAT epsilon    = atom->epsilon[type_ij];
#endif

            if (rsq < cutforcesq) {
                MD_FLOAT sr2   = num1 / rsq;
                MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
                MD_FLOAT force = num48 * sr6 * (sr6 - num05) * sr2 * epsilon;
                fix += delx * force;
                fiy += dely * force;
                fiz += delz * force;

                atom_fx(jatom) -= delx * force;
                atom_fy(jatom) -= dely * force;
                atom_fz(jatom) -= delz * force;
            }
        }

        atom_fx(iatom) += fix;
        atom_fy(iatom) += fiy;
        atom_fz(iatom) += fiz;
    }

    DEBUG_MESSAGE("computeForceGhostShell end\n");
}
