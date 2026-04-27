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
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < nlocal; i++) {
            int numneighs = neighbor->numneigh[i];
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            MD_FLOAT fix  = 0;
            MD_FLOAT fiy  = 0;
            MD_FLOAT fiz  = 0;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            const int type_i = atom->type[i];
#endif

            for (int k = 0; k < numneighs; k++) {
                int j = neighs(neighbor->neighbors, i, k, nlocal, neighbor->maxneighs);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
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

            atom_fx(i) += fix;
            atom_fy(i) += fiy;
            atom_fz(i) += fiz;

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
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int i = 0; i < nlocal; i++) {
            int numneighs = neighbor->numneigh[i];
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            MD_FLOAT fix  = 0;
            MD_FLOAT fiy  = 0;
            MD_FLOAT fiz  = 0;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            const int type_i = atom->type[i];
#endif

// Pragma required to vectorize the inner loop
#pragma omp simd reduction(+ : fix, fiy, fiz)
            for (int k = 0; k < numneighs; k++) {
                int j = neighs(neighbor->neighbors, i, k, nlocal, neighbor->maxneighs);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
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

                    // We do not need to update forces for ghost atoms
                    // We need to update forces for ghost atoms if shell_method  or half
                    // stencil is required
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

#if LJ_COMB_RULE != LJ_COMB_SINGLE
        const int type_i = atom->type[iatom];
#endif

        for (int k = 0; k < numneigh; k++) {
            int jatom = neighs(neighbor->neighshell, i, k, Nshell, neighbor->maxneighs);
            MD_FLOAT delx = xtmp - atom_x(jatom);
            MD_FLOAT dely = ytmp - atom_y(jatom);
            MD_FLOAT delz = ztmp - atom_z(jatom);
            MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
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
