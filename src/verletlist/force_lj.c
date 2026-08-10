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

#ifdef NBLIST_PAIRLIST
#include <allocate.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

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
#elif LJ_COMB_RULE == LJ_COMB_FULL
            const int type_i = atom->type[i];
#endif

            for (int k = 0; k < numneighs; k++) {
                int j         = neighs(neighbor->neighbors, i, k, nlocal, neighbor);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
                const MD_FLOAT sigma6  = sigma3_i * atom->sigma3[j];
                const MD_FLOAT epsilon = sqrt_eps_i * atom->sqrt_epsilon[j];
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
#elif LJ_COMB_RULE == LJ_COMB_FULL
            const int type_i = atom->type[i];
#endif

// Pragma required to vectorize the inner loop
#pragma omp simd reduction(+ : fix, fiy, fiz)
            for (int k = 0; k < numneighs; k++) {
                int j         = neighs(neighbor->neighbors, i, k, nlocal, neighbor);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
                const MD_FLOAT sigma6  = sigma3_i * atom->sigma3[j];
                const MD_FLOAT epsilon = sqrt_eps_i * atom->sqrt_epsilon[j];
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
        if (param->method == eightShell) {
#pragma omp single
            computeForceGhostShell(param, atom, neighbor);
        }
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

#ifdef NBLIST_PAIRLIST
// Per-thread scratch force buffer for the pair-list kernels: nthreads
// contiguous [nall*3] slices. Grown lazily (freed+reallocated, contents are
// never preserved across a resize since every step re-zeroes it) instead of
// being threaded through growAtom/Atom, since it's purely a kernel-internal
// reduction scratch, not atom state.
static MD_FLOAT* pairlist_fbuf     = NULL;
static size_t pairlist_fbuf_cap    = 0;

static MD_FLOAT* getPairlistForceBuffer(int nall, int nthreads)
{
    size_t needed = (size_t)nthreads * (size_t)nall * 3;
    if (needed > pairlist_fbuf_cap) {
        if (pairlist_fbuf) free(pairlist_fbuf);
        size_t alloc_size = needed > 0 ? needed : 1;
        pairlist_fbuf      = (MD_FLOAT*)allocate(ALIGNMENT, alloc_size * sizeof(MD_FLOAT));
        pairlist_fbuf_cap  = needed;
    }
    return pairlist_fbuf;
}

// Reference (scalar) pair-list kernel: iterates the flat (i,j) pair list
// built by flattenPairList() instead of a per-i neighbor list. Because
// different threads may process pairs touching the same i or j, each thread
// accumulates into its own private force-buffer slice (no cross-thread
// write races regardless of which atoms a pair touches), then a final pass
// reduces all thread slices into the real force array.
double computeForceLJPairListRef(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJPairListRef begin\n");

    int nlocal = atom->Nlocal;
    int nghost = atom->Nghost;
    int nall   = nlocal + nghost;
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

    for (int i = 0; i < nall; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    MD_FLOAT* fbuf  = getPairlistForceBuffer(MAX(1, nall), nthreads);
    int npairs      = neighbor->npairs;

    double timeStart = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        MD_FLOAT* tf = fbuf + (size_t)tid * nall * 3;
        for (int k = 0; k < nall * 3; k++) {
            tf[k] = 0.0;
        }

#pragma omp for schedule(static)
        for (int p = 0; p < npairs; p++) {
            int i = neighbor->pair_i[p];
            int j = neighbor->pair_j[p];

            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            MD_FLOAT delx = xtmp - atom_x(j);
            MD_FLOAT dely = ytmp - atom_y(j);
            MD_FLOAT delz = ztmp - atom_z(j);
            MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE == LJ_COMB_GEOM
            const MD_FLOAT sigma6  = atom->sigma3[i] * atom->sigma3[j];
            const MD_FLOAT epsilon = atom->sqrt_epsilon[i] * atom->sqrt_epsilon[j];
#elif LJ_COMB_RULE == LJ_COMB_NONE
            const int type_i          = atom->type[i];
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
                MD_FLOAT fx    = delx * force;
                MD_FLOAT fy    = dely * force;
                MD_FLOAT fz    = delz * force;

                tf[i * 3 + 0] += fx;
                tf[i * 3 + 1] += fy;
                tf[i * 3 + 2] += fz;

                if ((param->half_neigh && j < nlocal) || param->method) {
                    tf[j * 3 + 0] -= fx;
                    tf[j * 3 + 1] -= fy;
                    tf[j * 3 + 2] -= fz;
                }
            }
        }

#pragma omp for schedule(static)
        for (int i = 0; i < nall; i++) {
            MD_FLOAT sx = 0.0, sy = 0.0, sz = 0.0;
            for (int t = 0; t < nthreads; t++) {
                MD_FLOAT* tfi = fbuf + (size_t)t * nall * 3;
                sx += tfi[i * 3 + 0];
                sy += tfi[i * 3 + 1];
                sz += tfi[i * 3 + 2];
            }
            atom_fx(i) += sx;
            atom_fy(i) += sy;
            atom_fz(i) += sz;
        }

        LIKWID_MARKER_STOP("force");
    }

    addStat(stats->total_force_neighs, npairs);

    if (param->method == eightShell) {
        computeForceGhostShell(param, atom, neighbor);
    }

    double timeStop = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJPairListRef end\n");
    return timeStop - timeStart;
}
#endif
