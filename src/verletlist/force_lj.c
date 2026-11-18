/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <energy.h>
#include <kernel_variant.h>
#include <likwid-marker.h>
#include <neighbor.h>
#include <pair_scalar.h>
#include <parameter.h>
#include <stats.h>
#include <thread_forces.h>
#include <timing.h>
#include <util.h>

void computeForceGhostShell(Parameter*, Atom*, Neighbor*);

#if defined(__GNUC__) || defined(__clang__)
#define VL_ALWAYS_INLINE static inline __attribute__((always_inline))
#else
#define VL_ALWAYS_INLINE static inline
#endif

// GCC (tested 13.3 on SVE2) does not vectorize the full-rule neighbor loop when it honors
// "#pragma omp simd" (-fopenmp or -fopenmp-simd): the conditional sigma6/epsilon loads
// are then "not supported" as masked gathers and the loop stays scalar, with a
// data-dependent branch. Without the pragma GCC vectorizes it with SVE gathers. Clang
// needs the pragma to vectorize it.
#if defined(__GNUC__) && !defined(__clang__) && LJ_COMB_RULE == LJ_COMB_FULL
#define MD_OMP_SIMD_REDUCTION_F
#else
#define MD_OMP_SIMD_REDUCTION_F _Pragma("omp simd reduction(+ : fix, fiy, fiz, vlj, vc)")
#endif

// Read-only data of neighborLoop(), copied into restrict locals for vectorization
typedef struct {
    const MD_FLOAT* x;
    const MD_FLOAT* y;
    const MD_FLOAT* z;
    const MD_FLOAT* sqrt_epsilon;
    const MD_FLOAT* sigma3;
    const int* neighbors;
    const int* numneigh;
    const MD_FLOAT* q;
    const int* mol;
    CoulombScalar cs;
    MD_FLOAT cutforcesq; // single and geometric combination rules
    MD_FLOAT sigma6;     // single combination rule
    MD_FLOAT epsilon;
} LoopData;

static inline LoopData loopData(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    LoopData d;
    d.x = atom->x;
#ifdef ATOM_POSITION_AOS
    d.y = atom->x;
    d.z = atom->x;
#else
    d.y = atom->y;
    d.z = atom->z;
#endif
    d.sqrt_epsilon = atom->sqrt_epsilon;
    d.sigma3       = atom->sigma3;
    d.neighbors    = neighbor->neighbors;
    d.numneigh     = neighbor->numneigh_inner;
    d.q            = atom->q;
    d.mol          = atom->mol;
    d.cs           = coulombScalarInit(param);
    d.cutforcesq   = param->cutforce * param->cutforce;
    d.sigma6       = param->sigma6;
    d.epsilon      = param->epsilon;
    return d;
}

// Neighbor loop of atom i; energy, excl and reaction are constants at the call sites
VL_ALWAYS_INLINE void neighborLoop(const int energy,
    const int excl,
    const int reaction,
    const LoopData* d,
    const Atom* atom,
    const Neighbor* neighbor,
    const int i,
    const int nlocal,
    const int all_ghosts,
    MD_FLOAT* restrict fx_ptr,
    MD_FLOAT* restrict fy_ptr,
    MD_FLOAT* restrict fz_ptr,
    MD_FLOAT* fi,
    double* elj,
    double* ec,
    Stats* stats)
{
    const MD_FLOAT* restrict x_ptr            = d->x;
    const MD_FLOAT* restrict y_ptr            = d->y;
    const MD_FLOAT* restrict z_ptr            = d->z;
    const MD_FLOAT* restrict sqrt_epsilon_ptr = d->sqrt_epsilon;
    const MD_FLOAT* restrict sigma3_ptr       = d->sigma3;
    const int* restrict neighbors_ptr         = d->neighbors;
    const MD_FLOAT* restrict q_ptr            = d->q;
    const int* restrict mol_ptr               = d->mol;
    const CoulombScalar cs                    = d->cs;
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    const MD_FLOAT cutforcesq = d->cutforcesq;
    const MD_FLOAT sigma6     = d->sigma6;
    const MD_FLOAT epsilon    = d->epsilon;
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    const MD_FLOAT cutforcesq = d->cutforcesq;
#endif
    const int numneighs = d->numneigh[i];
    MD_FLOAT xtmp       = atom_access_x(x_ptr, i);
    MD_FLOAT ytmp       = atom_access_y(y_ptr, i);
    MD_FLOAT ztmp       = atom_access_z(z_ptr, i);
    MD_FLOAT fix        = 0;
    MD_FLOAT fiy        = 0;
    MD_FLOAT fiz        = 0;
    MD_FLOAT vlj        = 0;
    MD_FLOAT vc         = 0;

#if LJ_COMB_RULE == LJ_COMB_GEOM
    const MD_FLOAT sqrt_eps_i = sqrt_epsilon_ptr[i];
    const MD_FLOAT sigma3_i   = sigma3_ptr[i];
#elif LJ_COMB_RULE == LJ_COMB_FULL
    const int type_i = atom->type[i];
#endif
    const MD_FLOAT qi = cs.epsfac * q_ptr[i];
    const int mol_i   = mol_ptr[i];

    MD_OMP_SIMD_REDUCTION_F
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
            MD_FLOAT force = pairForceScalar(rsq,
                sigma6,
                epsilon,
                qi,
                q_ptr[j],
                mol_i,
                mol_ptr[j],
                excl,
                cs);
            fix += delx * force;
            fiy += dely * force;
            fiz += delz * force;
            const int react = reaction && (j < nlocal || all_ghosts);
            if (react) {
                atom_access_x(fx_ptr, j) -= delx * force;
                atom_access_y(fy_ptr, j) -= dely * force;
                atom_access_z(fz_ptr, j) -= delz * force;
            }
            if (energy) {
                MD_FLOAT pvlj, pvc;
                pairEnergyScalar(rsq,
                    sigma6,
                    epsilon,
                    qi,
                    q_ptr[j],
                    mol_i,
                    mol_ptr[j],
                    excl,
                    cs,
                    &pvlj,
                    &pvc);
                const MD_FLOAT w = react ? (MD_FLOAT)1.0 : (MD_FLOAT)0.5;
                vlj += w * pvlj;
                vc += w * pvc;
            }
        }
    }

    fi[0] = fix;
    fi[1] = fiy;
    fi[2] = fiz;
    if (energy) {
        *elj += vlj;
        *ec += vc;
    }
#ifdef COMPUTE_STATS
    // Separate loop: a shared counter in the SIMD loop above loses increments
    long long int within = 0;
    for (int k = 0; k < numneighs; k++) {
        const int j         = neighs(neighbors_ptr, i, k, nlocal, neighbor);
        const MD_FLOAT delx = xtmp - atom_access_x(x_ptr, j);
        const MD_FLOAT dely = ytmp - atom_access_y(y_ptr, j);
        const MD_FLOAT delz = ztmp - atom_access_z(z_ptr, j);
#if LJ_COMB_RULE == LJ_COMB_FULL
        const MD_FLOAT cutforcesq =
            atom->cutforcesq[type_i * atom->ntypes + atom->type[j]];
#endif
        within += (delx * delx + dely * dely + delz * delz < cutforcesq);
    }
    addStat(stats->atoms_within_cutoff, within);
    addStat(stats->atoms_outside_cutoff, numneighs - within);
    if (!reaction && numneighs % VECTOR_WIDTH > 0) {
        addStat(stats->atoms_outside_cutoff, VECTOR_WIDTH - (numneighs % VECTOR_WIDTH));
    }
#endif
    addStat(stats->total_force_neighs, numneighs);
    addStat(stats->total_force_iters, (numneighs + VECTOR_WIDTH - 1) / VECTOR_WIDTH);
    (void)atom;
    (void)stats;
    (void)sqrt_epsilon_ptr;
    (void)sigma3_ptr;
}

double KERNEL_FN(computeForceLJFullNeigh)(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJFullNeigh begin\n");
    const int nlocal  = atom->Nlocal;
    const int energy  = param->compute_energy;
    const LoopData ld = loopData(param, atom, neighbor);
    double elj = 0.0, ec = 0.0;

    for (int i = 0; i < nlocal; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    double timeStart = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime) reduction(+ : elj, ec)
        for (int i = 0; i < nlocal; i++) {
            MD_FLOAT fi[3];
            // KERNEL_EXCL: copy of this file for systems with molecules
            if (energy) {
                neighborLoop(1,
                    KERNEL_EXCL,
                    0,
                    &ld,
                    atom,
                    neighbor,
                    i,
                    nlocal,
                    0,
                    NULL,
                    NULL,
                    NULL,
                    fi,
                    &elj,
                    &ec,
                    stats);
            } else {
                neighborLoop(0,
                    KERNEL_EXCL,
                    0,
                    &ld,
                    atom,
                    neighbor,
                    i,
                    nlocal,
                    0,
                    NULL,
                    NULL,
                    NULL,
                    fi,
                    &elj,
                    &ec,
                    stats);
            }
            atom_fx(i) += fi[0];
            atom_fy(i) += fi[1];
            atom_fz(i) += fi[2];
        }
        LIKWID_MARKER_STOP("force");
    }

    if (energy) {
        storeEnergies(param,
            elj,
            ec,
            localChargeSquaredSum(nlocal, atom->type, atom->q_per_type));
    }

    double timeStop = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJFullNeigh end\n");
    return timeStop - timeStart;
}

double KERNEL_FN(computeForceLJHalfNeigh)(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJHalfNeigh begin\n");
    const int nlocal  = atom->Nlocal;
    const int nghost  = atom->Nghost;
    const int energy  = param->compute_energy;
    const LoopData ld = loopData(param, atom, neighbor);
    // Ghost reaction forces only with a shell method (else both processes compute them)
    const int all_ghosts = param->method != 0;
    double elj = 0.0, ec = 0.0;

    for (int i = 0; i < nlocal + nghost; i++) {
        atom_fx(i) = 0.0;
        atom_fy(i) = 0.0;
        atom_fz(i) = 0.0;
    }

    if (energy) {
        storeEnergies(param, 0.0, 0.0, 0.0);
    }

    const int nf     = threadForcesReserve(atom);
    double timeStart = getTimeStamp();

#pragma omp parallel
    {
        // Reaction forces go to the force arrays of this thread, see thread_forces.h
        const ThreadForces tf = threadForces(atom, nf);
#ifdef ATOM_POSITION_AOS
        MD_FLOAT* fx_ptr = tf.fx;
        MD_FLOAT* fy_ptr = tf.fx;
        MD_FLOAT* fz_ptr = tf.fx;
#else
        MD_FLOAT* fx_ptr = tf.fx;
        MD_FLOAT* fy_ptr = tf.fy;
        MD_FLOAT* fz_ptr = tf.fz;
#endif
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime) reduction(+ : elj, ec)
        for (int i = 0; i < nlocal; i++) {
            MD_FLOAT fi[3];
            // KERNEL_EXCL: copy of this file for systems with molecules
            if (energy) {
                neighborLoop(1,
                    KERNEL_EXCL,
                    1,
                    &ld,
                    atom,
                    neighbor,
                    i,
                    nlocal,
                    all_ghosts,
                    fx_ptr,
                    fy_ptr,
                    fz_ptr,
                    fi,
                    &elj,
                    &ec,
                    stats);
            } else {
                neighborLoop(0,
                    KERNEL_EXCL,
                    1,
                    &ld,
                    atom,
                    neighbor,
                    i,
                    nlocal,
                    all_ghosts,
                    fx_ptr,
                    fy_ptr,
                    fz_ptr,
                    fi,
                    &elj,
                    &ec,
                    stats);
            }
            atom_access_x(fx_ptr, i) += fi[0];
            atom_access_y(fy_ptr, i) += fi[1];
            atom_access_z(fz_ptr, i) += fi[2];
        }

        if (param->method == eightShell) {
            // The ghost-ghost pairs of the eight-shell method, by one thread into the
            // atom forces (the implicit barrier of the loop above precedes it)
#pragma omp single
            computeForceGhostShell(param, atom, neighbor);
        }

        threadForcesReduce(atom, nf);
        LIKWID_MARKER_STOP("force");
    }

    if (energy) {
        const double shell_lj   = param->energy_lj;
        const double shell_coul = param->energy_coul;
        storeEnergies(param,
            elj,
            ec,
            localChargeSquaredSum(nlocal, atom->type, atom->q_per_type));
        param->energy_lj += shell_lj;
        param->energy_coul += shell_coul;
    }

    double timeStop = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJHalfNeigh end\n");
    return timeStop - timeStart;
}

#ifndef KERNEL_EXCLUSIONS
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
    const CoulombScalar cs = coulombScalarInit(param);
    const int excl         = atom->molsize > 1; // molecular exclusions

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
        const MD_FLOAT qi = cs.epsfac * atom->q[iatom];
        const int mol_i   = atom->mol[iatom];

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
                MD_FLOAT force = pairForceScalar(rsq,
                    sigma6,
                    epsilon,
                    qi,
                    atom->q[jatom],
                    mol_i,
                    atom->mol[jatom],
                    excl,
                    cs);
                if (param->compute_energy) {
                    MD_FLOAT vlj, vc;
                    pairEnergyScalar(rsq,
                        sigma6,
                        epsilon,
                        qi,
                        atom->q[jatom],
                        mol_i,
                        atom->mol[jatom],
                        excl,
                        cs,
                        &vlj,
                        &vc);
                    param->energy_lj += vlj;
                    param->energy_coul += vc;
                }
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
#endif
