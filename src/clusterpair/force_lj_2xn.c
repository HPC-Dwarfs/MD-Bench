/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* SIMD 2xN cluster-pair kernels: 2 i-atoms against N = SIMD width j-atoms, so an
 * i-cluster of 2 atoms is 2 rows (the 4xN kernels with half the rows, for narrow
 * registers of NEON and SVE). Half and full neighbor lists. The pair math is shared with
 * the other SIMD kernels, see pair_rows.h. */
#include <stdio.h>

#include <atom.h>
#include <energy.h>
#include <force.h>
#include <kernel_variant.h>
#include <likwid-marker.h>
#include <math.h>
#include <neighbor.h>
#include <pair_rows.h>
#include <parameter.h>
#include <simd.h>
#include <stats.h>
#include <thread_forces.h>
#include <timing.h>
#include <util.h>

// Only the copies of the selected kernel have code
#if !defined(USE_REFERENCE_KERNEL) && defined(CLUSTERPAIR_KERNEL_2XN)

void RULE_VARIANT_NAME(computeForceGhostShell)(Parameter*, Atom*, Neighbor*);

static inline IRows load2xnIRows(const Atom* atom, const int ci, const PairConstants* pc)
{
    IRows ir;
    const MD_FLOAT* ci_x = &atom->cl_x[CI_VECTOR3_BASE_INDEX(ci)];
    const int sca        = CI_SCALAR_BASE_INDEX(ci);

    for (int r = 0; r < CLUSTER_M; r++) {
        ir.x.v[r] = simd_real_broadcast(ci_x[CL_X_INDEX_3D(r)]);
        ir.y.v[r] = simd_real_broadcast(ci_x[CL_Y_INDEX_3D(r)]);
        ir.z.v[r] = simd_real_broadcast(ci_x[CL_Z_INDEX_3D(r)]);
#if COULOMB_RULE != COULOMB_NONE
        ir.q.v[r] = simd_real_mul(pc->cs.epsfac,
            simd_real_broadcast(atom->cl_q[sca + r]));
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
        ir.a6.v[r]  = simd_real_broadcast(atom->cl_sqrt6c6[sca + r]);
        ir.b12.v[r] = simd_real_broadcast(atom->cl_sqrt12c12[sca + r]);
#elif LJ_COMB_RULE == LJ_COMB_FULL
        ir.type.v[r] = simd_i32_broadcast(atom->cl_t[sca + r] * atom->ntypes);
#endif
#if KERNEL_EXCL
        ir.mol.v[r] = simd_i32_broadcast(atom->cl_mol[sca + r]);
#endif
        ir.id.v[r] = simd_i32_broadcast(atom->cl_id[sca + r]);
    }

    (void)sca;
    (void)pc;
    return ir;
}

static inline JData load2xnJData(const Atom* atom, const int cj)
{
    JData jd;
    MD_FLOAT* cj_x = &atom->cl_x[CJ_VECTOR3_BASE_INDEX(cj)];
    const int sca  = CJ_SCALAR_BASE_INDEX(cj);

    jd.x = simd_real_load(&cj_x[CL_X_INDEX_3D(0)]);
    jd.y = simd_real_load(&cj_x[CL_Y_INDEX_3D(0)]);
    jd.z = simd_real_load(&cj_x[CL_Z_INDEX_3D(0)]);
#if COULOMB_RULE != COULOMB_NONE
    jd.q = simd_real_load(&atom->cl_q[sca]);
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
    jd.a6  = simd_real_load(&atom->cl_sqrt6c6[sca]);
    jd.b12 = simd_real_load(&atom->cl_sqrt12c12[sca]);
#elif LJ_COMB_RULE == LJ_COMB_FULL
    jd.type = simd_i32_load(&atom->cl_t[sca]);
#endif
#if KERNEL_EXCL
    jd.mol = simd_i32_load(&atom->cl_mol[sca]);
#endif
    jd.id = simd_i32_load(&atom->cl_id[sca]);

    (void)sca;
    return jd;
}

/* Rows of one entry of the diagonal mask table (the entries are selected by which
 * of the i-cluster halves is the j-cluster). */
#define DIAG_MASK_ROWS CLUSTER_M

/* Masks of the pairs to skip in a pair of clusters that is, or overlaps with, the
 * cluster itself, see the 4xN kernels */
static inline RowMasks load2xnDiagMasks(
    const unsigned int* table, const int ci, const int cj)
{
    RowMasks diag;
    unsigned int idx;
#if CLUSTER_M == CLUSTER_N
    idx = (unsigned int)(cj == CJ0_FROM_CI(ci)) * DIAG_MASK_ROWS;
#elif CLUSTER_M < CLUSTER_N
    idx = (unsigned int)((cj << 1) + 0 == ci) * (2 * DIAG_MASK_ROWS) +
          (unsigned int)((cj << 1) + 1 == ci) * DIAG_MASK_ROWS;
#else
    idx = (unsigned int)(cj == CJ0_FROM_CI(ci)) * (2 * DIAG_MASK_ROWS) +
          (unsigned int)(cj == CJ1_FROM_CI(ci)) * DIAG_MASK_ROWS;
#endif
    for (int r = 0; r < CLUSTER_M; r++) {
        diag.v[r] = simd_mask_from_u32(table[idx + r]);
    }
    return diag;
}

/* Reaction force on the atoms of the j-cluster */
static inline void decrementJForces(MD_FLOAT* cj_f, const PairForces* p)
{
    const MD_SIMD_FLOAT tx = rowsSum(CLUSTER_M, p->tx);
    const MD_SIMD_FLOAT ty = rowsSum(CLUSTER_M, p->ty);
    const MD_SIMD_FLOAT tz = rowsSum(CLUSTER_M, p->tz);
    simd_real_store(&cj_f[CL_X_INDEX_3D(0)],
        simd_real_sub(simd_real_load(&cj_f[CL_X_INDEX_3D(0)]), tx));
    simd_real_store(&cj_f[CL_Y_INDEX_3D(0)],
        simd_real_sub(simd_real_load(&cj_f[CL_Y_INDEX_3D(0)]), ty));
    simd_real_store(&cj_f[CL_Z_INDEX_3D(0)],
        simd_real_sub(simd_real_load(&cj_f[CL_Z_INDEX_3D(0)]), tz));
}

static inline void incrementIForces(
    MD_FLOAT* ci_f, const Rows* fx, const Rows* fy, const Rows* fz)
{
    simd_real_incr_reduced_sum_j2(&ci_f[CL_X_INDEX_3D(0)], fx->v[0], fx->v[1]);
    simd_real_incr_reduced_sum_j2(&ci_f[CL_Y_INDEX_3D(0)], fy->v[0], fy->v[1]);
    simd_real_incr_reduced_sum_j2(&ci_f[CL_Z_INDEX_3D(0)], fz->v[0], fz->v[1]);
}

// One i-cluster of the half list: nrlj LJ rows, energy flag; constants at the call sites
MD_ALWAYS_INLINE void iClusterHalf2xn(Parameter* param,
    Atom* atom,
    Neighbor* neighbor,
    Stats* stats,
    const PairConstants* pc,
    MD_FLOAT* f,
    const int ci,
    const int nrlj,
    const int energy,
    double* ec,
    double* elj)
{
    const int nbM                    = atom->Nclusters_local;
    const int nlocal_clusters        = CJ1_FROM_CI(atom->Nclusters_local);
    const int numneighs              = neighbor->numneigh_inner[ci];
    const int numneighs_masked       = neighbor->numneigh_masked[ci];
    const int numneighs_inner_masked = neighbor->numneigh_inner_masked[ci];
    const int unmasked_inner_end     = numneighs_masked +
                                       (numneighs - numneighs_inner_masked);
    MD_FLOAT* ci_f                   = &f[CI_VECTOR3_BASE_INDEX(ci)];

    const IRows ir = load2xnIRows(atom, ci, pc);
    Rows fx        = rowsZero();
    Rows fy        = rowsZero();
    Rows fz        = rowsZero();
    // e_half: pairs also computed on the process of the ghost j-cluster, they count half
    RowEnergies e      = rowEnergiesZero();
    RowEnergies e_half = rowEnergiesZero();

    for (int k = 0; k < numneighs_inner_masked; k++) {
        prefetchAhead(atom,
            neighbor,
            ci,
            nbM,
            k,
            numneighs,
            numneighs_masked,
            numneighs_inner_masked,
            f);
        const int cj        = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
        const JData jd      = load2xnJData(atom, cj);
        const RowMasks diag = load2xnDiagMasks(atom->masks_2xn_hn, ci, cj);
        const PairForces p =
            pairForces(CLUSTER_M, nrlj, energy, &ir, &jd, &diag, KERNEL_EXCL, pc, atom);

        accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
        if (cj < nlocal_clusters || param->method) {
            decrementJForces(&f[CJ_VECTOR3_BASE_INDEX(cj)], &p);
            if (energy) accumulateEnergies(CLUSTER_M, &e, &p);
        } else if (energy) {
            accumulateEnergies(CLUSTER_M, &e_half, &p);
        }
    }

    for (int k = numneighs_masked; k < unmasked_inner_end; k++) {
        prefetchAhead(atom,
            neighbor,
            ci,
            nbM,
            k - numneighs_masked + numneighs_inner_masked,
            numneighs,
            numneighs_masked,
            numneighs_inner_masked,
            f);
        const int cj   = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
        const JData jd = load2xnJData(atom, cj);
        const PairForces p =
            pairForces(CLUSTER_M, nrlj, energy, &ir, &jd, NULL, 0, pc, atom);

        accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
        if (cj < nlocal_clusters || param->method) {
            decrementJForces(&f[CJ_VECTOR3_BASE_INDEX(cj)], &p);
            if (energy) accumulateEnergies(CLUSTER_M, &e, &p);
        } else if (energy) {
            accumulateEnergies(CLUSTER_M, &e_half, &p);
        }
    }

    incrementIForces(ci_f, &fx, &fy, &fz);
    if (energy) {
        reduceEnergies(CLUSTER_M, &e, 1.0, ec, elj);
        reduceEnergies(CLUSTER_M, &e_half, 0.5, ec, elj);
    }

    addStat(stats->calculated_forces, 1);
    addStat(stats->num_neighs, numneighs);
    addStat(stats->force_iters,
        (long long int)((double)numneighs * CLUSTER_M / CLUSTER_N));
}

double VARIANT_NAME(computeForceLJ2xnHalfNeigh)(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_2xn begin\n");
    const PairConstants pc = pairConstants(param);
    const int energy       = param->compute_energy;
    double ec = 0.0, elj = 0.0;
#if LJ_COMB_RULE == LJ_COMB_FULL
    setupTypePairC6C12(atom);
#endif

    const int nf = threadForcesReserve(atom);
    zeroForces(atom, 1);

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");
        MD_FLOAT* f = threadForces(atom, nf);

#pragma omp for schedule(runtime) reduction(+ : ec, elj)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            const int half_lj = param->half_lj && atom->iclusters[ci].half_lj;
            if (!energy) {
                if (!half_lj) {
                    iClusterHalf2xn(param,
                        atom,
                        neighbor,
                        stats,
                        &pc,
                        f,
                        ci,
                        CLUSTER_M,
                        0,
                        &ec,
                        &elj);
                } else {
                    iClusterHalf2xn(param,
                        atom,
                        neighbor,
                        stats,
                        &pc,
                        f,
                        ci,
                        CLUSTER_M / 2,
                        0,
                        &ec,
                        &elj);
                }
            } else {
                if (!half_lj) {
                    iClusterHalf2xn(param,
                        atom,
                        neighbor,
                        stats,
                        &pc,
                        f,
                        ci,
                        CLUSTER_M,
                        1,
                        &ec,
                        &elj);
                } else {
                    iClusterHalf2xn(param,
                        atom,
                        neighbor,
                        stats,
                        &pc,
                        f,
                        ci,
                        CLUSTER_M / 2,
                        1,
                        &ec,
                        &elj);
                }
            }
        }
        threadForcesReduce(atom, nf);
        LIKWID_MARKER_STOP("force");
    }

    if (energy) {
        storeEnergies(param,
            elj,
            ec,
            localChargeSquaredSum(atom->Nlocal, atom->type, atom->q_per_type));
    }

    // Single threaded: it writes cl_f of the shell clusters without a thread buffer
    if (param->method == eightShell) {
        RULE_VARIANT_NAME(computeForceGhostShell)(param, atom, neighbor);
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_2xn end\n");
    return E - S;
}

/* One i-cluster of the full-list kernel, see iClusterHalf2xn() */
MD_ALWAYS_INLINE void iClusterFull2xn(Atom* atom,
    Neighbor* neighbor,
    Stats* stats,
    const PairConstants* pc,
    const int ci,
    const int nrlj,
    const int energy,
    double* ec,
    double* elj)
{
    const int nbM                    = atom->Nclusters_local;
    const int numneighs              = neighbor->numneigh_inner[ci];
    const int numneighs_masked       = neighbor->numneigh_masked[ci];
    const int numneighs_inner_masked = neighbor->numneigh_inner_masked[ci];
    const int unmasked_inner_end     = numneighs_masked +
                                       (numneighs - numneighs_inner_masked);
    MD_FLOAT* ci_f                   = &atom->cl_f[CI_VECTOR3_BASE_INDEX(ci)];

    const IRows ir = load2xnIRows(atom, ci, pc);
    Rows fx        = rowsZero();
    Rows fy        = rowsZero();
    Rows fz        = rowsZero();
    RowEnergies e  = rowEnergiesZero();

    for (int k = 0; k < numneighs_inner_masked; k++) {
        prefetchAhead(atom,
            neighbor,
            ci,
            nbM,
            k,
            numneighs,
            numneighs_masked,
            numneighs_inner_masked,
            NULL);
        const int cj        = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
        const JData jd      = load2xnJData(atom, cj);
        const RowMasks diag = load2xnDiagMasks(atom->masks_2xn_fn, ci, cj);
        const PairForces p =
            pairForces(CLUSTER_M, nrlj, energy, &ir, &jd, &diag, KERNEL_EXCL, pc, atom);

        accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
        if (energy) accumulateEnergies(CLUSTER_M, &e, &p);
    }

    for (int k = numneighs_masked; k < unmasked_inner_end; k++) {
        prefetchAhead(atom,
            neighbor,
            ci,
            nbM,
            k - numneighs_masked + numneighs_inner_masked,
            numneighs,
            numneighs_masked,
            numneighs_inner_masked,
            NULL);
        const int cj   = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
        const JData jd = load2xnJData(atom, cj);
        const PairForces p =
            pairForces(CLUSTER_M, nrlj, energy, &ir, &jd, NULL, 0, pc, atom);

        accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
        if (energy) accumulateEnergies(CLUSTER_M, &e, &p);
    }

    incrementIForces(ci_f, &fx, &fy, &fz);
    // A full list has every pair twice
    if (energy) reduceEnergies(CLUSTER_M, &e, 0.5, ec, elj);

    addStat(stats->calculated_forces, 1);
    addStat(stats->num_neighs, numneighs);
    addStat(stats->force_iters, (long long int)((double)numneighs));
}

double VARIANT_NAME(computeForceLJ2xnFullNeigh)(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_2xn begin\n");
    const PairConstants pc = pairConstants(param);
    const int energy       = param->compute_energy;
    double ec = 0.0, elj = 0.0;
#if LJ_COMB_RULE == LJ_COMB_FULL
    setupTypePairC6C12(atom);
#endif

    zeroForces(atom, 0);

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime) reduction(+ : ec, elj)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            const int half_lj = param->half_lj && atom->iclusters[ci].half_lj;
            if (!energy) {
                if (!half_lj) {
                    iClusterFull2xn(atom,
                        neighbor,
                        stats,
                        &pc,
                        ci,
                        CLUSTER_M,
                        0,
                        &ec,
                        &elj);
                } else {
                    iClusterFull2xn(atom,
                        neighbor,
                        stats,
                        &pc,
                        ci,
                        CLUSTER_M / 2,
                        0,
                        &ec,
                        &elj);
                }
            } else {
                if (!half_lj) {
                    iClusterFull2xn(atom,
                        neighbor,
                        stats,
                        &pc,
                        ci,
                        CLUSTER_M,
                        1,
                        &ec,
                        &elj);
                } else {
                    iClusterFull2xn(atom,
                        neighbor,
                        stats,
                        &pc,
                        ci,
                        CLUSTER_M / 2,
                        1,
                        &ec,
                        &elj);
                }
            }
        }

        LIKWID_MARKER_STOP("force");
    }

    if (energy) {
        storeEnergies(param,
            elj,
            ec,
            localChargeSquaredSum(atom->Nlocal, atom->type, atom->q_per_type));
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_2xn end\n");
    return E - S;
}

#endif
