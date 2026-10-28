/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* SIMD 2xNN cluster-pair kernels: a SIMD register holds 2 i-atoms in its two halves
 * against the N = SIMD width / 2 atoms of a j-cluster in both halves, so an i-cluster of
 * 4 atoms is 2 rows. Half and full neighbor lists. The pair math is shared with the other
 * SIMD kernels, see pair_rows.h. */
#include <stdio.h>

#include <atom.h>
#include <force.h>
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

#ifndef USE_REFERENCE_KERNEL

void computeForceGhostShell(Parameter*, Atom*, Neighbor*);

/* A row holds two i-atoms, one in each half of the register, so the 4 atoms of the
 * i-cluster are CLUSTER_M / 2 rows. Row r holds the i-atoms 2r and 2r + 1, one in each
 * half of the register */
static inline IRows load2xnnIRows(const Atom* atom, const int ci, const PairConstants* pc)
{
    IRows ir;
    const MD_FLOAT* ci_x = &atom->cl_x[CI_VECTOR3_BASE_INDEX(ci)];
    const int sca        = CI_SCALAR_BASE_INDEX(ci);

    for (int r = 0; r < (CLUSTER_M / 2); r++) {
        ir.x.v[r] = simd_real_load_h_dual(&ci_x[CL_X_INDEX_3D(2 * r)]);
        ir.y.v[r] = simd_real_load_h_dual(&ci_x[CL_Y_INDEX_3D(2 * r)]);
        ir.z.v[r] = simd_real_load_h_dual(&ci_x[CL_Z_INDEX_3D(2 * r)]);
#if COULOMB_RULE != COULOMB_NONE
        ir.q.v[r] = simd_real_mul(pc->cs.epsfac,
            simd_real_load_h_dual(&atom->cl_q[sca + 2 * r]));
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
        ir.a6.v[r]  = simd_real_load_h_dual(&atom->cl_sqrt6c6[sca + 2 * r]);
        ir.b12.v[r] = simd_real_load_h_dual(&atom->cl_sqrt12c12[sca + 2 * r]);
#elif LJ_COMB_RULE == LJ_COMB_FULL
        ir.type.v[r] = simd_i32_load_h_dual_scaled(&atom->cl_t[sca + 2 * r],
            atom->ntypes);
#endif
#if USE_EXCLUSIONS
        ir.mol.v[r] = simd_i32_load_h_dual_scaled(&atom->cl_mol[sca + 2 * r], 1);
#endif
    }

    (void)sca;
    (void)pc;
    return ir;
}

/* The j-cluster is duplicated into both halves of the register */
static inline JData load2xnnJData(const Atom* atom, const int cj)
{
    JData jd;
    MD_FLOAT* cj_x = &atom->cl_x[CJ_VECTOR3_BASE_INDEX(cj)];
    const int sca  = CJ_SCALAR_BASE_INDEX(cj);

    jd.x = simd_real_load_h_duplicate(&cj_x[CL_X_INDEX_3D(0)]);
    jd.y = simd_real_load_h_duplicate(&cj_x[CL_Y_INDEX_3D(0)]);
    jd.z = simd_real_load_h_duplicate(&cj_x[CL_Z_INDEX_3D(0)]);
#if COULOMB_RULE != COULOMB_NONE
    jd.q = simd_real_load_h_duplicate(&atom->cl_q[sca]);
#endif
#if LJ_COMB_RULE == LJ_COMB_GEOM
    jd.a6  = simd_real_load_h_duplicate(&atom->cl_sqrt6c6[sca]);
    jd.b12 = simd_real_load_h_duplicate(&atom->cl_sqrt12c12[sca]);
#elif LJ_COMB_RULE == LJ_COMB_FULL
    jd.type = simd_i32_load_h_duplicate(&atom->cl_t[sca]);
#endif
#if USE_EXCLUSIONS
    jd.mol = simd_i32_load_h_duplicate(&atom->cl_mol[sca]);
#endif

    (void)sca;
    return jd;
}

/* Rows of one entry of the diagonal mask table (the entries are selected by which
 * of the i-cluster halves is the j-cluster). */
#define DIAG_MASK_ROWS (CLUSTER_M / 2)

/* Masks of the pairs to skip in a pair of clusters that is, or overlaps with, the
 * cluster itself, see the 4xN kernels */
static inline RowMasks load2xnnDiagMasks(
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
    for (int r = 0; r < (CLUSTER_M / 2); r++) {
        diag.v[r] = simd_mask_from_u32(table[idx + r]);
    }
    return diag;
}

/* Reaction force on the atoms of the j-cluster */
static inline void decrementJForces(MD_FLOAT* cj_f, const PairForces* p)
{
    simd_real_h_decr3(cj_f,
        rowsSum((CLUSTER_M / 2), p->tx),
        rowsSum((CLUSTER_M / 2), p->ty),
        rowsSum((CLUSTER_M / 2), p->tz));
}

static inline void incrementIForces(
    MD_FLOAT* ci_f, const Rows* fx, const Rows* fy, const Rows* fz)
{
    simd_real_h_dual_incr_reduced_sum(&ci_f[CL_X_INDEX_3D(0)], fx->v[0], fx->v[1]);
    simd_real_h_dual_incr_reduced_sum(&ci_f[CL_Y_INDEX_3D(0)], fy->v[0], fy->v[1]);
    simd_real_h_dual_incr_reduced_sum(&ci_f[CL_Z_INDEX_3D(0)], fz->v[0], fz->v[1]);
}

double computeForceLJ2xnnHalfNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_2xnn begin\n");
    const int nbM             = atom->Nclusters_local;
    const PairConstants pc    = pairConstants(param);
    const int nlocal_clusters = CJ0_FROM_CI(atom->Nclusters_local);
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

#pragma omp for schedule(runtime)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            const int numneighs              = neighbor->numneigh_inner[ci];
            const int numneighs_masked       = neighbor->numneigh_masked[ci];
            const int numneighs_inner_masked = neighbor->numneigh_inner_masked[ci];
            const int unmasked_inner_end     = numneighs_masked +
                                           (numneighs - numneighs_inner_masked);
            MD_FLOAT* ci_f = &f[CI_VECTOR3_BASE_INDEX(ci)];

            const IRows ir = load2xnnIRows(atom, ci, &pc);
            Rows fx        = rowsZero();
            Rows fy        = rowsZero();
            Rows fz        = rowsZero();

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
                const JData jd      = load2xnnJData(atom, cj);
                const RowMasks diag = load2xnnDiagMasks(atom->masks_2xnn_hn, ci, cj);
                const PairForces p =
                    pairForces((CLUSTER_M / 2), &ir, &jd, &diag, 1, &pc, atom);

                accumulateForces((CLUSTER_M / 2), &fx, &fy, &fz, &p);
                if (cj < nlocal_clusters || param->method) {
                    decrementJForces(&f[CJ_VECTOR3_BASE_INDEX(cj)], &p);
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
                const JData jd = load2xnnJData(atom, cj);
                const PairForces p =
                    pairForces((CLUSTER_M / 2), &ir, &jd, NULL, 0, &pc, atom);

                accumulateForces((CLUSTER_M / 2), &fx, &fy, &fz, &p);
                if (cj < nlocal_clusters || param->method) {
                    decrementJForces(&f[CJ_VECTOR3_BASE_INDEX(cj)], &p);
                }
            }

            incrementIForces(ci_f, &fx, &fy, &fz);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters,
                (long long int)((double)numneighs * CLUSTER_M / CLUSTER_N));
        }
        threadForcesReduce(atom, nf);
        LIKWID_MARKER_STOP("force");
    }

    // Single threaded: it writes cl_f of the shell clusters without a thread buffer
    if (param->method == eightShell) {
        computeForceGhostShell(param, atom, neighbor);
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_2xnn end\n");
    return E - S;
}

double computeForceLJ2xnnFullNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_2xnn begin\n");
    const int nbM          = atom->Nclusters_local;
    const PairConstants pc = pairConstants(param);
#if LJ_COMB_RULE == LJ_COMB_FULL
    setupTypePairC6C12(atom);
#endif

    zeroForces(atom, 0);

    double S = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for schedule(runtime)
        for (int ci = 0; ci < atom->Nclusters_local; ci++) {
            const int numneighs              = neighbor->numneigh_inner[ci];
            const int numneighs_masked       = neighbor->numneigh_masked[ci];
            const int numneighs_inner_masked = neighbor->numneigh_inner_masked[ci];
            const int unmasked_inner_end     = numneighs_masked +
                                           (numneighs - numneighs_inner_masked);
            MD_FLOAT* ci_f = &atom->cl_f[CI_VECTOR3_BASE_INDEX(ci)];

            const IRows ir = load2xnnIRows(atom, ci, &pc);
            Rows fx        = rowsZero();
            Rows fy        = rowsZero();
            Rows fz        = rowsZero();

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
                const JData jd      = load2xnnJData(atom, cj);
                const RowMasks diag = load2xnnDiagMasks(atom->masks_2xnn_fn, ci, cj);
                const PairForces p =
                    pairForces((CLUSTER_M / 2), &ir, &jd, &diag, 1, &pc, atom);

                accumulateForces((CLUSTER_M / 2), &fx, &fy, &fz, &p);
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
                const JData jd = load2xnnJData(atom, cj);
                const PairForces p =
                    pairForces((CLUSTER_M / 2), &ir, &jd, NULL, 0, &pc, atom);

                accumulateForces((CLUSTER_M / 2), &fx, &fy, &fz, &p);
            }

            incrementIForces(ci_f, &fx, &fy, &fz);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters, (long long int)((double)numneighs));
        }

        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_2xnn end\n");
    return E - S;
}

#endif
