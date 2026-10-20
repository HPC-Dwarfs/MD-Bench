/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* SIMD 4xN cluster-pair kernels: 4 i-atoms against N = SIMD width j-atoms, with half and
 * full neighbor lists. The pair math is shared with the other SIMD kernels, see
 * pair_rows.h. */
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
#include <timing.h>
#include <util.h>

#ifndef USE_REFERENCE_KERNEL

void computeForceGhostShell(Parameter*, Atom*, Neighbor*);

static inline IRows load4xnIRows(const Atom* atom, const int ci, const PairConstants* pc)
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
#if USE_EXCLUSIONS
        ir.mol.v[r] = simd_i32_broadcast(atom->cl_mol[sca + r]);
#endif
    }

    (void)sca;
    (void)pc;
    return ir;
}

static inline JData load4xnJData(const Atom* atom, const int cj)
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
#if USE_EXCLUSIONS
    jd.mol = simd_i32_load(&atom->cl_mol[sca]);
#endif

    (void)sca;
    return jd;
}

/* Rows of one entry of the diagonal mask table (the entries are selected by which
 * of the i-cluster halves is the j-cluster). */
#define DIAG_MASK_ROWS CLUSTER_M

/* Masks of the pairs to skip in a pair of clusters that is, or overlaps with, the
 * cluster itself: a pair of an atom with itself or, in a half list, one that is counted
 * from the other side as well. The table depends on the shape of the clusters. */
static inline RowMasks load4xnDiagMasks(
    const Atom* atom, const unsigned int* table, const int ci, const int cj)
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
    (void)atom;
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
    simd_real_incr_reduced_sum(&ci_f[CL_X_INDEX_3D(0)],
        fx->v[0],
        fx->v[1],
        fx->v[2],
        fx->v[3]);
    simd_real_incr_reduced_sum(&ci_f[CL_Y_INDEX_3D(0)],
        fy->v[0],
        fy->v[1],
        fy->v[2],
        fy->v[3]);
    simd_real_incr_reduced_sum(&ci_f[CL_Z_INDEX_3D(0)],
        fz->v[0],
        fz->v[1],
        fz->v[2],
        fz->v[3]);
}

double computeForceLJ4xnHalfNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_4xn begin\n");
    const int nbM             = atom->Nclusters_local;
    const PairConstants pc    = pairConstants(param);
    const int nlocal_clusters = CJ1_FROM_CI(atom->Nclusters_local);
#if LJ_COMB_RULE == LJ_COMB_FULL
    setupTypePairC6C12(atom);
#endif

    zeroForces(atom, 1);

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

            const IRows ir = load4xnIRows(atom, ci, &pc);
            Rows fx        = rowsZero();
            Rows fy        = rowsZero();
            Rows fz        = rowsZero();

            for (int k = 0; k < numneighs_inner_masked; k++) {
                if (k + KERNEL_PREFETCH_DIST < numneighs) {
                    prefetchJ(atom,
                        neighs(neighbor->neighbors,
                            ci,
                            k + KERNEL_PREFETCH_DIST,
                            nbM,
                            neighbor),
                        1);
                }
                const int cj        = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
                const JData jd      = load4xnJData(atom, cj);
                const RowMasks diag = load4xnDiagMasks(atom, atom->masks_4xn_hn, ci, cj);
                const PairForces p = pairForces(CLUSTER_M, &ir, &jd, &diag, 1, &pc, atom);

                accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
                if (cj < nlocal_clusters || param->method) {
                    decrementJForces(&atom->cl_f[CJ_VECTOR3_BASE_INDEX(cj)], &p);
                }
            }

            for (int k = numneighs_masked; k < unmasked_inner_end; k++) {
                if (k + KERNEL_PREFETCH_DIST < numneighs) {
                    prefetchJ(atom,
                        neighs(neighbor->neighbors,
                            ci,
                            k + KERNEL_PREFETCH_DIST,
                            nbM,
                            neighbor),
                        1);
                }
                const int cj       = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
                const JData jd     = load4xnJData(atom, cj);
                const PairForces p = pairForces(CLUSTER_M, &ir, &jd, NULL, 0, &pc, atom);

                accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
                if (cj < nlocal_clusters || param->method) {
                    decrementJForces(&atom->cl_f[CJ_VECTOR3_BASE_INDEX(cj)], &p);
                }
            }

            incrementIForces(ci_f, &fx, &fy, &fz);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters,
                (long long int)((double)numneighs * CLUSTER_M / CLUSTER_N));
        }
        if (param->method == eightShell) computeForceGhostShell(param, atom, neighbor);
        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_4xn end\n");
    return E - S;
}

double computeForceLJ4xnFullNeigh(
    Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJ_4xn begin\n");
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

            const IRows ir = load4xnIRows(atom, ci, &pc);
            Rows fx        = rowsZero();
            Rows fy        = rowsZero();
            Rows fz        = rowsZero();

            for (int k = 0; k < numneighs_inner_masked; k++) {
                if (k + KERNEL_PREFETCH_DIST < numneighs) {
                    prefetchJ(atom,
                        neighs(neighbor->neighbors,
                            ci,
                            k + KERNEL_PREFETCH_DIST,
                            nbM,
                            neighbor),
                        0);
                }
                const int cj        = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
                const JData jd      = load4xnJData(atom, cj);
                const RowMasks diag = load4xnDiagMasks(atom, atom->masks_4xn_fn, ci, cj);
                const PairForces p = pairForces(CLUSTER_M, &ir, &jd, &diag, 1, &pc, atom);

                accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
            }

            for (int k = numneighs_masked; k < unmasked_inner_end; k++) {
                if (k + KERNEL_PREFETCH_DIST < numneighs) {
                    prefetchJ(atom,
                        neighs(neighbor->neighbors,
                            ci,
                            k + KERNEL_PREFETCH_DIST,
                            nbM,
                            neighbor),
                        0);
                }
                const int cj       = neighs(neighbor->neighbors, ci, k, nbM, neighbor);
                const JData jd     = load4xnJData(atom, cj);
                const PairForces p = pairForces(CLUSTER_M, &ir, &jd, NULL, 0, &pc, atom);

                accumulateForces(CLUSTER_M, &fx, &fy, &fz, &p);
            }

            incrementIForces(ci_f, &fx, &fy, &fz);

            addStat(stats->calculated_forces, 1);
            addStat(stats->num_neighs, numneighs);
            addStat(stats->force_iters, (long long int)((double)numneighs));
        }

        LIKWID_MARKER_STOP("force");
    }

    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJ_4xn end\n");
    return E - S;
}

#endif
