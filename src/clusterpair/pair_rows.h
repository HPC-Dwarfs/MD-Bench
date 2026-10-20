/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* Pair interaction of a cluster pair for the SIMD cluster-pair kernels (4xN, 2xNN, 2xN).
 *
 * A kernel keeps the i-cluster in "rows": each row is a SIMD register pair of i-atom
 * data against the j-atom data of a vector (one atom broadcast in the 4xN and 2xN
 * kernels, two atoms in the two halves of the register in the 2xNN kernel). The kernels
 * only differ in how rows and j-data are loaded and how forces are reduced, the math is
 * the same and lives here:
 *
 *  - one 1/sqrt per pair gives 1/r and 1/r^2, and the cut-off (and, in the masked loop,
 *    the diagonal) mask is applied to 1/r only, which zeroes every term of the pair;
 *  - LJ is evaluated from pre-multiplied 6*C6 and 12*C12 as
 *    (12*C12 * r^-12 - 6*C6 * r^-6) / r^2, with the C6/C12 of a pair obtained as the
 *    product of per-atom values for the geometric combination rule;
 *  - electrostatics (reaction field or real-space Ewald) and molecular exclusions as in
 *    GROMACS' NBNxM kernels, see coulomb.h;
 *  - everything is computed phase by phase over the rows (all rows finish a phase before
 *    the next one starts), which hands the instruction scheduler independent work to fill
 *    the latency of every phase. The kernels written this way are 15-25% faster than the
 *    row after row formulation.
 *
 * The functions are always inlined, the number of rows and the exclusion flag are
 * compile-time constants at every call site.
 */
#ifndef __PAIR_ROWS_H_
#define __PAIR_ROWS_H_

#include <allocate.h>
#include <atom.h>
#include <coulomb.h>
#include <force.h>
#include <kernel_math.h>
#include <parameter.h>
#include <simd.h>

#ifndef USE_REFERENCE_KERNEL

#define ROWS_MAX 4

#if defined(__GNUC__) || defined(__clang__)
#define MD_ALWAYS_INLINE static inline __attribute__((always_inline))
#else
#define MD_ALWAYS_INLINE static inline
#endif

typedef struct {
    MD_SIMD_FLOAT v[ROWS_MAX];
} Rows;

typedef struct {
    MD_SIMD_MASK v[ROWS_MAX];
} RowMasks;

typedef struct {
    MD_SIMD_INT v[ROWS_MAX];
} RowInts;

/* i-cluster data of every row */
typedef struct {
    Rows x, y, z;
    Rows q;       // charge times the Coulomb prefactor
    Rows a6, b12; // sqrt(6 C6) and sqrt(12 C12), geometric combination rule
    RowInts type; // type * ntypes, full combination rule
    RowInts mol;  // molecule index, exclusions
} IRows;

/* j-cluster data of the vector */
typedef struct {
    MD_SIMD_FLOAT x, y, z;
    MD_SIMD_FLOAT q;
    MD_SIMD_FLOAT a6, b12;
    MD_SIMD_INT type;
    MD_SIMD_INT mol;
} JData;

/* Values that stay the same for the whole kernel call */
typedef struct {
    MD_SIMD_FLOAT cutsq;
    MD_SIMD_FLOAT c6, c12; // 6 C6 and 12 C12 of the single atom type rule
    CoulombSimd cs;
} PairConstants;

/* Forces of the pairs, already masked, as force * distance */
typedef struct {
    Rows tx, ty, tz;
} PairForces;

static inline PairConstants pairConstants(const Parameter* param)
{
    PairConstants pc;
    const MD_FLOAT cutsq = param->cutforce * param->cutforce;
    pc.cutsq             = simd_real_broadcast(cutsq);
    // 6 C6 = 24 eps sigma^6, 12 C12 = 48 eps sigma^12
    pc.c6  = simd_real_broadcast((MD_FLOAT)24.0 * param->epsilon * param->sigma6);
    pc.c12 = simd_real_broadcast(
        (MD_FLOAT)48.0 * param->epsilon * param->sigma6 * param->sigma6);
    pc.cs = coulombSimdInit((Parameter*)param);
    return pc;
}

MD_ALWAYS_INLINE Rows rowsZero(void)
{
    Rows r;
    for (int i = 0; i < ROWS_MAX; i++) {
        r.v[i] = simd_real_zero();
    }
    return r;
}

MD_ALWAYS_INLINE Rows rowsAdd(const int nr, Rows a, Rows b)
{
    for (int i = 0; i < nr; i++) {
        a.v[i] = simd_real_add(a.v[i], b.v[i]);
    }
    return a;
}

/* Sum over the rows */
MD_ALWAYS_INLINE MD_SIMD_FLOAT rowsSum(const int nr, Rows a)
{
    MD_SIMD_FLOAT sum = a.v[0];
    for (int i = 1; i < nr; i++) {
        sum = simd_real_add(sum, a.v[i]);
    }
    return sum;
}

/* Adds the forces of a cluster pair to the i-forces of the rows */
MD_ALWAYS_INLINE void accumulateForces(
    const int nr, Rows* fx, Rows* fy, Rows* fz, const PairForces* p)
{
    for (int r = 0; r < nr; r++) {
        fx->v[r] = simd_real_add(fx->v[r], p->tx.v[r]);
        fy->v[r] = simd_real_add(fy->v[r], p->ty.v[r]);
        fz->v[r] = simd_real_add(fz->v[r], p->tz.v[r]);
    }
}

/* force (positive for repulsion) of the pairs of nr rows of the i-cluster with the
 * j-vector; diag is the mask of the pairs to skip because they are the pair of an atom
 * with itself or already counted from the other side (NULL for none) and excl selects
 * the molecular exclusions */
MD_ALWAYS_INLINE PairForces pairForces(const int nr,
    const IRows* ir,
    const JData* jd,
    const RowMasks* diag,
    const int excl,
    const PairConstants* pc,
    const Atom* atom)
{
    Rows dx, dy, dz, rsq, cutsq, c6, c12, rinv, rinvsq, six, frlj, force;
    RowMasks wm;
    PairForces out;

    for (int r = 0; r < nr; r++) {
        dx.v[r] = simd_real_sub(ir->x.v[r], jd->x);
        dy.v[r] = simd_real_sub(ir->y.v[r], jd->y);
        dz.v[r] = simd_real_sub(ir->z.v[r], jd->z);
    }
    for (int r = 0; r < nr; r++) {
        rsq.v[r] = simd_real_fma(dx.v[r],
            dx.v[r],
            simd_real_fma(dy.v[r], dy.v[r], simd_real_mul(dz.v[r], dz.v[r])));
    }

#if LJ_COMB_RULE == LJ_COMB_FULL
    for (int r = 0; r < nr; r++) {
        const MD_SIMD_INT tvec = simd_i32_add(ir->type.v[r], jd->type);
        cutsq.v[r] = simd_real_gather(tvec, atom->cutforcesq, sizeof(MD_FLOAT));
        c6.v[r]    = simd_real_gather(tvec, atom->ljc6x6, sizeof(MD_FLOAT));
        c12.v[r]   = simd_real_gather(tvec, atom->ljc12x12, sizeof(MD_FLOAT));
    }
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    for (int r = 0; r < nr; r++) {
        cutsq.v[r] = pc->cutsq;
        c6.v[r]    = simd_real_mul(ir->a6.v[r], jd->a6);
        c12.v[r]   = simd_real_mul(ir->b12.v[r], jd->b12);
    }
#else
    for (int r = 0; r < nr; r++) {
        cutsq.v[r] = pc->cutsq;
        c6.v[r]    = pc->c6;
        c12.v[r]   = pc->c12;
    }
#endif

    for (int r = 0; r < nr; r++) {
        wm.v[r] = simd_mask_cond_lt(rsq.v[r], cutsq.v[r]);
        if (diag != NULL) {
            wm.v[r] = simd_mask_and(diag->v[r], wm.v[r]);
        }
    }

#if COULOMB_RULE == COULOMB_NONE
    // Without electrostatics only 1/r^2 is needed, which is cheaper than via 1/r
    for (int r = 0; r < nr; r++) {
        rinvsq.v[r] = simd_real_select_by_mask(SIMD_RECIPROCAL(rsq.v[r]), wm.v[r]);
    }
#else
    for (int r = 0; r < nr; r++) {
        rinv.v[r] = simd_real_select_by_mask(SIMD_RSQRT(rsq.v[r]), wm.v[r]);
    }
    for (int r = 0; r < nr; r++) {
        rinvsq.v[r] = simd_real_mul(rinv.v[r], rinv.v[r]);
    }
#endif

    RowMasks interact;
#if USE_EXCLUSIONS
    if (excl) {
        for (int r = 0; r < nr; r++) {
            interact.v[r] = simd_mask_not(simd_mask_i32_cond_eq(ir->mol.v[r], jd->mol));
        }
    }
#endif
    for (int r = 0; r < nr; r++) {
        six.v[r] = simd_real_mul(rinvsq.v[r], simd_real_mul(rinvsq.v[r], rinvsq.v[r]));
#if USE_EXCLUSIONS
        if (excl) {
            six.v[r] = simd_real_select_by_mask(six.v[r], interact.v[r]);
        }
#endif
    }

    for (int r = 0; r < nr; r++) {
        frlj.v[r] = simd_real_mul(six.v[r],
            simd_real_sub(simd_real_mul(c12.v[r], six.v[r]), c6.v[r]));
    }

#if COULOMB_RULE != COULOMB_NONE
    Rows rinv_coul;
    for (int r = 0; r < nr; r++) {
#if USE_EXCLUSIONS
        rinv_coul.v[r] = excl ? simd_real_select_by_mask(rinv.v[r], interact.v[r])
                              : rinv.v[r];
#else
        rinv_coul.v[r] = rinv.v[r];
#endif
    }
    Rows fcoul;
#if COULOMB_RULE == COULOMB_RF
    for (int r = 0; r < nr; r++) {
        const MD_SIMD_FLOAT qq = simd_real_select_by_mask(
            simd_real_mul(ir->q.v[r], jd->q),
            wm.v[r]);
        fcoul.v[r] = simd_real_mul(qq,
            simd_real_sub(simd_real_mul(rinv_coul.v[r], rinvsq.v[r]), pc->cs.k_rf2));
    }
#else
    Rows brsq, ewc;
    for (int r = 0; r < nr; r++) {
        brsq.v[r] = simd_real_mul(pc->cs.beta2,
            simd_real_select_by_mask(rsq.v[r], wm.v[r]));
    }
    for (int r = 0; r < nr; r++) {
        ewc.v[r] = simd_real_mul(pc->cs.beta, pmeForceCorrectionSimd(brsq.v[r]));
    }
    for (int r = 0; r < nr; r++) {
        fcoul.v[r] = simd_real_mul(simd_real_mul(ir->q.v[r], jd->q),
            simd_real_mul(simd_real_fma(ewc.v[r], brsq.v[r], rinv_coul.v[r]),
                rinvsq.v[r]));
    }
#endif
    for (int r = 0; r < nr; r++) {
        force.v[r] = simd_real_fma(frlj.v[r], rinvsq.v[r], fcoul.v[r]);
    }
#else
    for (int r = 0; r < nr; r++) {
        force.v[r] = simd_real_mul(frlj.v[r], rinvsq.v[r]);
    }
#endif

    // The padding atoms of a cluster are far away but finite (CLUSTER_PAD_COORD), the
    // force of a masked pair is exactly zero
    for (int r = 0; r < nr; r++) {
        out.tx.v[r] = simd_real_mul(dx.v[r], force.v[r]);
        out.ty.v[r] = simd_real_mul(dy.v[r], force.v[r]);
        out.tz.v[r] = simd_real_mul(dz.v[r], force.v[r]);
    }

    (void)atom;
    (void)excl;
    (void)pc;
    return out;
}

/* Clears the forces of the local clusters and, for kernels that write reaction forces
 * into them, of the ghost clusters */
static inline void zeroForces(Atom* atom, const int ghosts)
{
    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        MD_FLOAT* ci_f = &atom->cl_f[CI_VECTOR3_BASE_INDEX(ci)];
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            ci_f[CL_X_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Y_INDEX_3D(cii)] = 0.0;
            ci_f[CL_Z_INDEX_3D(cii)] = 0.0;
        }
    }

    if (ghosts) {
        for (int cg = atom->ncj; cg < atom->ncj + atom->Nclusters_ghost; cg++) {
            MD_FLOAT* cj_f = &atom->cl_f[CJ_VECTOR3_BASE_INDEX(cg)];
            for (int cjj = 0; cjj < atom->jclusters[cg].natoms; cjj++) {
                cj_f[CL_X_INDEX_3D(cjj)] = 0.0;
                cj_f[CL_Y_INDEX_3D(cjj)] = 0.0;
                cj_f[CL_Z_INDEX_3D(cjj)] = 0.0;
            }
        }
    }
}

/* j-cluster data prefetch, a few list entries ahead: the list is an indirection into
 * data that does not fit the L1 cache. Define KERNEL_NO_PREFETCH to disable. */
#ifndef KERNEL_PREFETCH_DIST
#define KERNEL_PREFETCH_DIST 3
#endif

MD_ALWAYS_INLINE void prefetchJ(const Atom* atom, const int cj, const int write_f)
{
#ifndef KERNEL_NO_PREFETCH
    const char* px = (const char*)&atom->cl_x[CJ_VECTOR3_BASE_INDEX(cj)];
    __builtin_prefetch(px, 0, 3);
    __builtin_prefetch(px + 64, 0, 3);
    if (write_f) {
        const char* pf = (const char*)&atom->cl_f[CJ_VECTOR3_BASE_INDEX(cj)];
        __builtin_prefetch(pf, 1, 3);
        __builtin_prefetch(pf + 64, 1, 3);
    }
#if LJ_COMB_RULE == LJ_COMB_GEOM
    __builtin_prefetch(&atom->cl_sqrt6c6[CJ_SCALAR_BASE_INDEX(cj)], 0, 3);
    __builtin_prefetch(&atom->cl_sqrt12c12[CJ_SCALAR_BASE_INDEX(cj)], 0, 3);
#endif
#if COULOMB_RULE != COULOMB_NONE
    __builtin_prefetch(&atom->cl_q[CJ_SCALAR_BASE_INDEX(cj)], 0, 3);
#endif
#else
    (void)atom;
    (void)cj;
    (void)write_f;
#endif
}

/* Type-pair tables of 6*C6 and 12*C12 of the full combination rule, built on first use */
#if LJ_COMB_RULE == LJ_COMB_FULL
static inline void setupTypePairC6C12(Atom* atom)
{
    if (atom->ljc6x6 != NULL) {
        return;
    }
    const int n2   = atom->ntypes * atom->ntypes;
    atom->ljc6x6   = (MD_FLOAT*)allocate(ALIGNMENT, n2 * sizeof(MD_FLOAT));
    atom->ljc12x12 = (MD_FLOAT*)allocate(ALIGNMENT, n2 * sizeof(MD_FLOAT));
    for (int i = 0; i < n2; i++) {
        const MD_FLOAT s6 = atom->sigma6[i];
        atom->ljc6x6[i]   = (MD_FLOAT)24.0 * atom->epsilon[i] * s6;
        atom->ljc12x12[i] = (MD_FLOAT)48.0 * atom->epsilon[i] * s6 * s6;
    }
}
#endif

#endif /* !USE_REFERENCE_KERNEL */

#endif
