/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <allocate.h>
#include <atom.h>
#include <force.h>
#include <neighbor.h>
#include <pbc.h>
#include <util.h>

static int NmaxGhost;

#ifdef CUDA_TARGET
UpdatePbcFunction updatePbc      = updatePbcCUDA;
UpdatePbcFunction updateAtomsPbc = updateAtomsPbcCPU;
// UpdatePbcFunction updateAtomsPbc = updateAtomsPbcCUDA;
#else
UpdatePbcFunction updatePbc      = updatePbcCPU;
UpdatePbcFunction updateAtomsPbc = updateAtomsPbcCPU;
#endif

static void growPbc(Atom*);

/* exported subroutines */
void initPbc(Atom* atom) {
    NmaxGhost        = 0;
    atom->border_map = NULL;
    atom->PBCx       = NULL;
    atom->PBCy       = NULL;
    atom->PBCz       = NULL;
}

/* update coordinates of ghost atoms */
/* uses mapping created in setupPbc */
void updatePbcCPU(Atom* atom, Parameter* param, bool firstUpdate) {
    DEBUG_MESSAGE("updatePbc start\n");
    int ncj       = get_ncj_from_nci(atom->Nclusters_local);
    MD_FLOAT xprd = param->xprd;
    MD_FLOAT yprd = param->yprd;
    MD_FLOAT zprd = param->zprd;

    for (int cg = 0; cg < atom->Nclusters_ghost; cg++) {
        const int cj    = ncj + cg;
        int cjScaBase   = CJ_SCALAR_BASE_INDEX(cj);
        int cjVecBase   = CJ_VECTOR_BASE_INDEX(cj);
        int bmapScaBase = CJ_SCALAR_BASE_INDEX(atom->border_map[cg]);
        int bmapVecBase = CJ_VECTOR_BASE_INDEX(atom->border_map[cg]);
        int* cjT        = &atom->cl_t[cjScaBase];
        MD_FLOAT* cjX   = &atom->cl_x[cjVecBase];
        MD_FLOAT* bmapX = &atom->cl_x[bmapVecBase];
        int* bmapT      = &atom->cl_t[bmapScaBase];
        MD_FLOAT bbminx = INFINITY, bbmaxx = -INFINITY;
        MD_FLOAT bbminy = INFINITY, bbmaxy = -INFINITY;
        MD_FLOAT bbminz = INFINITY, bbmaxz = -INFINITY;

        for (int cjj = 0; cjj < atom->jclusters[cj].natoms; cjj++) {
            MD_FLOAT xtmp = bmapX[CL_X_INDEX(cjj)] + atom->PBCx[cg] * xprd;
            MD_FLOAT ytmp = bmapX[CL_Y_INDEX(cjj)] + atom->PBCy[cg] * yprd;
            MD_FLOAT ztmp = bmapX[CL_Z_INDEX(cjj)] + atom->PBCz[cg] * zprd;

            cjX[CL_X_INDEX(cjj)] = xtmp;
            cjX[CL_Y_INDEX(cjj)] = ytmp;
            cjX[CL_Z_INDEX(cjj)] = ztmp;
            cjT[cjj]               = bmapT[cjj];

            if (firstUpdate) {
                // TODO: To create the bounding boxes faster, we can use SIMD operations
                if (bbminx > xtmp) {
                    bbminx = xtmp;
                }
                if (bbmaxx < xtmp) {
                    bbmaxx = xtmp;
                }
                if (bbminy > ytmp) {
                    bbminy = ytmp;
                }
                if (bbmaxy < ytmp) {
                    bbmaxy = ytmp;
                }
                if (bbminz > ztmp) {
                    bbminz = ztmp;
                }
                if (bbmaxz < ztmp) {
                    bbmaxz = ztmp;
                }
            }
        }

        if (firstUpdate) {
            for (int cjj = atom->jclusters[cj].natoms; cjj < CLUSTER_N; cjj++) {
                cjX[CL_X_INDEX(cjj)] = INFINITY;
                cjX[CL_Y_INDEX(cjj)] = INFINITY;
                cjX[CL_Z_INDEX(cjj)] = INFINITY;
                cjT[cjj]               = 0;
            }

            atom->jclusters[cj].bbminx = bbminx;
            atom->jclusters[cj].bbmaxx = bbmaxx;
            atom->jclusters[cj].bbminy = bbminy;
            atom->jclusters[cj].bbmaxy = bbmaxy;
            atom->jclusters[cj].bbminz = bbminz;
            atom->jclusters[cj].bbmaxz = bbmaxz;
        }
    }

    DEBUG_MESSAGE("updatePbc end\n");
}

/* relocate atoms that have left domain according
 * to periodic boundary conditions */
void updateAtomsPbcCPU(Atom* atom, Parameter* param, bool dummy) {
    MD_FLOAT xprd = param->xprd;
    MD_FLOAT yprd = param->yprd;
    MD_FLOAT zprd = param->zprd;

    for (int i = 0; i < atom->Nlocal; i++) {
        if (atom_x(i) < 0.0) {
            atom_x(i) += xprd;
        } else if (atom_x(i) >= xprd) {
            atom_x(i) -= xprd;
        }

        if (atom_y(i) < 0.0) {
            atom_y(i) += yprd;
        } else if (atom_y(i) >= yprd) {
            atom_y(i) -= yprd;
        }

        if (atom_z(i) < 0.0) {
            atom_z(i) += zprd;
        } else if (atom_z(i) >= zprd) {
            atom_z(i) -= zprd;
        }
    }
}

/* setup periodic boundary conditions by
 * defining ghost atoms around domain
 * only creates mapping and coordinate corrections
 * that are then enforced in updatePbc */
#define ADDGHOST(dx, dy, dz)                                                             \
    Nghost++;                                                                            \
    const int cg               = ncj + Nghost;                                           \
    const int cj_natoms        = atom->jclusters[cj].natoms;                             \
    atom->border_map[Nghost]   = cj;                                                     \
    atom->PBCx[Nghost]         = dx;                                                     \
    atom->PBCy[Nghost]         = dy;                                                     \
    atom->PBCz[Nghost]         = dz;                                                     \
    atom->jclusters[cg].natoms = cj_natoms;                                              \
    Nghost_atoms += cj_natoms;                                                           \
    int cj_sca_base = CJ_SCALAR_BASE_INDEX(cj);                                          \
    int cg_sca_base = CJ_SCALAR_BASE_INDEX(cg);                                          \
    for (int cjj = 0; cjj < cj_natoms; cjj++) {                                          \
        atom->cl_t[cg_sca_base + cjj] = atom->cl_t[cj_sca_base + cjj];                   \
    }

/* internal subroutines */
void growPbc(Atom* atom) {
    int nold = NmaxGhost;
    NmaxGhost += DELTA;

    atom->border_map = (int*)reallocate(atom->border_map,
        ALIGNMENT,
        NmaxGhost * sizeof(int),
        nold * sizeof(int));
    atom->PBCx       = (int*)
        reallocate(atom->PBCx, ALIGNMENT, NmaxGhost * sizeof(int), nold * sizeof(int));
    atom->PBCy = (int*)
        reallocate(atom->PBCy, ALIGNMENT, NmaxGhost * sizeof(int), nold * sizeof(int));
    atom->PBCz = (int*)
        reallocate(atom->PBCz, ALIGNMENT, NmaxGhost * sizeof(int), nold * sizeof(int));
}

void setupPbc(Atom* atom, Parameter* param) {
    DEBUG_MESSAGE("setupPbc start\n");
    MD_FLOAT xprd     = param->xprd;
    MD_FLOAT yprd     = param->yprd;
    MD_FLOAT zprd     = param->zprd;
    MD_FLOAT cutNeigh = param->cutneigh;
    int jfac          = MAX(1, CLUSTER_N / CLUSTER_M);
    int ncj           = get_ncj_from_nci(atom->Nclusters_local);
    int Nghost        = -1;
    int Nghost_atoms  = 0;

    for (int cj = 0; cj < ncj; cj++) {
        if (atom->jclusters[cj].natoms > 0) {
            if (atom->Nclusters_local + (Nghost + 7) * jfac >= atom->Nclusters_max) {
                growClusters(atom, param->super_clustering);
            }

            if ((Nghost + 7) * jfac >= NmaxGhost) {
                growPbc(atom);
            }

            MD_FLOAT bbminx = atom->jclusters[cj].bbminx;
            MD_FLOAT bbmaxx = atom->jclusters[cj].bbmaxx;
            MD_FLOAT bbminy = atom->jclusters[cj].bbminy;
            MD_FLOAT bbmaxy = atom->jclusters[cj].bbmaxy;
            MD_FLOAT bbminz = atom->jclusters[cj].bbminz;
            MD_FLOAT bbmaxz = atom->jclusters[cj].bbmaxz;

            /* Setup ghost atoms */
            /* 6 planes */
            if (bbminx < cutNeigh) {
                ADDGHOST(+1, 0, 0);
            }
            if (bbmaxx >= (xprd - cutNeigh)) {
                ADDGHOST(-1, 0, 0);
            }
            if (bbminy < cutNeigh) {
                ADDGHOST(0, +1, 0);
            }
            if (bbmaxy >= (yprd - cutNeigh)) {
                ADDGHOST(0, -1, 0);
            }
            if (bbminz < cutNeigh) {
                ADDGHOST(0, 0, +1);
            }
            if (bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(0, 0, -1);
            }
            /* 8 corners */
            if (bbminx < cutNeigh && bbminy < cutNeigh && bbminz < cutNeigh) {
                ADDGHOST(+1, +1, +1);
            }
            if (bbminx < cutNeigh && bbmaxy >= (yprd - cutNeigh) && bbminz < cutNeigh) {
                ADDGHOST(+1, -1, +1);
            }
            if (bbminx < cutNeigh && bbminy < cutNeigh && bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(+1, +1, -1);
            }
            if (bbminx < cutNeigh && bbmaxy >= (yprd - cutNeigh) &&
                bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(+1, -1, -1);
            }
            if (bbmaxx >= (xprd - cutNeigh) && bbminy < cutNeigh && bbminz < cutNeigh) {
                ADDGHOST(-1, +1, +1);
            }
            if (bbmaxx >= (xprd - cutNeigh) && bbmaxy >= (yprd - cutNeigh) &&
                bbminz < cutNeigh) {
                ADDGHOST(-1, -1, +1);
            }
            if (bbmaxx >= (xprd - cutNeigh) && bbminy < cutNeigh &&
                bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(-1, +1, -1);
            }
            if (bbmaxx >= (xprd - cutNeigh) && bbmaxy >= (yprd - cutNeigh) &&
                bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(-1, -1, -1);
            }
            /* 12 edges */
            if (bbminx < cutNeigh && bbminz < cutNeigh) {
                ADDGHOST(+1, 0, +1);
            }
            if (bbminx < cutNeigh && bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(+1, 0, -1);
            }
            if (bbmaxx >= (xprd - cutNeigh) && bbminz < cutNeigh) {
                ADDGHOST(-1, 0, +1);
            }
            if (bbmaxx >= (xprd - cutNeigh) && bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(-1, 0, -1);
            }
            if (bbminy < cutNeigh && bbminz < cutNeigh) {
                ADDGHOST(0, +1, +1);
            }
            if (bbminy < cutNeigh && bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(0, +1, -1);
            }
            if (bbmaxy >= (yprd - cutNeigh) && bbminz < cutNeigh) {
                ADDGHOST(0, -1, +1);
            }
            if (bbmaxy >= (yprd - cutNeigh) && bbmaxz >= (zprd - cutNeigh)) {
                ADDGHOST(0, -1, -1);
            }
            if (bbminy < cutNeigh && bbminx < cutNeigh) {
                ADDGHOST(+1, +1, 0);
            }
            if (bbminy < cutNeigh && bbmaxx >= (xprd - cutNeigh)) {
                ADDGHOST(-1, +1, 0);
            }
            if (bbmaxy >= (yprd - cutNeigh) && bbminx < cutNeigh) {
                ADDGHOST(+1, -1, 0);
            }
            if (bbmaxy >= (yprd - cutNeigh) && bbmaxx >= (xprd - cutNeigh)) {
                ADDGHOST(-1, -1, 0);
            }
        }
    }

    if (ncj + (Nghost + 1) * jfac >= atom->Nclusters_max) {
        growClusters(atom, param->super_clustering);
    }

    // Add dummy cluster at the end
    int cjScaBase = CJ_SCALAR_BASE_INDEX(ncj + Nghost + 1);
    int cjVecBase = CJ_VECTOR_BASE_INDEX(ncj + Nghost + 1);
    int* cjT      = &atom->cl_t[cjScaBase];
    MD_FLOAT* cjX = &atom->cl_x[cjVecBase];

    for (int cjj = 0; cjj < CLUSTER_N; cjj++) {
        cjX[CL_X_INDEX(cjj)] = INFINITY;
        cjX[CL_Y_INDEX(cjj)] = INFINITY;
        cjX[CL_Z_INDEX(cjj)] = INFINITY;
        cjT[cjj]               = 0;
    }

    // Increase by one to make it the ghost atom count
    atom->dummy_cj        = ncj + Nghost + 1;
    atom->Nghost          = Nghost_atoms;
    atom->Nclusters_ghost = Nghost + 1;

    // Update created ghost clusters positions
    updatePbcCPU(atom, param, 1);
    DEBUG_MESSAGE("setupPbc end\n");
}
