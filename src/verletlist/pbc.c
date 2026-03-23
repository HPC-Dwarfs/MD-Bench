/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <allocate.h>
#include <atom.h>
#include <pbc.h>

#define DELTA 20000

int NmaxGhost;
int *PBCx, *PBCy, *PBCz;

#ifdef CUDA_TARGET
UpdatePbcFunction updatePbc      = updatePbcCUDA;
UpdatePbcFunction updateAtomsPbc = updateAtomsPbcCUDA;
#else
UpdatePbcFunction updatePbc      = updatePbcCPU;
UpdatePbcFunction updateAtomsPbc = updateAtomsPbcCPU;
#endif

static void growPbc(Atom*);

/* exported subroutines */
void initPbc(Atom* atom)
{
    NmaxGhost        = 0;
    atom->border_map = NULL;
    PBCx             = NULL;
    PBCy             = NULL;
    PBCz             = NULL;
}

/* update coordinates of ghost atoms */
/* uses mapping created in setupPbc */
void updatePbcCPU(Atom* atom, Parameter* param, bool doReneighbor)
{
    int* borderMap = atom->border_map;
    int nlocal     = atom->Nlocal;
    MD_FLOAT xprd  = param->xprd;
    MD_FLOAT yprd  = param->yprd;
    MD_FLOAT zprd  = param->zprd;

    for (int i = 0; i < atom->Nghost; i++) {
        atom_x(nlocal + i) = atom_x(borderMap[i]) + PBCx[i] * xprd;
        atom_y(nlocal + i) = atom_y(borderMap[i]) + PBCy[i] * yprd;
        atom_z(nlocal + i) = atom_z(borderMap[i]) + PBCz[i] * zprd;
        // Copy per-atom LJ parameters for geometric combination rule
        atom->sqrt_epsilon[nlocal + i] = atom->sqrt_epsilon[borderMap[i]];
        atom->sigma3[nlocal + i] = atom->sigma3[borderMap[i]];
    }
}

/* relocate atoms that have left domain according
 * to periodic boundary conditions */
void updateAtomsPbcCPU(Atom* atom, Parameter* param, bool doReneighbor)
{
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
    border_map[Nghost]                = i;                                               \
    PBCx[Nghost]                      = dx;                                              \
    PBCy[Nghost]                      = dy;                                              \
    PBCz[Nghost]                      = dz;                                              \
    atom->type[atom->Nlocal + Nghost] = atom->type[i];                                   \
    atom->sqrt_epsilon[atom->Nlocal + Nghost] = atom->sqrt_epsilon[i];                  \
    atom->sigma3[atom->Nlocal + Nghost] = atom->sigma3[i]

void setupPbc(Atom* atom, Parameter* param)
{
    int* border_map   = atom->border_map;
    MD_FLOAT xprd     = param->xprd;
    MD_FLOAT yprd     = param->yprd;
    MD_FLOAT zprd     = param->zprd;
    MD_FLOAT cutneigh = param->cutneigh;
    int Nghost        = -1;

    for (int i = 0; i < atom->Nlocal; i++) {
        if (atom->Nlocal + Nghost + 7 >= atom->Nmax) {
            growAtom(atom);
        }

        if (Nghost + 7 >= NmaxGhost) {
            growPbc(atom);
            border_map = atom->border_map;
        }

        MD_FLOAT x = atom_x(i);
        MD_FLOAT y = atom_y(i);
        MD_FLOAT z = atom_z(i);

        /* Setup ghost atoms */
        /* 6 planes */
        if (param->pbc_x != 0) {
            if (x < cutneigh) {
                ADDGHOST(+1, 0, 0);
            }
            if (x >= (xprd - cutneigh)) {
                ADDGHOST(-1, 0, 0);
            }
        }

        if (param->pbc_y != 0) {
            if (y < cutneigh) {
                ADDGHOST(0, +1, 0);
            }
            if (y >= (yprd - cutneigh)) {
                ADDGHOST(0, -1, 0);
            }
        }

        if (param->pbc_z != 0) {
            if (z < cutneigh) {
                ADDGHOST(0, 0, +1);
            }
            if (z >= (zprd - cutneigh)) {
                ADDGHOST(0, 0, -1);
            }
        }

        /* 8 corners */
        if (param->pbc_x != 0 && param->pbc_y != 0 && param->pbc_z != 0) {
            if (x < cutneigh && y < cutneigh && z < cutneigh) {
                ADDGHOST(+1, +1, +1);
            }
            if (x < cutneigh && y >= (yprd - cutneigh) && z < cutneigh) {
                ADDGHOST(+1, -1, +1);
            }
            if (x < cutneigh && y < cutneigh && z >= (zprd - cutneigh)) {
                ADDGHOST(+1, +1, -1);
            }
            if (x < cutneigh && y >= (yprd - cutneigh) && z >= (zprd - cutneigh)) {
                ADDGHOST(+1, -1, -1);
            }
            if (x >= (xprd - cutneigh) && y < cutneigh && z < cutneigh) {
                ADDGHOST(-1, +1, +1);
            }
            if (x >= (xprd - cutneigh) && y >= (yprd - cutneigh) && z < cutneigh) {
                ADDGHOST(-1, -1, +1);
            }
            if (x >= (xprd - cutneigh) && y < cutneigh && z >= (zprd - cutneigh)) {
                ADDGHOST(-1, +1, -1);
            }
            if (x >= (xprd - cutneigh) && y >= (yprd - cutneigh) &&
                z >= (zprd - cutneigh)) {
                ADDGHOST(-1, -1, -1);
            }
        }

        /* 12 edges */
        if (param->pbc_x != 0 && param->pbc_z != 0) {
            if (x < cutneigh && z < cutneigh) {
                ADDGHOST(+1, 0, +1);
            }
            if (x < cutneigh && z >= (zprd - cutneigh)) {
                ADDGHOST(+1, 0, -1);
            }
            if (x >= (xprd - cutneigh) && z < cutneigh) {
                ADDGHOST(-1, 0, +1);
            }
            if (x >= (xprd - cutneigh) && z >= (zprd - cutneigh)) {
                ADDGHOST(-1, 0, -1);
            }
        }

        if (param->pbc_y != 0 && param->pbc_z != 0) {
            if (y < cutneigh && z < cutneigh) {
                ADDGHOST(0, +1, +1);
            }
            if (y < cutneigh && z >= (zprd - cutneigh)) {
                ADDGHOST(0, +1, -1);
            }
            if (y >= (yprd - cutneigh) && z < cutneigh) {
                ADDGHOST(0, -1, +1);
            }
            if (y >= (yprd - cutneigh) && z >= (zprd - cutneigh)) {
                ADDGHOST(0, -1, -1);
            }
        }

        if (param->pbc_x != 0 && param->pbc_y != 0) {
            if (y < cutneigh && x < cutneigh) {
                ADDGHOST(+1, +1, 0);
            }
            if (y < cutneigh && x >= (xprd - cutneigh)) {
                ADDGHOST(-1, +1, 0);
            }
            if (y >= (yprd - cutneigh) && x < cutneigh) {
                ADDGHOST(+1, -1, 0);
            }
            if (y >= (yprd - cutneigh) && x >= (xprd - cutneigh)) {
                ADDGHOST(-1, -1, 0);
            }
        }
    }
    // increase by one to make it the ghost atom count
    atom->Nghost = Nghost + 1;
}

/* internal subroutines */
void growPbc(Atom* atom)
{
    int nold = NmaxGhost;
    NmaxGhost += DELTA;

    atom->border_map = (int*)reallocate(atom->border_map,
        ALIGNMENT,
        NmaxGhost * sizeof(int),
        nold * sizeof(int));
    PBCx = (int*)reallocate(PBCx, ALIGNMENT, NmaxGhost * sizeof(int), nold * sizeof(int));
    PBCy = (int*)reallocate(PBCy, ALIGNMENT, NmaxGhost * sizeof(int), nold * sizeof(int));
    PBCz = (int*)reallocate(PBCz, ALIGNMENT, NmaxGhost * sizeof(int), nold * sizeof(int));
}
