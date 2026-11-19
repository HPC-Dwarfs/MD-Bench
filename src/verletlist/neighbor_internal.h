/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// State and helpers of neighbor.c shared with the list-build template neighbor_build.c
#ifndef __NEIGHBOR_INTERNAL_H_
#define __NEIGHBOR_INTERNAL_H_

#include <math.h>

#include <atom.h>
#include <kernel_variant.h>
#include <neighbor.h>
#include <parameter.h>
#include <util.h>

extern MD_FLOAT bininvx, bininvy, bininvz;
extern int mbinx, mbiny, mbinz;
extern int* bincount;
extern int* bins;
extern int atoms_per_bin;
extern MD_FLOAT cutneighsq;
extern MD_FLOAT cutneigh_inner_sq;
extern int dcut_enabled;
extern int* is_inner_buf;
extern int nstencil;
extern int* stencil;
extern MD_FLOAT binsizex, binsizey, binsizez;
extern int pad_x, pad_y, pad_z;
extern int me;
extern int method;
extern int half_stencil;

int ghostZone(Atom*, int);
int eightZone(Atom*, int);
int halfZone(Atom*, int);
DECLARE_RULE_VARIANT_FUNCTIONS(void, buildNeighborLists, Atom*, Neighbor*)
DECLARE_RULE_VARIANT_FUNCTIONS(void, neighborGhost, Atom*, Neighbor*)

static inline int coord2bin(MD_FLOAT xin, MD_FLOAT yin, MD_FLOAT zin)
{
    int ix, iy, iz;
    MD_FLOAT eps = 1e-9;
    MD_FLOAT xlo = 0.0;
    MD_FLOAT ylo = 0.0;
    MD_FLOAT zlo = 0.0;
    xlo          = fabs(xlo - pad_x * binsizex) + eps;
    ylo          = fabs(ylo - pad_y * binsizey) + eps;
    zlo          = fabs(zlo - pad_z * binsizez) + eps;
    ix           = (int)((xin + xlo) * bininvx);
    iy           = (int)((yin + ylo) * bininvy);
    iz           = (int)((zin + zlo) * bininvz);

    return (iz * mbiny * mbinx + iy * mbinx + ix);

    /*
    if (xin >= xprd) {
        ix = (int)((xin - xprd) * bininvx) + nbinx - mbinxlo;
    } else if (xin >= 0.0) {
        ix = (int)(xin * bininvx) - mbinxlo;
    } else {
        ix = (int)(xin * bininvx) - mbinxlo - 1;
    }

    if (yin >= yprd) {
        iy = (int)((yin - yprd) * bininvy) + nbiny - mbinylo;
    } else if (yin >= 0.0) {
        iy = (int)(yin * bininvy) - mbinylo;
    } else {
        iy = (int)(yin * bininvy) - mbinylo - 1;
    }

    if (zin >= zprd) {
        iz = (int)((zin - zprd) * bininvz) + nbinz - mbinzlo;
    } else if (zin >= 0.0) {
        iz = (int)(zin * bininvz) - mbinzlo;
    } else {
        iz = (int)(zin * bininvz) - mbinzlo - 1;
    }

    return (iz * mbiny * mbinx + iy * mbinx + ix + 1);
    */
}

static inline int skipNeigh(Atom* atom, int i, int j)
{
    if (i > j && j < atom->Nlocal) {
        return 1;
    } else if (atom_z(i) > atom_z(j) && j >= atom->Nlocal) {
        return 1;
    } else if (Equal(atom_z(i), atom_z(j)) && atom_y(i) > atom_y(j) &&
               j >= atom->Nlocal) {
        return 1;
    } else if (Equal(atom_z(i), atom_z(j)) && Equal(atom_y(i), atom_y(j)) &&
               atom_x(i) > atom_x(j) && j >= atom->Nlocal) {
        return 1;
    } else {
        return 0;
    }
}

#endif
