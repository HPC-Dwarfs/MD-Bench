/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <box.h>
#include <parameter.h>
#ifdef _MPI
#include <mpi.h>
#endif

#ifndef __ATOM_H_
#define __ATOM_H_

typedef struct {
    MD_FLOAT *x, *y, *z;
    MD_FLOAT *vx, *vy, *vz;
    MD_FLOAT *fx, *fy, *fz;
    int* border_map;
    int* type;
    MD_FLOAT* epsilon;
    MD_FLOAT* sigma6;
    MD_FLOAT* cutforcesq;
    MD_FLOAT* cutneighsq;
    // Per-atom LJ params for geometric mixing (avoids indirect gathers)
    MD_FLOAT* sqrt_epsilon;  // sqrt(epsilon) per atom
    MD_FLOAT* sigma3;        // sigma^3 per atom (sigma6 = sigma3_i * sigma3_j)
} DeviceAtom;

typedef struct {
    int Natoms, Nlocal, Nghost, Nmax;
    MD_FLOAT *x, *y, *z;
    MD_FLOAT *vx, *vy, *vz;
    MD_FLOAT *fx, *fy, *fz;
    int* border_map;
    int* type;
    int ntypes;
    MD_FLOAT* epsilon;
    MD_FLOAT* sigma6;
    MD_FLOAT* cutforcesq;
    MD_FLOAT* cutneighsq;
    // Per-atom LJ params for geometric mixing (avoids indirect gathers)
    MD_FLOAT* sqrt_epsilon;  // sqrt(epsilon) per atom
    MD_FLOAT* sigma3;        // sigma^3 per atom (sigma6 = sigma3_i * sigma3_j)
    // Per-type LJ params for filling per-atom arrays
    MD_FLOAT* sqrt_epsilon_per_type;
    MD_FLOAT* sigma3_per_type;

    // Device data
    DeviceAtom d_atom;
    // Info Subdomain
    Box mybox;
} Atom;

extern void initAtom(Atom*);
extern void createAtom(Atom*, Parameter*);
extern int readAtom(Atom*, Parameter*);
extern int readAtom_pdb(Atom*, Parameter*);
extern int readAtom_gro(Atom*, Parameter*);
extern int readAtom_dmp(Atom*, Parameter*);
extern int readAtom_in(Atom*, Parameter*);
extern void writeAtom(Atom*, Parameter*);
extern void growAtom(Atom*);
extern void freeAtom(Atom*);

int packGhost(Atom*, int, MD_FLOAT*, int*);
int unpackGhost(Parameter*, Atom*, int, MD_FLOAT*);
int packExchange(Atom*, int, MD_FLOAT*);
int unpackExchange(Atom*, int, MD_FLOAT*);
void packForward(Atom*, int, int*, MD_FLOAT*, int*);
void unpackForward(Atom*, int, int, MD_FLOAT*);
void packReverse(Atom*, int, int, MD_FLOAT*);
void unpackReverse(Atom*, int, int*, MD_FLOAT*);
void pbc(Atom*);
void copy(Atom*, int, int);

#ifdef ATOM_POSITION_AOS
#define POS_DATA_LAYOUT "AoS"
#define atom_x(i)       atom->x[(i)*3 + 0]
#define atom_y(i)       atom->x[(i)*3 + 1]
#define atom_z(i)       atom->x[(i)*3 + 2]
#define atom_vx(i)      atom->vx[(i)*3 + 0]
#define atom_vy(i)      atom->vx[(i)*3 + 1]
#define atom_vz(i)      atom->vx[(i)*3 + 2]
#define atom_fx(i)      atom->fx[(i)*3 + 0]
#define atom_fy(i)      atom->fx[(i)*3 + 1]
#define atom_fz(i)      atom->fx[(i)*3 + 2]
#else
#define POS_DATA_LAYOUT "SoA"
#define atom_x(i)       atom->x[i]
#define atom_y(i)       atom->y[i]
#define atom_z(i)       atom->z[i]
#define atom_vx(i)      atom->vx[i]
#define atom_vy(i)      atom->vy[i]
#define atom_vz(i)      atom->vz[i]
#define atom_fx(i)      atom->fx[i]
#define atom_fy(i)      atom->fy[i]
#define atom_fz(i)      atom->fz[i]
#endif

#define buf_x(i) buf[3 * (i)]
#define buf_y(i) buf[3 * (i) + 1]
#define buf_z(i) buf[3 * (i) + 2]

#endif // __ATOM_H_
