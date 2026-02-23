/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <parameter.h>
#ifdef _MPI
#include <mpi.h>
#endif

#ifndef __NEIGHBOR_H_
#define __NEIGHBOR_H_

#ifdef NBLIST_AOS
#define NBLIST_DATA_LAYOUT "AoS"
#define neighs(nblist,i,j,M,N) nblist[(i) * N + (j)]
#else
#define NBLIST_DATA_LAYOUT "SoA"
#define neighs(nblist,i,j,M,N) nblist[(j) * M + (i)]
#endif

typedef struct {
    int* neighbors;
    int* numneigh;
} DeviceNeighbor;

typedef struct {
    int every;
    int ncalls;
    int maxneighs;
    int half_neigh;
    int* neighbors;
    int* numneigh;

    // Device data
    DeviceNeighbor d_neighbor;
    // MPI
    int half_stencil;
    int Nshell;         // # of atoms in listShell
    int* numNeighShell; // # of neighs for each atom in listShell
    int* neighshell;    // list of neighs for each atom in listShell
    int* listshell;     // Atoms to compute the force

} Neighbor;

typedef struct {
    MD_FLOAT xprd;
    MD_FLOAT yprd;
    MD_FLOAT zprd;
    MD_FLOAT bininvx;
    MD_FLOAT bininvy;
    MD_FLOAT bininvz;
    int mbinxlo;
    int mbinylo;
    int mbinzlo;
    int nbinx;
    int nbiny;
    int nbinz;
    int mbinx;
    int mbiny;
    int mbinz;
    // Multigpu
    int pad_x;
    int pad_y;
    int pad_z;
    MD_FLOAT binsizex;
    MD_FLOAT binsizey;
    MD_FLOAT binsizez;
} Neighbor_params;

typedef struct {
    int* bincount;
    int* bins;
    int mbins;
    int atoms_per_bin;
} Binning;

typedef void (*BuildNeighborFunction)(Atom*, Neighbor*);
extern BuildNeighborFunction buildNeighbor;

extern void initNeighbor(Neighbor*, Parameter*);
extern void setupNeighbor(Parameter*);
extern void binatoms(Atom*);
extern void sortAtom(Atom*);
extern void buildNeighborCPU(Atom*, Neighbor*);
#ifdef CUDA_TARGET
#ifdef __cplusplus
extern "C"
#endif
extern void buildNeighborCUDA(Atom*, Neighbor*);
#endif
#endif //__NEIGHBOR_H_
