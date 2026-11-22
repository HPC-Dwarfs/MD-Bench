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
#define NBLIST_DATA_LAYOUT      "AoS"
#define neighs(l, i, j, M, nbr) (l)[(i) * (nbr)->maxneighs + (j)]
#elif defined(NBLIST_CSR)
#define NBLIST_DATA_LAYOUT      "CSR"
#define neighs(l, i, j, M, nbr) (l)[(nbr)->neigh_start[(i)] + (j)]
#else
/* SoA ("neighbor"-major). GPU lists built for T = 2^tpa_log2 lanes per atom (TPA) keep the
 * neighbors of an atom in groups of T (LAMMPS GPU package layout), so the lanes of an atom and
 * the atoms of a warp read consecutive words; tpa_log2 = 0 is the plain layout j * M + i. */
#define NBLIST_DATA_LAYOUT "SoA"
#define neighs(l, i, j, M, nbr)                                                             \
    (l)[((((j) >> (nbr)->tpa_log2) * (M)) << (nbr)->tpa_log2) + ((i) << (nbr)->tpa_log2) + \
        ((j) & ((1 << (nbr)->tpa_log2) - 1))]
#endif
/* Shell list and build-phase scratch always use padded AOS layout */
#define neighshell(nblist, i, j, N)       nblist[(i) * (N) + (j)]
#define neighs_padded(nblist, i, j, M, N) nblist[(i) * (N) + (j)]

/* GPU force kernels: lanes per atom (TPA), fixed per list type as LAMMPS' `package gpu tpa`.
 * Given as log2 because they also set the layout of the SOA lists (see neighs). */
#define VL_TPA_LOG2_HALF 3 // 8 lanes per atom with half lists
#define VL_TPA_LOG2_FULL 4 // 16 lanes per atom with full lists

typedef struct {
    int* neighbors;
    int* numneigh;
    int* numneigh_inner;
    int maxneighs;
    int tpa_log2; // lanes per atom (TPA) the list is laid out for, see neighs
} DeviceNeighbor;

typedef struct {
    int every;
    int ncalls;
    int maxneighs;
    int tpa_log2; // lanes per atom (TPA) the list is laid out for, see neighs; 0 on the CPU
    int half_neigh;
    int* neighbors;
    int* neigh_start;
    int* numneigh;
    int* numneigh_inner;

    // Device data
    DeviceNeighbor d_neighbor;
    // MPI
    int half_stencil;
    int Nshell;         // # of atoms in listShell
    int* numNeighShell; // # of neighs for each atom in listShell
    int* neighshell;    // list of neighs for each atom in listShell
    int* listshell;     // Atoms to compute the force

    // LJ combination rule: the variant of the list build (per-type or global cut-off)
    int lj_comb_rule;
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
typedef void (*PruneNeighborFunction)(Parameter*, Atom*, Neighbor*);
extern BuildNeighborFunction buildNeighbor;
extern PruneNeighborFunction pruneNeighbor;

extern void initNeighbor(Neighbor*, Parameter*);
extern void setupNeighbor(Parameter*);
extern void binatoms(Atom*);
extern void sortAtom(Atom*);
extern void buildNeighborCPU(Atom*, Neighbor*);
extern void pruneNeighborCPU(Parameter*, Atom*, Neighbor*);
#ifdef CUDA_TARGET
#ifdef __cplusplus
extern "C"
#endif
    extern void
    buildNeighborCUDA(Atom*, Neighbor*);
#ifdef __cplusplus
extern "C"
#endif
    extern void pruneNeighborCUDA(Parameter*, Atom*, Neighbor*);
#ifdef __cplusplus
extern "C"
#endif
    extern void sortAtomsCUDA(Atom*);
#endif
#endif //__NEIGHBOR_H_
