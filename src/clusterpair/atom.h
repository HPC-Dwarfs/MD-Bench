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

/*
    A little workaround to get also newer versions of ICX
    to compile the code with a valid definition for INFINITY.
    Also should support compiling with -Ofast, which falls
    back to finit-math-only and renders the INFINITY macro
    as undefinde value.
*/
#if !defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
#include <math.h>
#if defined(MD_FLOAT) && MD_FLOAT == double
#include <float.h>
#define INF DBL_MAX
#else
#include <float.h>
#define INF FLT_MAX
#endif
#else
#include <math.h>
#define INF INFINITY
#endif

#define DELTA 100000

typedef struct {
    int natoms;
    MD_FLOAT bbminx, bbmaxx;
    MD_FLOAT bbminy, bbmaxy;
    MD_FLOAT bbminz, bbmaxz;
} Cluster;

typedef struct {
    int nclusters;
    MD_FLOAT bbminx, bbmaxx;
    MD_FLOAT bbminy, bbmaxy;
    MD_FLOAT bbminz, bbmaxz;
} SuperCluster;

typedef struct {
    int Natoms, Nlocal, Nghost, Nmax;
    int Nclusters_local, Nclusters_ghost, Nclusters_max, NmaxGhost, ncj;
    MD_FLOAT *x, *y, *z;
    MD_FLOAT *vx, *vy, *vz;
    int* border_map;
    int* type;
    int ntypes;
    MD_FLOAT* epsilon;
    MD_FLOAT* sigma6;
    MD_FLOAT* cutforcesq;
    MD_FLOAT* cutneighsq;
    // Per-type LJ parameters for filling cluster arrays
    MD_FLOAT* sqrt_epsilon_per_type;
    MD_FLOAT* sigma3_per_type;
    int *PBCx, *PBCy, *PBCz;
    // Data in cluster format
    MD_FLOAT* cl_x;
    MD_FLOAT* cl_v;
    MD_FLOAT* cl_f;
    int* cl_t;
    // Per-cluster LJ parameters for geometric combination (Gromacs-style optimization)
    MD_FLOAT* cl_sqrt_epsilon; // sqrt(epsilon) per atom in cluster layout
    MD_FLOAT* cl_sigma3;       // sigma^3 per atom in cluster layout
    Cluster *iclusters, *jclusters;
    SuperCluster* siclusters;
    int* cluster_bin;
    int dummy_cj;
    MD_UINT* exclusion_filter;
    MD_FLOAT* diagonal_4xn_j_minus_i;
    MD_FLOAT* diagonal_2xnn_j_minus_i;
    unsigned int masks_2xnn_hn[8];
    unsigned int masks_2xnn_fn[8];
    unsigned int masks_4xn_hn[16];
    unsigned int masks_4xn_fn[16];
    // Info Subdomain
    Box mybox;
} Atom;

extern int get_ncj_from_nci(int nci);
extern void initAtom(Atom*);
extern void initMasks(Atom*);
extern void createAtom(Atom*, Parameter*);
extern int readAtom(Atom*, Parameter*);
extern int readAtomPdb(Atom*, Parameter*);
extern int readAtomGro(Atom*, Parameter*);
extern int readAtomDmp(Atom*, Parameter*);
extern void growAtom(Atom*);
extern void freeAtom(Atom*);
extern void growClusters(Atom*, int);

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

#ifdef CUDA_TARGET
#ifdef __cplusplus
extern "C"
#endif
    extern void
    growClustersCUDA(Atom*);
#endif

#ifdef ATOM_POSITION_AOS
#define POS_DATA_LAYOUT "AoS"
#define atom_x(i)       atom->x[(i) * 3 + 0]
#define atom_y(i)       atom->x[(i) * 3 + 1]
#define atom_z(i)       atom->x[(i) * 3 + 2]
/*
#   define atom_vx(i)          atom->vx[(i) * 3 + 0]
#   define atom_vy(i)          atom->vx[(i) * 3 + 1]
#   define atom_vz(i)          atom->vx[(i) * 3 + 2]
#   define atom_fx(i)          atom->fx[(i) * 3 + 0]
#   define atom_fy(i)          atom->fx[(i) * 3 + 1]
#   define atom_fz(i)          atom->fx[(i) * 3 + 2]
*/
#else
#define POS_DATA_LAYOUT "SoA"
#define atom_x(i)       atom->x[i]
#define atom_y(i)       atom->y[i]
#define atom_z(i)       atom->z[i]
#endif

// TODO: allow to switch velocites and forces to AoS
#define atom_vx(i) atom->vx[i]
#define atom_vy(i) atom->vy[i]
#define atom_vz(i) atom->vz[i]
#define atom_fx(i) atom->fx[i]
#define atom_fy(i) atom->fy[i]
#define atom_fz(i) atom->fz[i]

#endif
