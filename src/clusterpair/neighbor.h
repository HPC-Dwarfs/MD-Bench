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
// GPU kernels only receive a flat neighbor buffer plus maxneighs (no Neighbor*),
// so the neighs() macro's AOS/CSR branches (which index through nbr->maxneighs /
// nbr->neigh_start) cannot be used there; only the SoA branch works.
#if defined(CUDA_TARGET) && (defined(NBLIST_AOS) || defined(NBLIST_CSR))
#error "NBLIST_DATA_LAYOUT must be SOA (or auto) for GPU targets (NVCC/HIPCC)"
#endif
// Interaction masks from GROMACS, things to remember (maybe these confused just me):
//   1. These are not "exclusion" masks as the name suggests in GROMACS, but rather
//      interaction masks (1 = interaction, 0 = no interaction)
//   2. These are inverted (maybe because that is how you use in AVX2/AVX512 masking),
//      so read them from right to left (least significant to most significant bit)
// All interaction mask is the same for all kernels
#define NBNXN_INTERACTION_MASK_ALL 0xffffffffU
// 4x4 kernel diagonal mask
#define NBNXN_INTERACTION_MASK_DIAG 0x08ceU
// 4x2 kernel diagonal masks
#define NBNXN_INTERACTION_MASK_DIAG_J2_0 0x0002U
#define NBNXN_INTERACTION_MASK_DIAG_J2_1 0x002fU
// 4x8 kernel diagonal masks
#define NBNXN_INTERACTION_MASK_DIAG_J8_0 0xf0f8fcfeU
#define NBNXN_INTERACTION_MASK_DIAG_J8_1 0x0080c0e0U

#ifdef NBLIST_AOS
#define NBLIST_DATA_LAYOUT      "AoS"
#define neighs(l, i, j, M, nbr) (l)[(i) * (nbr)->maxneighs + (j)]
#elif defined(NBLIST_CSR)
#define NBLIST_DATA_LAYOUT      "CSR"
#define neighs(l, i, j, M, nbr) (l)[(nbr)->neigh_start[(i)] + (j)]
#else
#define NBLIST_DATA_LAYOUT      "SoA"
#define neighs(l, i, j, M, nbr) (l)[(j) * (M) + (i)]
#endif
/* Shell list and build-phase scratch always use padded AOS layout */
#define neighshell(nblist, i, j, N)       nblist[(i) * (N) + (j)]
#define neighs_padded(nblist, i, j, M, N) nblist[(i) * (N) + (j)]

typedef struct {
    int every;
    int ncalls;
    int maxneighs;
    int* numneigh;
    int* numneigh_masked;
    int* numneigh_inner;
    int* numneigh_inner_masked;
    int half_neigh;
    int* neighbors;
    int* neigh_start;
    unsigned int* neighbors_imask;

#ifdef NBLIST_PAIRLIST
    // Explicit flat (ci,cj) pair list, flattened from the per-ci list above
    // right after pruning, preserving masked-prefix/unmasked-suffix order.
    // Padded to a multiple of VECTOR_WIDTH; padding lanes have
    // pair_ci = pair_cj = -1 and pair_mask = 0 (skipped explicitly).
    int npairs;
    int nchunks;
    int* pair_ci;
    int* pair_cj;
    unsigned int* pair_mask;
#endif

    // MPI
    int Nshell;
    int* numNeighShell;
    int* neighshell;
    int* listshell;
} Neighbor;

typedef void (*BuildNeighborFunction)(Atom*, Neighbor*);
typedef void (*PruneNeighborFunction)(Parameter*, Atom*, Neighbor*);
typedef void (*BuildClustersFunction)(Atom*);
extern BuildNeighborFunction buildNeighbor;
extern PruneNeighborFunction pruneNeighbor;
extern BuildClustersFunction buildClusters;

extern void initNeighbor(Neighbor*, Parameter*);
extern void setupNeighbor(Parameter*, Atom*);
extern void binatoms(Atom*);
extern void buildNeighborCPU(Atom*, Neighbor*);
extern void buildNeighborSuperclusters(Atom*, Neighbor*);
extern void pruneNeighborCPU(Parameter*, Atom*, Neighbor*);
extern void pruneNeighborSuperclusters(Parameter*, Atom*, Neighbor*);
extern void buildClustersCPU(Atom*);
extern void buildSuperclusters(Atom*);
extern void defineJClusters(Parameter*, Atom*);
extern void binJClusters(Parameter*, Atom*);
extern void updateSingleAtoms(Parameter*, Atom*);

#ifdef CUDA_TARGET
#ifdef __cplusplus
extern "C"
#endif
    extern void
    growNeighborCUDA(Atom*, Neighbor*);
#endif

#endif
