/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#ifndef __TIMERS_H_
#define __TIMERS_H_

// typedef enum { TOTAL = 0, NEIGH, FORCE, NUMTIMER } timertype;

typedef enum {
    TOTAL = 0,
    FORCE,
    NEIGH,
    FORWARD,
    REVERSE,
    UPDATE,
    BALANCE,
    SETUP,
    REST,
    // Fine-grained per-routine breakdown, only printed with -v/--verbose.
    // These are sub-components of the buckets above (e.g. NEIGH_PBC +
    // NEIGH_BUILD + NEIGH_PRUNE [+ NEIGH_SORT | NEIGH_CLUSTERS + NEIGH_BIN]
    // sum to NEIGH), always accumulated since the extra getTimeStamp() calls
    // are negligible next to a simulation step.
    NEIGH_PBC,      // setupPbc()/updatePbc() (verletlist) or ghostNeighbor() (MPI)
    NEIGH_BUILD,    // buildNeighbor()
    NEIGH_PRUNE,    // pruneNeighbor()/pruneNeighborCUDA()
    NEIGH_SORT,     // sortAtom() (verletlist, SORT_ATOMS only)
    NEIGH_CLUSTERS, // buildClusters()+defineJClusters() (clusterpair only)
    NEIGH_BIN,      // binJClusters() (clusterpair only)
    INTEGRATE_INITIAL, // initialIntegrate()
    INTEGRATE_FINAL,   // finalIntegrate()
    THERMO,         // computeThermo()
    RENEIGH_CHECK,  // needsReneigh() (displacement_reneigh only)
    NUMTIMER
} timertype;

#endif
