/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <likwid-marker.h>
#include <math.h>

#include <allocate.h>
#include <atom.h>
#include <eam.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>

double computeForceEam(Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats) {
    double S = getTimeStamp();
    LIKWID_MARKER_START("force");
    fprintf(stderr, "computeForceEam(): function not implemented for Cluster Pair algorithm!");
    exit(-1);
    LIKWID_MARKER_STOP("force");
    double E = getTimeStamp();
    return E - S;
}
