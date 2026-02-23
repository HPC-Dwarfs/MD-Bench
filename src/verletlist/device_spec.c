/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <device.h>
#include <util.h>

#ifdef CUDA_TARGET

void initDevice(Parameter* param, Atom* atom, Neighbor* neighbor) {
    DEBUG_MESSAGE("initDevice start\n");

    DeviceAtom* d_atom         = &(atom->d_atom);
    DeviceNeighbor* d_neighbor = &(neighbor->d_neighbor);

    d_atom->epsilon = (MD_FLOAT*)allocateGPU(
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    d_atom->sigma6 = (MD_FLOAT*)allocateGPU(
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    d_atom->cutneighsq = (MD_FLOAT*)allocateGPU(
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    d_atom->cutforcesq = (MD_FLOAT*)allocateGPU(
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    d_neighbor->neighbors = (int*)allocateGPU(
        sizeof(int) * atom->Nmax * neighbor->maxneighs);
    d_neighbor->numneigh = (int*)allocateGPU(sizeof(int) * atom->Nmax);

    memcpyToGPU(d_atom->x, atom->x, sizeof(MD_FLOAT) * atom->Nmax * 3);
    memcpyToGPU(d_atom->vx, atom->vx, sizeof(MD_FLOAT) * atom->Nmax * 3);
    memcpyToGPU(d_atom->sigma6,
        atom->sigma6,
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    memcpyToGPU(d_atom->epsilon,
        atom->epsilon,
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    memcpyToGPU(d_atom->cutneighsq,
        atom->cutneighsq,
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    memcpyToGPU(d_atom->cutforcesq,
        atom->cutforcesq,
        sizeof(MD_FLOAT) * atom->ntypes * atom->ntypes);
    memcpyToGPU(d_atom->type, atom->type, sizeof(int) * atom->Nmax);

    DEBUG_MESSAGE("initDevice stop\n");
}

#endif
