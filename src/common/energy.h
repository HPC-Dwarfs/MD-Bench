/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// Energy output of the force kernels (-energy): pair sums plus Coulomb self energy
#ifndef __ENERGY_H_
#define __ENERGY_H_

#include <stdio.h>

#ifdef _MPI
#include <mpi.h>
#endif

#include <coulomb_scalar.h>
#include <parameter.h>

// Sum of q^2 of the local atoms (q per atom type)
static inline double localChargeSquaredSum(
    const int nlocal, const int* type, const MD_FLOAT* q_per_type)
{
    double sum = 0.0;
    if (q_per_type == NULL) {
        return sum;
    }
    for (int i = 0; i < nlocal; i++) {
        const double q = q_per_type[type[i]];
        sum += q * q;
    }
    return sum;
}

// Stores the energies of a force call of this process in param
static inline void storeEnergies(Parameter* param, double elj, double ecoul, double sumq2)
{
    const CoulombScalar cs = coulombScalarInit(param);
    param->energy_lj       = elj;
    param->energy_coul     = ecoul - param->epsfac * sumq2 * coulombSelfEnergy(&cs);
}

// Energies of the whole system (summed over the MPI processes)
static inline void globalEnergies(const Parameter* param, double* elj, double* ecoul)
{
    double e[2] = { param->energy_lj, param->energy_coul };
#ifdef _MPI
    MPI_Allreduce(MPI_IN_PLACE, e, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    *elj   = e[0];
    *ecoul = e[1];
}

#endif
