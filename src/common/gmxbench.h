/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
/* The water system of the GROMACS non-bonded benchmark (-gmxbench), shared by the
 * cluster-pair and verlet-list schemes. Include it after the scheme's atom.h: it fills
 * the Atom fields both schemes have (positions, velocities, types, the per-type LJ
 * parameters and charges, molsize). */
#ifndef __GMXBENCH_H_
#define __GMXBENCH_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <allocate.h>
#include <gmxcoords.h>
#include <parameter.h>
#include <util.h>

// Water model of the GROMACS non-bonded benchmark: atoms per molecule (O, H, H), atom
// types, LJ C6/C12 of the oxygen (as given to GROMACS) and partial charges
#define WATER_ATOMS      3
#define WATER_NTYPES     2
#define WATER_C6_OXYGEN  0.0026173456
#define WATER_C12_OXYGEN 2.634129e-06
#define WATER_Q_OXYGEN   -0.8476
#define WATER_Q_HYDROGEN 0.4238

// Random numbers skipped between two velocity components
#define GMX_RANDOM_SKIP 5

/* Builds the water system of the GROMACS non-bonded benchmark (nbnxm/benchmark):
 * 1000 3-site water molecules (O,H,H) in a 3.107 nm box, replicated param->size times
 * (a power of two). Units are GROMACS' (nm, kJ/mol, e). Only oxygen has LJ parameters
 * and all atoms of a molecule are mutually excluded. */
static int buildGmxBenchAtoms(Atom* atom, Parameter* param)
{
    int factor = param->size;

    if (factor < 1 || (factor & (factor - 1)) != 0) {
        fprintf(stderr, "size factor must be a power of two\n");
        exit(-1);
    }

    int fac[3] = { 1, 1, 1 };
    for (int dim = 0; factor > 1;) {
        fac[dim] *= 2;
        factor >>= 1;
        dim = (dim + 1) % 3;
    }

    const int baseN = (int)(sizeof(coordinates1000) / sizeof(coordinates1000[0]));
    const int total = baseN * fac[0] * fac[1] * fac[2];
    while (atom->Nmax < total) {
        growAtom(atom);
    }

    atom->Natoms = atom->Nlocal = total;

    param->xlo = param->ylo = param->zlo = 0.0;
    param->xhi                           = fac[0] * box1000[0][0];
    param->yhi                           = fac[1] * box1000[1][1];
    param->zhi                           = fac[2] * box1000[2][2];
    param->xprd                          = param->xhi;
    param->yprd                          = param->yhi;
    param->zprd                          = param->zhi;

    const int numAtomsInMolecule = WATER_ATOMS;
    const int typeOxygen         = 0;
    const int typeHydrogen       = 1;
    const int ntypes             = WATER_NTYPES;
    atom->ntypes                 = ntypes;
    atom->molsize                = numAtomsInMolecule;

    int idx = 0;
    for (int x = 0; x < fac[0]; ++x) {
        for (int y = 0; y < fac[1]; ++y) {
            for (int z = 0; z < fac[2]; ++z) {
                const MD_FLOAT sx = (MD_FLOAT)x * (MD_FLOAT)box1000[0][0];
                const MD_FLOAT sy = (MD_FLOAT)y * (MD_FLOAT)box1000[1][1];
                const MD_FLOAT sz = (MD_FLOAT)z * (MD_FLOAT)box1000[2][2];

                for (int p = 0; p < baseN; ++p, ++idx) {
                    atom_x(idx) = (MD_FLOAT)coordinates1000[p][0] + sx;
                    atom_y(idx) = (MD_FLOAT)coordinates1000[p][1] + sy;
                    atom_z(idx) = (MD_FLOAT)coordinates1000[p][2] + sz;

                    int n = idx + 1;
                    for (int m = 0; m < GMX_RANDOM_SKIP; m++) {
                        myrandom(&n);
                    }
                    atom_vx(idx) = myrandom(&n);
                    for (int m = 0; m < GMX_RANDOM_SKIP; m++) {
                        myrandom(&n);
                    }
                    atom_vy(idx) = myrandom(&n);
                    for (int m = 0; m < GMX_RANDOM_SKIP; m++) {
                        myrandom(&n);
                    }
                    atom_vz(idx) = myrandom(&n);

                    atom->type[idx] = (p % numAtomsInMolecule == 0) ? typeOxygen
                                                                    : typeHydrogen;
                }
            }
        }
    }

    // Put all atoms in the box individually, as the GROMACS benchmark does. The domain
    // box (atom->mybox) is set up later, by setup() or the MPI grid
    const MD_FLOAT prd[3] = { param->xprd, param->yprd, param->zprd };
    for (int i = 0; i < total; i++) {
        MD_FLOAT* r[3] = { &atom_x(i), &atom_y(i), &atom_z(i) };
        for (int dim = 0; dim < 3; dim++) {
            if (*r[dim] < 0.0) *r[dim] += prd[dim];
            if (*r[dim] >= prd[dim]) *r[dim] -= prd[dim];
        }
    }

    // SPC/E-like parameters: LJ on oxygen only, given as C6/C12 in the GROMACS
    // benchmark. The benchmark hands them to the NBNxM atom data, which expects the
    // pre-multiplied 6*C6 and 12*C12, so the force it evaluates is the one of a LJ
    // potential with C6/6 and C12/12. Use those effective values to get the same forces.
    // Convert to the epsilon/sigma representation used by the kernels:
    // sigma^6 = C12 / C6, epsilon = C6^2 / (4 C12), i.e. sqrt(eps) = C6 / (2 sqrt(C12))
    const double c6Oxygen       = WATER_C6_OXYGEN / 6.0;
    const double c12Oxygen      = WATER_C12_OXYGEN / 12.0;
    const double chargeOxygen   = WATER_Q_OXYGEN;
    const double chargeHydrogen = WATER_Q_HYDROGEN;
    atom->sqrt_epsilon_per_type = allocate(ALIGNMENT, ntypes * sizeof(MD_FLOAT));
    atom->sigma3_per_type       = allocate(ALIGNMENT, ntypes * sizeof(MD_FLOAT));
    atom->q_per_type            = allocate(ALIGNMENT, ntypes * sizeof(MD_FLOAT));
    atom->sqrt_epsilon_per_type[typeOxygen]   = (MD_FLOAT)(c6Oxygen /
                                                           (2.0 * sqrt(c12Oxygen)));
    atom->sigma3_per_type[typeOxygen]         = (MD_FLOAT)sqrt(c12Oxygen / c6Oxygen);
    atom->sqrt_epsilon_per_type[typeHydrogen] = (MD_FLOAT)0.0;
    atom->sigma3_per_type[typeHydrogen]       = (MD_FLOAT)0.0;
    atom->q_per_type[typeOxygen]              = (MD_FLOAT)chargeOxygen;
    atom->q_per_type[typeHydrogen]            = (MD_FLOAT)chargeHydrogen;

    atom->epsilon    = allocate(ALIGNMENT, ntypes * ntypes * sizeof(MD_FLOAT));
    atom->sigma6     = allocate(ALIGNMENT, ntypes * ntypes * sizeof(MD_FLOAT));
    atom->cutforcesq = allocate(ALIGNMENT, ntypes * ntypes * sizeof(MD_FLOAT));
    atom->cutneighsq = allocate(ALIGNMENT, ntypes * ntypes * sizeof(MD_FLOAT));
    computeTypePairLJParameters(ntypes,
        atom->sqrt_epsilon_per_type,
        atom->sigma3_per_type,
        atom->epsilon,
        atom->sigma6);
    // One cut-off for all type pairs. Pairs with hydrogen have C6 = C12 = 0 and only
    // interact through their charges
    for (int k = 0; k < ntypes * ntypes; ++k) {
        atom->cutforcesq[k] = param->cutforce * param->cutforce;
        atom->cutneighsq[k] = param->cutneigh * param->cutneigh;
    }

    return total;
}

#endif
