/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <force.h>
#include <thermo.h>
#include <util.h>

#ifdef _MPI
#include <mpi.h>
static MPI_Datatype type_float = (sizeof(MD_FLOAT) == 4) ? MPI_FLOAT : MPI_DOUBLE;
#endif

static int* steparr;
static MD_FLOAT* tmparr;
static MD_FLOAT* engarr;
static MD_FLOAT* prsarr;
static MD_FLOAT mvv2e;
static MD_FLOAT dof_boltz;
static MD_FLOAT t_scale;
static MD_FLOAT p_scale;
static MD_FLOAT e_scale;
static MD_FLOAT t_act;
static MD_FLOAT p_act;
static MD_FLOAT e_act;
static int mstat;

/* exported subroutines */
void setupThermo(Parameter* param, int natoms)
{
    int maxstat = param->ntimes / param->nstat + 2;

    steparr = (int*)malloc(maxstat * sizeof(int));
    tmparr  = (MD_FLOAT*)malloc(maxstat * sizeof(MD_FLOAT));
    engarr  = (MD_FLOAT*)malloc(maxstat * sizeof(MD_FLOAT));
    prsarr  = (MD_FLOAT*)malloc(maxstat * sizeof(MD_FLOAT));

    if (param->force_field == FF_LJ) {
        mvv2e     = 1.0;
        dof_boltz = (natoms * 3 - 3);
        t_scale   = mvv2e / dof_boltz;
        p_scale   = 1.0 / 3 / param->xprd / param->yprd / param->zprd;
        e_scale   = 0.5;
    } else if (param->force_field == FF_EAM) {
        mvv2e     = 1.036427e-04;
        dof_boltz = (natoms * 3 - 3) * 8.617343e-05;
        t_scale   = mvv2e / dof_boltz;
        p_scale   = 1.602176e+06 / 3 / param->xprd / param->yprd / param->zprd;
        e_scale   = 524287.985533; // 16.0;
        param->dtforce /= mvv2e;
    }
}

void computeThermo(int iflag, Parameter* param, Atom* atom)
{
    MD_FLOAT t     = 0.0, p;
    MD_FLOAT t_sum = 0.0;
    int me         = 0;

#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif

    for (int i = 0; i < atom->Nlocal; i++) {
        t += (atom_vx(i) * atom_vx(i) + atom_vy(i) * atom_vy(i) +
                 atom_vz(i) * atom_vz(i)) *
             param->mass;
    }
#ifdef _MPI
    MPI_Reduce(&t, &t_sum, 1, type_float, MPI_SUM, 0, MPI_COMM_WORLD);
    t = t_sum;
#endif

    if (me == 0) {
        t         = t * t_scale;
        p         = (t * dof_boltz) * p_scale;
        int istep = iflag;

        if (iflag == -1) {
            istep = param->ntimes;
        }
        if (iflag == 0) {
            mstat = 0;
        }

        steparr[mstat] = istep;
        tmparr[mstat]  = t;
        prsarr[mstat]  = p;
        mstat++;
        fprintf(stdout, "%i\t%e\t%e\n", istep, t, p);
        fflush(stdout);
    }
}

void adjustThermo(Parameter* param, Atom* atom)
{
    /* zero center-of-mass motion */
    MD_FLOAT vxtot = 0.0;
    MD_FLOAT vytot = 0.0;
    MD_FLOAT vztot = 0.0;
    MD_FLOAT vtot[3];

    for (int i = 0; i < atom->Nlocal; i++) {
        vxtot += atom_vx(i);
        vytot += atom_vy(i);
        vztot += atom_vz(i);
    }

    vtot[0] = vxtot;
    vtot[1] = vytot;
    vtot[2] = vztot;

#ifdef _MPI
    MD_FLOAT v_sum[3];
    MPI_Allreduce(vtot, v_sum, 3, type_float, MPI_SUM, MPI_COMM_WORLD);
    vtot[0] = v_sum[0];
    vtot[1] = v_sum[1];
    vtot[2] = v_sum[2];
#endif
    vxtot = vtot[0] / atom->Natoms;
    vytot = vtot[1] / atom->Natoms;
    vztot = vtot[2] / atom->Natoms;

    for (int i = 0; i < atom->Nlocal; i++) {
        atom_vx(i) -= vxtot;
        atom_vy(i) -= vytot;
        atom_vz(i) -= vztot;
    }

    // t_act      = 0;
    MD_FLOAT t     = 0.0;
    MD_FLOAT t_sum = 0.0;

    for (int i = 0; i < atom->Nlocal; i++) {
        t += (atom_vx(i) * atom_vx(i) + atom_vy(i) * atom_vy(i) +
                 atom_vz(i) * atom_vz(i)) *
             param->mass;
    }
#ifdef _MPI
    MPI_Allreduce(&t, &t_sum, 1, type_float, MPI_SUM, MPI_COMM_WORLD);
    t = t_sum;
#endif

    t *= t_scale;
    MD_FLOAT factor = sqrt(param->temp / t);

    for (int i = 0; i < atom->Nlocal; i++) {
        atom_vx(i) *= factor;
        atom_vy(i) *= factor;
        atom_vz(i) *= factor;
    }
}
