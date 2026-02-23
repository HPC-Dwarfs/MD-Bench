/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdlib.h>
//---
#include <allocate.h>
#include <atom.h>
#include <xtc.h>

#ifdef XTC_OUTPUT
#include <gromacs/fileio/xtcio.h>

static struct t_fileio* xtc_file = NULL;
static rvec* x_buf               = NULL;
static rvec basis[3];

void xtc_init(const char* filename, Atom* atom, int timestep)
{
    basis[0][XX] = 1.0;
    basis[0][YY] = 0.0;
    basis[0][ZZ] = 0.0;
    basis[1][XX] = 0.0;
    basis[1][YY] = 1.0;
    basis[1][ZZ] = 0.0;
    basis[2][XX] = 0.0;
    basis[2][YY] = 0.0;
    basis[2][ZZ] = 1.0;

    xtc_file = open_xtc(filename, "w");
    x_buf    = (rvec*)allocate(ALIGNMENT, sizeof(rvec) * (atom->Nlocal + 1));
    xtc_write(atom, timestep, 1, 1);
}

void xtc_write(Atom* atom, int timestep, int write_pos, int write_vel)
{
    int i = 0;
    for (int ci = 0; ci < atom->Nclusters_local; ++ci) {
        int ci_vec_base = CI_VECTOR_BASE_INDEX(ci);
        MD_FLOAT* ci_x  = &atom->cl_x[ci_vec_base];
        for (int cii = 0; cii < atom->clusters[ci].natoms; ++cii) {
            x_buf[i][XX] = ci_x[CL_X_INDEX(cii)];
            x_buf[i][YY] = ci_x[CL_Y_INDEX(cii)];
            x_buf[i][ZZ] = ci_x[CL_Z_INDEX(cii)];
            i++;
        }
    }

    write_xtc(xtc_file,
        atom->Nlocal,
        timestep,
        0.0,
        (const rvec*)basis,
        (const rvec*)x_buf,
        1000);
}

void xtc_end()
{
    free(x_buf);
    close_xtc(xtc_file);
}
#endif
