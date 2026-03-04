/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <atom.h>
#include <force.h>
#include <parameter.h>
#include <util.h>
#if defined(CLUSTERPAIR) || !defined(USE_REFERENCE_KERNEL)
#include <simd.h>
#endif

void initParameter(Parameter* param) {
    param->input_file      = NULL;
    param->vtk_file        = NULL;
    param->xtc_file        = NULL;
    param->eam_file        = NULL;
    param->write_atom_file = NULL;
    param->force_field     = FF_LJ;
    param->epsilon         = 1.0;
    param->sigma           = 1.0;
    param->sigma6          = 1.0;
    param->rho             = 0.8442;
#ifdef ONE_ATOM_TYPE
    param->ntypes        = 1;
#else
    param->ntypes        = 4;
#endif
    param->ntimes        = 200;
    param->dt            = 0.005;
    param->nx            = 32;
    param->ny            = 32;
    param->nz            = 32;
    param->pbc_x         = 1;
    param->pbc_y         = 1;
    param->pbc_z         = 1;
    param->cutforce      = 2.5;
    param->skin          = 0.3;
    param->cutneigh      = param->cutforce + param->skin;
    param->temp          = 1.44;
    param->nstat         = 100;
    param->mass          = 1.0;
    param->dtforce       = 0.5 * param->dt;
    param->reneigh_every = 20;
    param->resort_every  = 400;
    param->prune_every   = 1000;
    param->x_out_every   = 20;
    param->v_out_every   = 5;
    param->half_neigh    = 0;
    param->proc_freq     = 2.4;
#ifdef CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
    param->super_clustering = 1;
#else
    param->super_clustering = 0;
#endif
    // MPI
    param->balance       = 0;
    param->method        = 0;
    param->balance_every = param->reneigh_every;
    param->setup         = 1;
}

void readParameter(Parameter* param, const char* filename) {
    FILE* fp = fopen(filename, "r");
    char line[MAXLINE];
    int i;

    if (!fp) {
        fprintf(stderr, "Could not open parameter file: %s\n", filename);
        exit(-1);
    }

    while (fgets(line, MAXLINE, fp) != NULL) {
        for (i = 0; line[i] != '\0' && line[i] != '#'; i++)
            ;
        line[i] = '\0';

        char* tok = strtok(line, " ");
        char* val = strtok(NULL, " ");

#define PARSE_PARAM(p, f)                                                                \
    if (strncmp(tok, #p, sizeof(#p) / sizeof(#p[0]) - 1) == 0) {                         \
        param->p = f(val);                                                               \
    }
#define PARSE_STRING(p) PARSE_PARAM(p, strdup)
#define PARSE_INT(p)    PARSE_PARAM(p, atoi)
#define PARSE_REAL(p)   PARSE_PARAM(p, atof)

        if (tok != NULL && val != NULL) {
            PARSE_PARAM(force_field, str2ff);
            PARSE_STRING(input_file);
            PARSE_STRING(eam_file);
            PARSE_STRING(vtk_file);
            PARSE_STRING(xtc_file);
            PARSE_REAL(epsilon);
            PARSE_REAL(sigma);
            PARSE_REAL(rho);
            PARSE_REAL(dt);
            PARSE_REAL(cutforce);
            PARSE_REAL(skin);
            PARSE_REAL(temp);
            PARSE_REAL(mass);
            PARSE_REAL(proc_freq);
            PARSE_INT(ntypes);
            PARSE_INT(ntimes);
            PARSE_INT(nx);
            PARSE_INT(ny);
            PARSE_INT(nz);
            PARSE_INT(pbc_x);
            PARSE_INT(pbc_y);
            PARSE_INT(pbc_z);
            PARSE_INT(nstat);
            PARSE_INT(reneigh_every);
            PARSE_INT(resort_every);
            PARSE_INT(prune_every);
            PARSE_INT(x_out_every);
            PARSE_INT(v_out_every);
            PARSE_INT(half_neigh);
            PARSE_INT(method);
            PARSE_INT(balance);
            PARSE_INT(balance_every);
            PARSE_INT(super_clustering);
        }
    }

    // Update dtforce
    param->dtforce = 0.5 * param->dt;

    // Update sigma6 parameter
    MD_FLOAT s2   = param->sigma * param->sigma;
    param->sigma6 = s2 * s2 * s2;

    // Update balance parameter, 10 could be change
    param->balance_every *= param->reneigh_every;
    fclose(fp);
}

void printParameter(Parameter* param) {
    fprintf(stdout, "Parameters:\n");
    if (param->input_file != NULL) {
        fprintf(stdout, "\tInput file: %s\n", param->input_file);
    }

    if (param->vtk_file != NULL) {
        fprintf(stdout, "\tVTK file: %s\n", param->vtk_file);
    }

    if (param->xtc_file != NULL) {
        fprintf(stdout, "\tXTC file: %s\n", param->xtc_file);
    }

    if (param->eam_file != NULL) {
        fprintf(stdout, "\tEAM file: %s\n", param->eam_file);
    }

    fprintf(stdout, "\tForce field: %s\n", ff2str(param->force_field));
#ifdef CLUSTER_M
    fprintf(stdout,
        "\tKernel: %s, MxN: %dx%d, Vector width: %d\n",
        KERNEL_NAME,
        CLUSTER_M,
        CLUSTER_N,
        VECTOR_WIDTH);
#else
    fprintf(stdout, "\tKernel: %s\n", KERNEL_NAME);
#endif

#ifdef CUDA_TARGET
    fprintf(stdout, "\tSIMD Intrinsics: CUDA\n");
    fprintf(stdout, "\tSuper-clustering: %s\n", (param->super_clustering) ? "yes" : "no");
#else
    fprintf(stdout, "\tSIMD Intrinsics: %s\n", SIMD_INTRINSICS);
#endif
    fprintf(stdout, "\tAtom data layout: %s\n", POS_DATA_LAYOUT);
    fprintf(stdout, "\tNeighbor-lists data layout: %s\n", NBLIST_DATA_LAYOUT);
    fprintf(stdout, "\tFloating-point precision: %s\n", PRECISION_STRING);
    fprintf(stdout,
        "\tUnit cells (nx, ny, nz): %d, %d, %d\n",
        param->nx,
        param->ny,
        param->nz);
    fprintf(stdout,
        "\tDomain box sizes (x, y, z): %e, %e, %e\n",
        param->xprd,
        param->yprd,
        param->zprd);
    fprintf(stdout,
        "\tPeriodic (x, y, z): %d, %d, %d\n",
        param->pbc_x,
        param->pbc_y,
        param->pbc_z);
    fprintf(stdout, "\tLattice size: %e\n", param->lattice);
    fprintf(stdout, "\tEpsilon: %e\n", param->epsilon);
    fprintf(stdout, "\tSigma: %e\n", param->sigma);
    fprintf(stdout, "\tTemperature: %e\n", param->temp);
    fprintf(stdout, "\tRHO: %e\n", param->rho);
    fprintf(stdout, "\tMass: %e\n", param->mass);
    fprintf(stdout, "\tNumber of types: %d\n", param->ntypes);
    fprintf(stdout, "\tNumber of timesteps: %d\n", param->ntimes);
    fprintf(stdout, "\tReport stats every (timesteps): %d\n", param->nstat);
    fprintf(stdout, "\tReneighbor every (timesteps): %d\n", param->reneigh_every);
#ifdef SORT_ATOMS
    fprintf(stdout, "\tResort atoms every (timesteps): %d\n", param->resort_every);
#else
    fprintf(stdout, "\tSort atoms: no\n");
#endif
#ifdef ONE_ATOM_TYPE
    fprintf(stdout, "\tSingle atom type: true\n");
#else
    fprintf(stdout, "\tSingle atom type: false\n");
#endif
    fprintf(stdout, "\tPrune every (timesteps): %d\n", param->prune_every);
    fprintf(stdout, "\tOutput positions every (timesteps): %d\n", param->x_out_every);
    fprintf(stdout, "\tOutput velocities every (timesteps): %d\n", param->v_out_every);
    fprintf(stdout, "\tDelta time (dt): %e\n", param->dt);
    fprintf(stdout, "\tCutoff radius: %e\n", param->cutforce);
    fprintf(stdout, "\tSkin: %e\n", param->skin);
    fprintf(stdout, "\tHalf neighbor lists: %d\n", param->half_neigh);
    fprintf(stdout, "\tProcessor frequency (GHz): %.4f\n", param->proc_freq);

    // ================ New MPI features =============
#ifdef _MPI
    char str[20];
    strcpy(str,
        (param->method == 1)   ? "Half Shell"
        : (param->method == 2) ? "Eight Shell"
        : (param->method == 3) ? "Half Stencil"
                               : "Full Shell");
    fprintf(stdout, "\tMethod: %s\n", str);
    strcpy(str,
        (param->balance == 1)   ? "mean RCB"
        : (param->balance == 2) ? "mean Time RCB"
        : (param->balance == 3) ? "Staggered"
                                : "cartesian");
    fprintf(stdout, "\tPartition: %s\n", str);
    if (param->balance)
        fprintf(stdout, "\tRebalancing every (timesteps): %d\n", param->balance_every);
#endif
    fflush(stdout);
}
