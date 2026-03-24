/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
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
// Fallback for SIMD_INTRINSICS when simd.h not included
#ifndef SIMD_INTRINSICS
#define SIMD_INTRINSICS "scalar"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _MPI
#include <mpi.h>
#endif

// Compute per-type sqrt(epsilon) and sigma^3 for geometric combination rule
void computePerTypeLJParameters(
    int ntypes, Parameter* param, MD_FLOAT* sqrt_epsilon, MD_FLOAT* sigma3)
{
    int param_ntypes = param->ntypes;
    for (int t = 0; t < ntypes; t++) {
        MD_FLOAT eps    = (t < param_ntypes && param->epsilon_per_type)
                              ? param->epsilon_per_type[t]
                              : param->epsilon;
        MD_FLOAT sig    = (t < param_ntypes && param->sigma_per_type)
                              ? param->sigma_per_type[t]
                              : param->sigma;
        sqrt_epsilon[t] = sqrt(eps);
        sigma3[t]       = sig * sig * sig;
    }
}

// Compute type-pair matrix from per-type arrays (for scalar kernels)
// epsilon[i*ntypes+j] = sqrt_epsilon[i] * sqrt_epsilon[j]
// sigma6[i*ntypes+j] = sigma3[i] * sigma3[j]
void computeTypePairLJParameters(int ntypes,
    MD_FLOAT* sqrt_epsilon_per_type,
    MD_FLOAT* sigma3_per_type,
    MD_FLOAT* epsilon,
    MD_FLOAT* sigma6)
{
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    // For single type, we only need one value (will be broadcast anyway)
    epsilon[0] = sqrt_epsilon_per_type[0] * sqrt_epsilon_per_type[0];
    sigma6[0]  = sigma3_per_type[0] * sigma3_per_type[0];
#else
    // Geometric combination: compute all pairs
    for (int i = 0; i < ntypes; i++) {
        for (int j = 0; j < ntypes; j++) {
            int idx      = i * ntypes + j;
            epsilon[idx] = sqrt_epsilon_per_type[i] * sqrt_epsilon_per_type[j];
            sigma6[idx]  = sigma3_per_type[i] * sigma3_per_type[j];
        }
    }
#endif
}

void initParameter(Parameter* param)
{
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
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    param->ntypes = 1;
#else
    param->ntypes = 4;
#endif
    param->epsilon_per_type = NULL;
    param->sigma_per_type   = NULL;
    param->ntimes           = 200;
    param->dt               = 0.005;
    param->nx               = 32;
    param->ny               = 32;
    param->nz               = 32;
    param->pbc_x            = 1;
    param->pbc_y            = 1;
    param->pbc_z            = 1;
    param->cutforce         = 2.5;
    param->skin             = 0.3;
    param->cutneigh         = param->cutforce + param->skin;
    param->temp             = 1.44;
    param->nstat            = 100;
    param->mass             = 1.0;
    param->dtforce          = 0.5 * param->dt;
    param->reneigh_every    = 20;
    param->resort_every     = 400;
    param->prune_every      = 1000;
    param->x_out_every      = 20;
    param->v_out_every      = 5;
    param->half_neigh       = 0;
    param->proc_freq        = 2.4;
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

void readParameter(Parameter* param, const char* filename)
{
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

            // Parse per-type epsilon
            if (strncmp(tok, "epsilon_type_", 13) == 0) {
                int type_id = atoi(tok + 13);
                if (type_id >= 0 && type_id < param->ntypes) {
                    if (param->epsilon_per_type == NULL) {
                        param->epsilon_per_type = (MD_FLOAT*)malloc(
                            param->ntypes * sizeof(MD_FLOAT));
                        for (int i = 0; i < param->ntypes; i++) {
                            param->epsilon_per_type[i] = param->epsilon;
                        }
                    }
                    param->epsilon_per_type[type_id] = atof(val);
                }
            }

            // Parse per-type sigma
            if (strncmp(tok, "sigma_type_", 11) == 0) {
                int type_id = atoi(tok + 11);
                if (type_id >= 0 && type_id < param->ntypes) {
                    if (param->sigma_per_type == NULL) {
                        param->sigma_per_type = (MD_FLOAT*)malloc(
                            param->ntypes * sizeof(MD_FLOAT));
                        for (int i = 0; i < param->ntypes; i++) {
                            param->sigma_per_type[i] = param->sigma;
                        }
                    }
                    param->sigma_per_type[type_id] = atof(val);
                }
            }
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

void printParameter(Parameter* param)
{
    fprintf(stdout, "SIMULATION PARAMETERS\n");
    fprintf(stdout,
        "-------------------------------------------------------------------------------"
        "\n");

    // Computational kernel
    fprintf(stdout, "  Computational Kernel:\n");
    fprintf(stdout,
        "    Force field:                       %s\n",
        ff2str(param->force_field));
#ifdef CLUSTER_M
    fprintf(stdout,
        "    Kernel:                            %s (MxN: %dx%d, Vector width: %d)\n",
        KERNEL_NAME,
        CLUSTER_M,
        CLUSTER_N,
        VECTOR_WIDTH);
#else
    fprintf(stdout, "    Kernel:                            %s\n", KERNEL_NAME);
#endif

#ifdef CUDA_TARGET
    fprintf(stdout, "    SIMD/Architecture:                 CUDA\n");
    fprintf(stdout,
        "    Super-clustering:                  %s\n",
        (param->super_clustering) ? "yes" : "no");
#else
    fprintf(stdout, "    SIMD/Architecture:                 %s\n", SIMD_INTRINSICS);
#endif
    fprintf(stdout, "    Atom data layout:                  %s\n", POS_DATA_LAYOUT);
    fprintf(stdout, "    Neighbor-list layout:              %s\n", NBLIST_DATA_LAYOUT);
    fprintf(stdout, "    FP precision:                      %s\n", PRECISION_STRING);

    // System configuration
    fprintf(stdout, "\n  System Configuration:\n");
    if (param->input_file != NULL) {
        fprintf(stdout, "    Input file:                        %s\n", param->input_file);
    }
    if (param->vtk_file != NULL) {
        fprintf(stdout, "    VTK file:                          %s\n", param->vtk_file);
    }
    if (param->xtc_file != NULL) {
        fprintf(stdout, "    XTC file:                          %s\n", param->xtc_file);
    }
    if (param->eam_file != NULL) {
        fprintf(stdout, "    EAM file:                          %s\n", param->eam_file);
    }
    fprintf(stdout,
        "    Unit cells (nx,ny,nz):             %d x %d x %d\n",
        param->nx,
        param->ny,
        param->nz);
    fprintf(stdout,
        "    Domain box sizes:                  %.2e x %.2e x %.2e\n",
        param->xprd,
        param->yprd,
        param->zprd);
    fprintf(stdout,
        "    Periodic boundary:                 %s %s %s\n",
        param->pbc_x ? "x" : "-",
        param->pbc_y ? "y" : "-",
        param->pbc_z ? "z" : "-");

    // Physical parameters
    fprintf(stdout, "\n  Physical Parameters:\n");
    fprintf(stdout, "    Lattice constant:                  %.6e\n", param->lattice);
    fprintf(stdout, "    Temperature:                       %.6e\n", param->temp);
    fprintf(stdout, "    Density:                           %.6e\n", param->rho);
    fprintf(stdout, "    Mass:                              %.6e\n", param->mass);
    fprintf(stdout, "    Epsilon:                           %.6e\n", param->epsilon);
    fprintf(stdout, "    Sigma:                             %.6e\n", param->sigma);
    fprintf(stdout, "    Number of types:                   %d\n", param->ntypes);
    fprintf(stdout, "    LJ combination rule:               %s\n", LJ_COMB_RULE_NAME);
#if LJ_COMB_RULE == LJ_COMB_GEOM
    if (param->epsilon_per_type != NULL) {
        fprintf(stdout, "    Per-type epsilon:                  ");
        for (int i = 0; i < param->ntypes; i++) {
            fprintf(stdout, "[%d]=%.4e ", i, param->epsilon_per_type[i]);
        }
        fprintf(stdout, "\n");
    }
    if (param->sigma_per_type != NULL) {
        fprintf(stdout, "    Per-type sigma:                    ");
        for (int i = 0; i < param->ntypes; i++) {
            fprintf(stdout, "[%d]=%.4e ", i, param->sigma_per_type[i]);
        }
        fprintf(stdout, "\n");
    }
#endif

    // Simulation parameters
    fprintf(stdout, "\n  Simulation Control:\n");
    fprintf(stdout, "    Timesteps:                         %d\n", param->ntimes);
    fprintf(stdout, "    Timestep (dt):                     %.6e\n", param->dt);
    fprintf(stdout, "    Cutoff radius:                     %.6e\n", param->cutforce);
    fprintf(stdout, "    Skin distance:                     %.6e\n", param->skin);
    fprintf(stdout,
        "    Half neighbor-lists:               %s\n",
        param->half_neigh ? "yes" : "no");
    fprintf(stdout,
        "    Reneighbor every:                  %d steps\n",
        param->reneigh_every);
    fprintf(stdout, "    Report stats every:                %d steps\n", param->nstat);
#ifdef SORT_ATOMS
    fprintf(stdout,
        "    Resort atoms every:                %d steps\n",
        param->resort_every);
#endif
#if LJ_COMB_RULE == LJ_COMB_SINGLE
    fprintf(stdout, "    LJ combination rule:               single\n");
#elif LJ_COMB_RULE == LJ_COMB_GEOM
    fprintf(stdout, "    LJ combination rule:               geometric\n");
#else
    fprintf(stdout, "    LJ combination rule:               none\n");
#endif
    fprintf(stdout,
        "    Prune every:                       %d steps\n",
        param->prune_every);
    fprintf(stdout,
        "    Output positions:                  every %d steps\n",
        param->x_out_every);
    fprintf(stdout,
        "    Output velocities:                 every %d steps\n",
        param->v_out_every);
    fprintf(stdout,
        "    Processor freq:                    %.2f GHz\n",
        param->proc_freq);

    // Parallel configuration
    fprintf(stdout, "\n  Parallel Configuration:\n");
#ifdef _MPI
    int nranks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    fprintf(stdout, "    MPI ranks:                         %d\n", nranks);
    char str[20];
    strcpy(str,
        (param->method == 1)   ? "Half Shell"
        : (param->method == 2) ? "Eight Shell"
        : (param->method == 3) ? "Half Stencil"
                               : "Full Shell");
    fprintf(stdout, "    MPI method:                        %s\n", str);
    strcpy(str,
        (param->balance == 1)   ? "mean RCB"
        : (param->balance == 2) ? "mean Time RCB"
        : (param->balance == 3) ? "Staggered"
                                : "cartesian");
    fprintf(stdout, "    Domain partition:                  %s\n", str);
    if (param->balance)
        fprintf(stdout,
            "    Rebalancing every:                 %d steps\n",
            param->balance_every);
#else
    fprintf(stdout, "    MPI ranks:                         1 (not compiled)\n");
#endif

#ifdef _OPENMP
    int nthreads  = 0;
    int chunkSize = 0;
    omp_sched_t schedKind;
    char schedType[10];
#pragma omp parallel
#pragma omp master
    {
        omp_get_schedule(&schedKind, &chunkSize);

        switch (schedKind) {
        case omp_sched_static:
            strcpy(schedType, "static");
            break;
        case omp_sched_dynamic:
            strcpy(schedType, "dynamic");
            break;
        case omp_sched_guided:
            strcpy(schedType, "guided");
            break;
        case omp_sched_auto:
            strcpy(schedType, "auto");
            break;
        case omp_sched_monotonic:
            strcpy(schedType, "auto");
            break;
        }

        nthreads = omp_get_max_threads();
    }

    fprintf(stdout, "    OpenMP threads:                    %d\n", nthreads);
    fprintf(stdout,
        "    OpenMP schedule:                   (%s,%d)\n",
        schedType,
        chunkSize);
#else
    fprintf(stdout, "    OpenMP threads:                    1 (not compiled)\n");
#endif

    fprintf(stdout,
        "-------------------------------------------------------------------------------"
        "\n");
    fflush(stdout);
}
