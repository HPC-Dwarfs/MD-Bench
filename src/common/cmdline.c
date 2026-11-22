/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmdline.h>
#include <comm.h>
#include <force.h>
#include <parameter.h>
#include <util.h>

#define HLINE                                                                            \
    "-------------------------------------------------------------------------------\n"

#ifdef CLUSTER_PAIR
#define SCHEME_NAME "cluster-pair"
#define MAX_METHOD  2
#else
#define SCHEME_NAME "verlet-list"
#define MAX_METHOD  3
#endif

static void cmdError(Comm* comm, const char* format, ...)
{
    if (comm->myproc == 0) {
        va_list args;
        va_start(args, format);
        fprintf(stderr, "Error: ");
        vfprintf(stderr, format, args);
        fprintf(stderr, "\n");
        va_end(args);
    }
    endComm(comm);
    exit(EXIT_FAILURE);
}

static void printHelp(void)
{
    printf("MD-Bench: A performance-oriented prototyping harness for MD algorithms "
           "(" SCHEME_NAME " scheme)\n");
    printf(HLINE);
    printf("Input and system:\n");
    printf("  -p / --param <file>       read parameters from a file (can be given more "
           "than once)\n");
    printf("  -f <string>               force field: lj (default), eam or lj_table\n");
    printf("  -i <file>                 input file with atom positions (dump/pdb/gro)\n");
    printf("  -e <file>                 EAM potential file\n");
    printf("  -nx/-ny/-nz <int>         unit cells of the lattice in x/y/z\n");
    printf("  -m / --mass <real>        atom mass\n");
    printf("  --ntypes <int>            number of atom types\n");
    printf("  -setup <int>              1 (default) creates the FCC lattice of the "
           "argon/copper setup, 0 skips it\n");
    printf("  -gmxbench                 water system of the GROMACS non-bonded benchmark "
           "(cut-off 1.0, no skin unless -r/-s are given)\n");
    printf("  -size <int>               replicate the -gmxbench system this many times "
           "(power of two)\n");
    printf("  -w <file>                 write the input atoms to a file (verlet-list)\n");
    printf("Simulation and neighbor lists:\n");
    printf("  -n / --nsteps <int>       number of time steps\n");
    printf("  -r / --radius <real>      cut-off radius\n");
    printf("  -s / --skin <real>        skin (Verlet buffer of the force list)\n");
    printf("  --outer-skin <real>       additive outer skin, > 0 enables the dual pair "
           "list (outer cut-off = cut-off + skin + outer skin)\n");
    printf("  --prune-every <int>       steps between prunings of the inner list\n");
    printf("  --reneigh-every <int>     steps between neighbor-list builds\n");
    printf("  --displacement-reneigh    rebuild when an atom moved more than skin/2 "
           "(--reneigh-every 0 disables the fixed interval)\n");
    printf("  --sort-every <int>        sort the atoms spatially every <int> steps, 0 "
           "disables it (verlet-list; default: off on the CPU, at every neighbor-list "
           "build on the GPU without MPI)\n");
    printf("  -half <int>               half (1) or full (0) neighbor lists\n");
    printf("  -method <int>             ghost interactions with MPI: 0 full shell, 1 "
           "half shell, 2 eighth shell, 3 half stencil (verlet-list)\n");
    printf("  -bal <int>                load balancing with MPI: 0 none, 1 RCB, 2 mean "
           "time RCB, 3 staggered\n");
    printf("Force kernel:\n");
    printf(
        "  --lj-comb <string>        LJ combination rule: full (default, type-pair "
        "tables), geometric (per-type parameters), single (one set of parameters, "
        "all types alike) or auto (single when all types are alike, else geometric)\n");
    printf("  -energy                   also compute the LJ and Coulomb energies "
           "(reported by --kernel-bench and --dump-forces)\n");
    printf("  -halflj                   LJ only for the first half of the i-clusters "
           "with at most half of their atoms having LJ (CPU cluster-pair)\n");
    printf("  --no-accurate-recip       use the 1/r and 1/r^2 hardware estimates of the "
           "SIMD kernels without the Newton-Raphson refinement, and the fast "
           "approximate division and rsqrt of the GPU kernels (as -use_fast_math) "
           "instead of the IEEE ones (default: accurate, as GROMACS; --accurate-recip "
           "selects it explicitly)\n");
    printf("  --no-prefetch             no j-cluster software prefetch (CPU "
           "cluster-pair SIMD kernels, GROMACS does not prefetch)\n");
    printf("  -ref                      GPU builds: the reference force kernels (one "
           "thread per atom pair, no packing, slicing or warp reductions, all forces "
           "with atomicAdd) instead of the optimized ones\n");
    printf("  --kernel-bench <int>      after the setup, only run the force kernel this "
           "many times on the fixed neighbor lists\n");
    printf("  --dump-forces <file>      write the forces of a single kernel call (atom "
           "id, fx, fy, fz) to a text file\n");
    printf("Output:\n");
    printf("  --freq <real>             processor frequency (GHz) for the statistics\n");
    printf("  --vtk <file>              VTK output for visualization\n");
    printf("  --xtc <file>              XTC trajectory (cluster-pair, build with "
           "XTC_OUTPUT=true)\n");
    printf("  -v / --verbose            per-routine timing breakdown in the report\n");
    printf("  -h / --help               this help\n");
    printf(HLINE);
}

void parseCommandLine(int argc, char** argv, Parameter* param, Comm* comm)
{
    int cutoffGiven = 0, skinGiven = 0;

// Value of the option argv[i], fails when it is missing
#define NEXT_ARG()                                                                       \
    ((i + 1 < argc) ? argv[++i] : (cmdError(comm, "%s needs a value", argv[i]), NULL))
#define IS(opt)       (strcmp(argv[i], opt) == 0)
#define IS2(opt, alt) (IS(opt) || IS(alt))

    for (int i = 1; i < argc; i++) {
        if (IS("-p") || IS2("--param", "--params")) {
            readParameter(param, NEXT_ARG());
        } else if (IS("-f")) {
            const char* ff = NEXT_ARG();
            if ((param->force_field = str2ff(ff)) < 0) {
                cmdError(comm, "invalid force field %s", ff);
            }
        } else if (IS("-i")) {
            param->input_file = strdup(NEXT_ARG());
        } else if (IS("-e")) {
            param->eam_file = strdup(NEXT_ARG());
        } else if (IS("-nx")) {
            param->nx = atoi(NEXT_ARG());
        } else if (IS("-ny")) {
            param->ny = atoi(NEXT_ARG());
        } else if (IS("-nz")) {
            param->nz = atoi(NEXT_ARG());
        } else if (IS2("-m", "--mass")) {
            param->mass = atof(NEXT_ARG());
        } else if (IS("--ntypes")) {
            param->ntypes = atoi(NEXT_ARG());
        } else if (IS("-setup")) {
            param->setup = atoi(NEXT_ARG());
        } else if (IS("-gmxbench")) {
            param->gmxbenchmark = true;
            if (param->input_file != NULL) {
                free(param->input_file);
                param->input_file = NULL;
            }
        } else if (IS("-size")) {
            param->size = atoi(NEXT_ARG());
        } else if (IS("-w")) {
            param->write_atom_file = strdup(NEXT_ARG());
        } else if (IS2("-n", "--nsteps")) {
            param->ntimes = atoi(NEXT_ARG());
        } else if (IS2("-r", "--radius")) {
            param->cutforce = atof(NEXT_ARG());
            cutoffGiven     = 1;
        } else if (IS2("-s", "--skin")) {
            param->skin = atof(NEXT_ARG());
            skinGiven   = 1;
        } else if (IS("--outer-skin")) {
            param->outer_skin = atof(NEXT_ARG());
        } else if (IS("--prune-every")) {
            param->prune_every = atoi(NEXT_ARG());
        } else if (IS("--reneigh-every")) {
            param->reneigh_every = atoi(NEXT_ARG());
        } else if (IS("--displacement-reneigh")) {
            param->displacement_reneigh = 1;
        } else if (IS("--sort-every")) {
            param->resort_every = atoi(NEXT_ARG());
        } else if (IS("-half")) {
            param->half_neigh = atoi(NEXT_ARG());
        } else if (IS("-method")) {
            param->method = atoi(NEXT_ARG());
        } else if (IS("-bal")) {
            param->balance = atoi(NEXT_ARG());
        } else if (IS("--lj-comb")) {
            const char* rule = NEXT_ARG();
            if ((param->lj_comb_rule = str2ljcomb(rule)) < 0) {
                cmdError(comm,
                    "invalid --lj-comb %s, use single, geometric, full or auto",
                    rule);
            }
        } else if (IS("-energy")) {
            param->compute_energy = true;
        } else if (IS("-halflj")) {
            param->half_lj = true;
        } else if (IS("--accurate-recip")) {
            param->accurate_reciprocal = 1;
        } else if (IS("--no-accurate-recip")) {
            param->accurate_reciprocal = 0;
        } else if (IS("--no-prefetch")) {
            param->kernel_prefetch = 0;
        } else if (IS("-ref")) {
            param->ref_kernel = 1;
        } else if (IS("--kernel-bench")) {
            param->kernel_bench = atoi(NEXT_ARG());
        } else if (IS("--dump-forces")) {
            param->dump_forces_file = strdup(NEXT_ARG());
        } else if (IS("--freq")) {
            param->proc_freq = atof(NEXT_ARG());
        } else if (IS("--vtk")) {
            param->vtk_file = strdup(NEXT_ARG());
        } else if (IS("--xtc")) {
            param->xtc_file = strdup(NEXT_ARG());
        } else if (IS2("-v", "--verbose")) {
            param->verbose = 1;
        } else if (IS2("-h", "--help")) {
            if (comm->myproc == 0) {
                printHelp();
            }
            endComm(comm);
            exit(EXIT_SUCCESS);
        } else {
            cmdError(comm, "unknown option %s, see -h", argv[i]);
        }
    }
#undef NEXT_ARG
#undef IS
#undef IS2

    if (param->method < 0 || param->method > MAX_METHOD) {
        cmdError(comm,
            "-method %d does not exist for the " SCHEME_NAME " scheme",
            param->method);
    }
    if (param->balance < 0 || param->balance > 3) {
        cmdError(comm, "-bal %d does not exist", param->balance);
    }
    if (param->balance > 0 && param->method == 1) {
        cmdError(comm, "the half shell method (-method 1) does not support -bal");
    }
    if (param->gmxbenchmark && param->kernel_bench == 0 && param->ntimes != 0) {
        // The water has no bonds or constraints and is in GROMACS units: no MD
        cmdError(comm,
            "-gmxbench is a force kernel benchmark, use --kernel-bench <n> or -n 0");
    }
    // The molecules of -gmxbench need the exclusions of the LJ kernels
    if (param->gmxbenchmark && param->force_field != FF_LJ) {
        cmdError(comm, "-gmxbench needs -f lj");
    }
#if !defined(CLUSTER_PAIR) && defined(__SIMD_KERNEL__)
    if (param->gmxbenchmark) {
        cmdError(comm,
            "the verlet-list SIMD kernels have no molecular exclusions for -gmxbench, "
            "build with USE_SIMD_KERNEL=false");
    }
#endif
#if !defined(CLUSTER_PAIR) || defined(CUDA_TARGET)
    if (param->half_lj) {
        cmdError(comm, "-halflj needs the CPU cluster-pair kernels");
    }
#endif
    // Auto (< 0): the GPU verlet list sorts at every neighbor-list build, else no sort
    if (param->resort_every < 0) {
#if defined(CUDA_TARGET) && !defined(CLUSTER_PAIR) && !defined(_MPI)
        param->resort_every = -1;
#else
        param->resort_every = 0;
#endif
    }
#if defined(CUDA_TARGET) && defined(_MPI)
    if (param->resort_every > 0) {
        cmdError(comm, "--sort-every is not supported by the GPU kernels with MPI");
    }
#endif
#ifdef CLUSTER_PAIR
    if (param->resort_every > 0) {
        cmdError(comm,
            "--sort-every is a verlet-list option, the cluster-pair scheme sorts the "
            "atoms into clusters at every neighbor-list build");
    }
    if (param->write_atom_file != NULL) {
        cmdError(comm, "-w is only implemented for the verlet-list scheme");
    }
#ifndef XTC_OUTPUT
    if (param->xtc_file != NULL) {
        cmdError(comm, "--xtc needs a build with XTC_OUTPUT=true");
    }
#endif
#else
    if (param->xtc_file != NULL) {
        cmdError(comm, "--xtc is only implemented for the cluster-pair scheme");
    }
#endif

    if (param->gmxbenchmark) {
        // The GROMACS benchmark builds its pair list with exactly the interaction
        // cut-off (no buffer)
        if (!cutoffGiven) {
            param->cutforce = 1.0;
        }
        if (!skinGiven) {
            param->skin = 0.0;
        }
    }

    param->outer_skin = MAX(param->outer_skin, 0.0);
    param->cutneigh   = param->cutforce + param->skin + param->outer_skin;
}
