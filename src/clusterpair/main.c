/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <likwid-marker.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <allocate.h>
#include <atom.h>
#include <balance.h>
#include <comm.h>
#include <device.h>
#include <eam.h>
#include <force.h>
#include <grid.h>
#include <integrate.h>
#include <neighbor.h>
#include <parameter.h>
#include <pbc.h>
#include <shell_methods.h>
#include <stats.h>
#include <thermo.h>
#include <timers.h>
#include <timing.h>
#include <util.h>
#include <vtk.h>
#include <xtc.h>

extern void copyDataToCUDADevice(Parameter*, Atom*, Neighbor*);
extern void copyDataFromCUDADevice(Parameter*, Atom*);
extern void cudaDeviceFree(Parameter*);

#define HLINE                                                                            \
    "----------------------------------------------------------------------------\n"

double setup(Parameter* param,
    Eam* eam,
    Atom* atom,
    Neighbor* neighbor,
    Stats* stats,
    Comm* comm,
    Grid* grid)
{
    if (param->force_field == FF_EAM) {
        initEam(param);
    }

    double timeStart, timeStop;
    param->lattice = pow((4.0 / param->rho), (1.0 / 3.0));
    param->xprd    = param->nx * param->lattice;
    param->yprd    = param->ny * param->lattice;
    param->zprd    = param->nz * param->lattice;

    timeStart = getTimeStamp();
    initAtom(atom);
    initForce(param);
    initPbc(atom);
    initStats(stats);
    initNeighbor(neighbor, param);
    if (param->input_file == NULL) {
        createAtom(atom, param);
    } else {
        readAtom(atom, param);
    }
    setupNeighbor(param, atom);
#ifdef _MPI
    setupGrid(grid, atom, param);
    setupComm(comm, param, grid);
    if (param->balance) {
        initialBalance(param, atom, neighbor, stats, comm, grid);
    }
#endif
    setupThermo(param, atom->Natoms);
    if (param->input_file == NULL) {
        adjustThermo(param, atom);
    }
    buildClusters(atom);
    defineJClusters(param, atom);
#ifdef _MPI
    ghostNeighbor(comm, atom, param);
#else
    setupPbc(atom, param);
#endif
    binJClusters(param, atom);
    buildNeighbor(atom, neighbor);
    initDevice(param, atom, neighbor);
    timeStop = getTimeStamp();
    return timeStop - timeStart;
}

double reneighbour(Comm* comm, Parameter* param, Atom* atom, Neighbor* neighbor)
{
    double timeStart, timeStop;
    timeStart = getTimeStamp();
    LIKWID_MARKER_START("reneighbour");
    // updateSingleAtoms(param, atom);
    // updateAtomsPbc(atom, param, false);
    buildClusters(atom);
    defineJClusters(param, atom);
#ifdef _MPI
    ghostNeighbor(comm, atom, param);
#else
    setupPbc(atom, param);
#endif
    binJClusters(param, atom);
    buildNeighbor(atom, neighbor);
    LIKWID_MARKER_STOP("reneighbour");
    timeStop = getTimeStamp();
    return timeStop - timeStart;
}

double updateAtoms(Comm* comm, Atom* atom, Parameter* param)
{
    double timeStart, timeStop;
    timeStart = getTimeStamp();
    updateSingleAtoms(param, atom);
#ifdef _MPI
    exchangeComm(comm, atom);
#else
    updateAtomsPbc(atom, param, false);
#endif
    timeStop = getTimeStamp();
    return timeStop - timeStart;
}

int main(int argc, char** argv)
{
    double timer[NUMTIMER];
    Eam eam;
    Atom atom;
    Neighbor neighbor;
    Stats stats;
    Parameter param;
    Comm comm;
    Grid grid;
    LIKWID_MARKER_INIT;
#pragma omp parallel
    {
        LIKWID_MARKER_REGISTER("force");
        // LIKWID_MARKER_REGISTER("reneighbour");
        // LIKWID_MARKER_REGISTER("pbc");
    }

    initComm(&argc, &argv, &comm);
    initParameter(&param);
    for (int i = 0; i < argc; i++) {
        if ((strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "--param") == 0)) {
            readParameter(&param, argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-f") == 0)) {
            if ((param.force_field = str2ff(argv[++i])) < 0) {
                fprintf(stderr, "Invalid force field!\n");
                exit(-1);
            }
            continue;
        }
        if ((strcmp(argv[i], "-i") == 0)) {
            param.input_file = strdup(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-e") == 0)) {
            param.eam_file = strdup(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-n") == 0) || (strcmp(argv[i], "--nsteps") == 0)) {
            param.ntimes = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-nx") == 0)) {
            param.nx = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-ny") == 0)) {
            param.ny = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-nz") == 0)) {
            param.nz = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-half") == 0)) {
            param.half_neigh = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-method") == 0)) {
            param.method = atoi(argv[++i]);
            if (param.method > 2 || param.method < 0) {
                fprintf_once(comm.myproc, stderr, "Method does not exist!\n");
                endComm(&comm);
                exit(0);
            }
            continue;
        }
        if ((strcmp(argv[i], "-bal") == 0)) {
            param.balance = atoi(argv[++i]);

            if (param.balance > 3 || param.balance < 0) {
                fprintf_once(comm.myproc, stderr, "Load balance does not exist!\n");
                endComm(&comm);
                exit(0);
            }

            continue;
        }
        if ((strcmp(argv[i], "-m") == 0) || (strcmp(argv[i], "--mass") == 0)) {
            param.mass = atof(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-r") == 0) || (strcmp(argv[i], "--radius") == 0)) {
            param.cutforce = atof(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "-s") == 0) || (strcmp(argv[i], "--skin") == 0)) {
            param.skin = atof(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "--freq") == 0)) {
            param.proc_freq = atof(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "--vtk") == 0)) {
            param.vtk_file = strdup(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "--xtc") == 0)) {
#ifndef XTC_OUTPUT
            fprintf(stderr,
                "XTC not available, set XTC_OUTPUT option in config.mk file and "
                "recompile MD-Bench!");
            exit(-1);
#else
            param.xtc_file = strdup(argv[++i]);
#endif
            continue;
        }

        if ((strcmp(argv[i], "-setup") == 0)) {
            param.setup = atoi(argv[++i]);
            continue;
        }

        if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0)) {
            printf("MD Bench: A minimalistic re-implementation of miniMD\n");
            printf(HLINE);
            printf("-p <string>:          file to read parameters from (can be specified "
                   "more than once)\n");
            printf("-f <string>:          force field (lj or eam), default lj\n");
            printf("-i <string>:          input file with atom positions (dump)\n");
            printf("-e <string>:          input file for EAM\n");
            printf("-n / --nsteps <int>:  set number of timesteps for simulation\n");
            printf("-nx/-ny/-nz <int>:    set linear dimension of systembox in x/y/z "
                   "direction\n");
            printf("-r / --radius <real>: set cutoff radius\n");
            printf("-s / --skin <real>:   set skin (verlet buffer)\n");
            printf("--freq <real>:        processor frequency (GHz)\n");
            printf("--vtk <string>:       VTK file for visualization\n");
            printf("--xtc <string>:       XTC file for visualization\n");
            printf(HLINE);
            exit(EXIT_SUCCESS);
        }
    }

    if (param.balance > 0 && param.method == 1) {
        fprintf_once(comm.myproc,
            stderr,
            "Half Shell is not supported with load balance!\n");
        endComm(&comm);
        exit(0);
    }

    param.cutneigh = param.cutforce + param.skin;
    timer[SETUP]   = setup(&param, &eam, &atom, &neighbor, &stats, &comm, &grid);

    if (comm.myproc == 0) {
        printParameter(&param);
    }

    fprintf_once(comm.myproc, stdout, "\n");
    fprintf_once(comm.myproc, stdout, "SIMULATION PROGRESS\n");
    fprintf_once(comm.myproc,
        stdout,
        "-------------------------------------------------------------------------------"
        "\n");
    fprintf_once(comm.myproc,
        stdout,
        "  %-10s %15s %15s\n",
        "Step",
        "Temperature",
        "Pressure");
    fflush(stdout);
    computeThermo(0, &param, &atom);
#if defined(MEM_TRACER) || defined(INDEX_TRACER)
    traceAddresses(&param, &atom, &neighbor, n + 1);
#endif

#ifdef CUDA_TARGET
    copyDataToCUDADevice(&param, &atom, &neighbor);
#endif
    barrierComm();
    timer[TOTAL]   = getTimeStamp();
    timer[FORCE]   = computeForce(&param, &atom, &neighbor, &stats);
    timer[NEIGH]   = 0.0;
    timer[FORWARD] = 0.0;
    timer[UPDATE]  = 0.0;
    timer[BALANCE] = 0.0;
    timer[REVERSE] = reverse(&comm, &atom, &param);

    if (param.vtk_file != NULL) {
        printvtk(param.vtk_file, &comm, &atom, &param, 0);
    }

    // TODO: modify xct
    if (param.xtc_file != NULL) {
        xtc_init(param.xtc_file, &atom, 0);
    }

    for (int n = 0; n < param.ntimes; n++) {
        initialIntegrate(&param, &atom);

        if ((n + 1) % param.reneigh_every) {
            if (!((n + 1) % param.prune_every)) {
                pruneNeighbor(&param, &atom, &neighbor);
            }

            timer[FORWARD] += forward(&comm, &atom, &param);
        } else {
#ifdef CUDA_TARGET
            copyDataFromCUDADevice(&param, &atom);
#endif
            timer[UPDATE] += updateAtoms(&comm, &atom, &param);
            if (param.balance && !((n + 1) % param.balance_every)) {
                timer[BALANCE] += dynamicBalance(&comm,
                    &grid,
                    &atom,
                    &param,
                    timer[FORCE]);
            }

            timer[NEIGH] += reneighbour(&comm, &param, &atom, &neighbor);
#ifdef CUDA_TARGET
            copyDataToCUDADevice(&param, &atom, &neighbor);
#endif
        }
#if defined(MEM_TRACER) || defined(INDEX_TRACER)
        traceAddresses(&param, &atom, &neighbor, n + 1);
#endif
        timer[FORCE] += computeForce(&param, &atom, &neighbor, &stats);
        timer[REVERSE] += reverse(&comm, &atom, &param);
        finalIntegrate(&param, &atom);

        if (!((n + 1) % param.nstat) && (n + 1) < param.ntimes) {
            computeThermo(n + 1, &param, &atom);
        }

        int writePos = !((n + 1) % param.x_out_every);
        int writeVel = !((n + 1) % param.v_out_every);
        if (writePos || writeVel) {
            if (param.vtk_file != NULL) {
#ifdef CUDA_TARGET
                copyDataFromCUDADevice(&param, &atom);
#endif
                printvtk(param.vtk_file, &comm, &atom, &param, n + 1);
            }

            // TODO: xtc file
            if (param.xtc_file != NULL) {
                xtc_write(&atom, n + 1, write_pos, write_vel);
            }
        }
    }

#ifdef CUDA_TARGET
    copyDataFromCUDADevice(&param, &atom);
#endif
    barrierComm();
    timer[TOTAL] = getTimeStamp() - timer[TOTAL];
    updateAtoms(&comm, &atom, &param);
    computeThermo(-1, &param, &atom);
    // TODO: xtc file
    if (param.xtc_file != NULL) {
        xtc_end();
    }

#ifdef CUDA_TARGET
    cudaDeviceFree(&param);
#endif
    timer[REST] = timer[TOTAL] - timer[FORCE] - timer[NEIGH] - timer[BALANCE] -
                  timer[FORWARD] - timer[REVERSE];
#ifdef _MPI
    double mint[NUMTIMER];
    double maxt[NUMTIMER];
    double sumt[NUMTIMER];
    int Nghost = atom.Nghost;
    MPI_Reduce(timer, mint, NUMTIMER, MPI_DOUBLE, MPI_MIN, 0, world);
    MPI_Reduce(timer, maxt, NUMTIMER, MPI_DOUBLE, MPI_MAX, 0, world);
    MPI_Reduce(timer, sumt, NUMTIMER, MPI_DOUBLE, MPI_SUM, 0, world);
    MPI_Reduce(&atom.Nghost, &Nghost, 1, MPI_INT, MPI_SUM, 0, world);
#else
    int Nghost   = atom.Nghost;
    double* mint = timer;
    double* maxt = timer;
    double* sumt = timer;
#endif

    if (comm.myproc == 0) {
        int n         = comm.numproc;
        double ns_day = (param.ntimes * param.dt * 1e-6 * 86400.0) / timer[TOTAL];
        fprintf_once(comm.myproc,
            stdout,
            "----------------------------------------------------------------------------"
            "---\n");
        fprintf(stdout, "\n");
        fprintf(stdout, "PERFORMANCE REPORT\n");
        fprintf(stdout,
            "----------------------------------------------------------------------------"
            "---\n");
        fprintf(stdout, "  Timing Breakdown\n");
        fprintf(stdout,
            "                          Avg (s)    Min (s)    Max (s)    %% Time    "
            "Imbalance\n");
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Force",
            sumt[FORCE] / n,
            mint[FORCE],
            maxt[FORCE],
            100.0 * sumt[FORCE] / (n * timer[TOTAL]),
            100.0 * (maxt[FORCE] - mint[FORCE]) / (sumt[FORCE] / n));
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Neighbor",
            sumt[NEIGH] / n,
            mint[NEIGH],
            maxt[NEIGH],
            100.0 * sumt[NEIGH] / (n * timer[TOTAL]),
            100.0 * (maxt[NEIGH] - mint[NEIGH]) / (sumt[NEIGH] / n));
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Rest",
            sumt[REST] / n,
            mint[REST],
            maxt[REST],
            100.0 * sumt[REST] / (n * timer[TOTAL]),
            100.0 * (maxt[REST] - mint[REST]) / (sumt[REST] / n));
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Integration",
            sumt[UPDATE] / n,
            mint[UPDATE],
            maxt[UPDATE],
            100.0 * sumt[UPDATE] / (n * timer[TOTAL]),
            100.0 * (maxt[UPDATE] - mint[UPDATE]) / (sumt[UPDATE] / n));
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Setup",
            sumt[SETUP] / n,
            mint[SETUP],
            maxt[SETUP],
            100.0 * sumt[SETUP] / (n * timer[TOTAL]),
            100.0 * (maxt[SETUP] - mint[SETUP]) / (sumt[SETUP] / n));
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Reverse comm",
            sumt[REVERSE] / n,
            mint[REVERSE],
            maxt[REVERSE],
            100.0 * sumt[REVERSE] / (n * timer[TOTAL]),
            100.0 * (maxt[REVERSE] - mint[REVERSE]) / (sumt[REVERSE] / n));
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Forward comm",
            sumt[FORWARD] / n,
            mint[FORWARD],
            maxt[FORWARD],
            100.0 * sumt[FORWARD] / (n * timer[TOTAL]),
            100.0 * (maxt[FORWARD] - mint[FORWARD]) / (sumt[FORWARD] / n));
        fprintf(stdout,
            "\n  System: %d atoms (%d ghost) | %d timesteps\n",
            atom.Natoms,
            Nghost,
            param.ntimes);
        fprintf(stdout,
            "  Performance: %.2fs total | %.2f atom updates/us | %.2f steps/s | %.2f "
            "ns/day\n",
            timer[TOTAL],
            (double)atom.Natoms * param.ntimes / (timer[TOTAL] * 1e6),
            param.ntimes / timer[TOTAL],
            ns_day);
        fprintf(stdout,
            "----------------------------------------------------------------------------"
            "---\n");
    }

#ifdef COMPUTE_STATS
    displayStatistics(&atom, &param, &stats, timer);
#endif

    endComm(&comm);
    LIKWID_MARKER_CLOSE;
    return EXIT_SUCCESS;
}
