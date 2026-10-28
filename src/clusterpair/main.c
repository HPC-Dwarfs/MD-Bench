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
#include <ljtable.h>
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
extern void pruneNeighborCUDA(Parameter*, Atom*, Neighbor*);

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
    if (param->force_field == FF_LJ_TABLE) {
        initLJTable(param);
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
        if (param->gmxbenchmark) {
#if LJ_COMB_RULE == LJ_COMB_SINGLE
            fprintf(stderr, "-gmxbench needs LJ_COMB_RULE=geometric or full\n");
            exit(EXIT_FAILURE);
#endif
            readAtomGmx(atom, param);
        } else {
            createAtom(atom, param);
        }
    } else {
        readAtom(atom, param);
    }
    // Before the atoms are distributed over the MPI ranks, so the ids are global
    initAtomIds(atom);
#if COULOMB_RULE == COULOMB_EWALD
    param->ewald_beta = computeEwaldCoefficient(param->cutforce, param->ewald_rtol);
#endif
#ifndef _MPI
    // Single domain: its box is the simulation box, which the half-shell selection of the
    // periodic ghost clusters (see initForce()) works on
    atom->mybox.xprd  = param->xprd;
    atom->mybox.yprd  = param->yprd;
    atom->mybox.zprd  = param->zprd;
    atom->mybox.lo[0] = atom->mybox.lo[1] = atom->mybox.lo[2] = 0.0;
    atom->mybox.hi[0]                                         = param->xprd;
    atom->mybox.hi[1]                                         = param->yprd;
    atom->mybox.hi[2]                                         = param->zprd;
#endif
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
    updatePbcCPU(atom, param, true);
#endif
    binJClusters(param, atom);
    buildNeighbor(atom, neighbor);
    initDevice(param, atom, neighbor);
    timeStop = getTimeStamp();
    return timeStop - timeStart;
}

double reneighbour(
    Comm* comm, Parameter* param, Atom* atom, Neighbor* neighbor, double* timer)
{
    double timeStart, timeStop, stepStart, stepStop;
    timeStart = getTimeStamp();
    LIKWID_MARKER_START("reneighbour");
    // updateSingleAtoms(param, atom);
    // updateAtomsPbc(atom, param, false);
    stepStart = getTimeStamp();
    buildClusters(atom);
    defineJClusters(param, atom);
    stepStop = getTimeStamp();
    timer[NEIGH_CLUSTERS] += stepStop - stepStart;

    stepStart = getTimeStamp();
#ifdef _MPI
    ghostNeighbor(comm, atom, param);
#else
    setupPbc(atom, param);
    updatePbcCPU(atom, param, true);
#endif
    stepStop = getTimeStamp();
    timer[NEIGH_PBC] += stepStop - stepStart;

    stepStart = getTimeStamp();
    binJClusters(param, atom);
    stepStop = getTimeStamp();
    timer[NEIGH_BIN] += stepStop - stepStart;

    stepStart = getTimeStamp();
    buildNeighbor(atom, neighbor);
    stepStop = getTimeStamp();
    timer[NEIGH_BUILD] += stepStop - stepStart;

    LIKWID_MARKER_STOP("reneighbour");
    timeStop = getTimeStamp();
    return timeStop - timeStart;
}

// Prints one indented sub-row of the timing breakdown table, nested under its
// parent aggregate row (only shown with -v/--verbose). Indent (6) + label
// width (18) match the parent rows' indent (4) + label width (20) so the
// numeric columns line up.
static void printTimerSubRow(const char* label,
    timertype idx,
    double* sumt,
    double* mint,
    double* maxt,
    double total,
    int nproc)
{
    double avg = sumt[idx] / nproc;
    fprintf(stdout,
        "      %-18s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
        label,
        avg,
        mint[idx],
        maxt[idx],
        100.0 * sumt[idx] / (nproc * total),
        avg > 0.0 ? 100.0 * (maxt[idx] - mint[idx]) / avg : 0.0);
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

// Writes the forces of the local atoms, ordered by their original index
static void dumpForces(Parameter* param, Atom* atom)
{
#ifdef CUDA_TARGET
    fprintf(stderr, "--dump-forces is not supported for GPU builds\n");
    return;
#endif
    MD_FLOAT* f = (MD_FLOAT*)allocate(ALIGNMENT, atom->Natoms * 3 * sizeof(MD_FLOAT));
    for (int i = 0; i < atom->Natoms * 3; i++) {
        f[i] = 0.0;
    }

    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        const int ci_vec_base = CI_VECTOR3_BASE_INDEX(ci);
        const int ci_sca_base = CI_SCALAR_BASE_INDEX(ci);
        for (int cii = 0; cii < atom->iclusters[ci].natoms; cii++) {
            const int id  = atom->cl_id[ci_sca_base + cii];
            f[id * 3 + 0] = atom->cl_f[ci_vec_base + CL_X_INDEX_3D(cii)];
            f[id * 3 + 1] = atom->cl_f[ci_vec_base + CL_Y_INDEX_3D(cii)];
            f[id * 3 + 2] = atom->cl_f[ci_vec_base + CL_Z_INDEX_3D(cii)];
        }
    }

#ifdef _MPI
    // Every atom is local to exactly one rank, the others contribute zeros
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Reduce(me == 0 ? MPI_IN_PLACE : f,
        f,
        atom->Natoms * 3,
        sizeof(MD_FLOAT) == 4 ? MPI_FLOAT : MPI_DOUBLE,
        MPI_SUM,
        0,
        MPI_COMM_WORLD);
    if (me != 0) {
        free(f);
        return;
    }
#endif

    FILE* fp = fopen(param->dump_forces_file, "w");
    if (fp == NULL) {
        perror("Could not open the force dump file");
        free(f);
        return;
    }

    fprintf(fp, "# natoms %d\n", atom->Natoms);
    for (int i = 0; i < atom->Natoms; i++) {
        fprintf(fp, "%d %.9g %.9g %.9g\n", i, f[i * 3 + 0], f[i * 3 + 1], f[i * 3 + 2]);
    }

    free(f);
    fclose(fp);
}

// Runs only the force kernel on the fixed neighbor lists, like the GROMACS non-bonded
// benchmark does, and reports the time per kernel call
static void runKernelBench(Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    double time = 0.0;
    for (int it = 0; it < param->kernel_bench; it++) {
        time += computeForce(param, atom, neighbor, stats);
    }

    // Same estimate of the number of useful pairs as the GROMACS benchmark uses
    const double volume   = (double)param->xprd * param->yprd * param->zprd;
    const double density  = atom->Natoms / volume;
    const double rc       = param->cutforce;
    const double pairsCut = density * 4.0 / 3.0 * M_PI * rc * rc * rc;
    const double useful   = atom->Natoms * 0.5 * (pairsCut + 1.0);
    const double perIter  = time / param->kernel_bench;

    long long clusterPairs = 0;
    for (int ci = 0; ci < atom->Nclusters_local; ci++) {
        clusterPairs += neighbor->numneigh_inner[ci];
    }

    fprintf(stdout, "\nKERNEL BENCHMARK\n");
    fprintf(stdout, "  Atoms:                 %d\n", atom->Natoms);
    fprintf(stdout, "  Cut-off:               %g\n", rc);
    fprintf(stdout, "  Neighbor lists:        %s\n", param->half_neigh ? "half" : "full");
    fprintf(stdout,
        "  Cluster pairs in list: %lld (%dx%d, %lld atom pairs)\n",
        clusterPairs,
        CLUSTER_M,
        CLUSTER_N,
        clusterPairs * CLUSTER_M * CLUSTER_N);
    fprintf(stdout, "  Iterations:            %d\n", param->kernel_bench);
    fprintf(stdout, "  Time per iteration:    %.3f us\n", perIter * 1e6);
    fprintf(stdout, "  Useful pairs per iter: %.0f\n", useful);
    fprintf(stdout,
        "  Computed pairs/us:     %.3f\n",
        (double)clusterPairs * CLUSTER_M * CLUSTER_N / (perIter * 1e6));
    fprintf(stdout, "  Useful pairs/us:       %.3f\n", useful / (perIter * 1e6));
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

    int cutoffGiven = 0, skinGiven = 0;
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
            cutoffGiven    = 1;
            continue;
        }
        if ((strcmp(argv[i], "-s") == 0) || (strcmp(argv[i], "--skin") == 0)) {
            param.skin = atof(argv[++i]);
            skinGiven  = 1;
            continue;
        }
        if ((strcmp(argv[i], "--outer-skin") == 0)) {
            param.outer_skin = atof(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "--prune-every") == 0)) {
            param.prune_every = atoi(argv[++i]);
            continue;
        }
        if ((strcmp(argv[i], "--reneigh-every") == 0)) {
            param.reneigh_every = atoi(argv[++i]);
            continue;
        }
        if (strcmp(argv[i], "--displacement-reneigh") == 0) {
            param.displacement_reneigh = 1;
            continue;
        }
        if ((strcmp(argv[i], "-v") == 0) || (strcmp(argv[i], "--verbose") == 0)) {
            param.verbose = 1;
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

        if ((strcmp(argv[i], "-gmxbench") == 0)) {
            param.gmxbenchmark = true;
            if (param.input_file != NULL) {
                free(param.input_file);
                param.input_file = NULL;
            }
            continue;
        }

        if ((strcmp(argv[i], "-size") == 0)) {
            param.size = atoi(argv[++i]);
            continue;
        }

        if ((strcmp(argv[i], "--kernel-bench") == 0)) {
            param.kernel_bench = atoi(argv[++i]);
            continue;
        }

        if ((strcmp(argv[i], "--dump-forces") == 0)) {
            param.dump_forces_file = strdup(argv[++i]);
            continue;
        }

        if ((strcmp(argv[i], "--ntypes") == 0)) {
            param.ntypes = atoi(argv[++i]);
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
            printf("-s / --skin <real>:   set skin (Verlet buffer for the "
                   "force-evaluation list)\n");
            printf("--outer-skin <real>:  additive outer skin (outer cutoff = cutforce + "
                   "skin + outer_skin; >0 enables double-cutoff)\n");
            printf("--prune-every <int>:  pruning frequency (steps between inner-list "
                   "refresh)\n");
            printf("--reneigh-every <int>: neighbor-list rebuild frequency (steps)\n");
            printf("--displacement-reneigh: rebuild neighbor list based on atom "
                   "displacement "
                   "(skin/2 threshold); use --reneigh-every 0 to disable fixed-interval "
                   "rebuild\n");
            printf("--freq <real>:        processor frequency (GHz)\n");
            printf("--vtk <string>:       VTK file for visualization\n");
            printf("--xtc <string>:       XTC file for visualization\n");
            printf("--ntypes <int>:       set number of atom types\n");
            printf("-gmxbench:            use the water system of the GROMACS non-bonded "
                   "benchmark (cut-off 1.0 and no skin unless -r/-s are given)\n");
            printf("-size <int>:          replicate the -gmxbench system this many times "
                   "(power of two)\n");
            printf(
                "--kernel-bench <int>: after setup, only run the force kernel this many "
                "times on the fixed neighbor lists and report its time\n");
            printf("--dump-forces <file>: write the forces of a single kernel pass "
                   "(atom id, fx, fy, fz) to a text file\n");
            printf("-v / --verbose:       print a per-routine timing breakdown in the "
                   "performance report\n");
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

    if (param.gmxbenchmark) {
        // The GROMACS benchmark builds its pair list with exactly the interaction
        // cut-off (no buffer)
        if (!cutoffGiven) {
            param.cutforce = 1.0;
        }
        if (!skinGiven) {
            param.skin = 0.0;
        }
    }

    param.outer_skin = MAX(param.outer_skin, 0.0);
    param.cutneigh   = param.cutforce + param.skin + param.outer_skin;
#ifdef CUDA_TARGET
    // Must run before setup(), which allocates pinned host memory
    // (growClusters()) before initDevice() is reached.
    initDeviceContext();
#endif
    timer[SETUP] = setup(&param, &eam, &atom, &neighbor, &stats, &comm, &grid);
    if (param.displacement_reneigh) {
        storeReferencePositions(&atom);
    }

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
    timer[TOTAL]             = getTimeStamp();
    timer[FORCE]             = computeForce(&param, &atom, &neighbor, &stats);
    timer[NEIGH]             = 0.0;
    timer[FORWARD]           = 0.0;
    timer[UPDATE]            = 0.0;
    timer[BALANCE]           = 0.0;
    timer[NEIGH_PBC]         = 0.0;
    timer[NEIGH_BUILD]       = 0.0;
    timer[NEIGH_PRUNE]       = 0.0;
    timer[NEIGH_CLUSTERS]    = 0.0;
    timer[NEIGH_BIN]         = 0.0;
    timer[INTEGRATE_INITIAL] = 0.0;
    timer[INTEGRATE_FINAL]   = 0.0;
    timer[THERMO]            = 0.0;
    timer[RENEIGH_CHECK]     = 0.0;
    timer[DEVICE_COPY]       = 0.0;
    timer[REVERSE]           = reverse(&comm, &atom, &param);

    if (param.dump_forces_file != NULL) {
        dumpForces(&param, &atom);
    }

    if (param.kernel_bench > 0) {
        runKernelBench(&param, &atom, &neighbor, &stats);
        endComm(&comm);
        LIKWID_MARKER_CLOSE;
        return EXIT_SUCCESS;
    }

    if (param.vtk_file != NULL) {
        printvtk(param.vtk_file, &comm, &atom, &param, 0);
    }

    // TODO: modify xct
    if (param.xtc_file != NULL) {
        xtc_init(param.xtc_file, &atom, 0);
    }

    for (int n = 0; n < param.ntimes; n++) {
        {
            double integrateStart = getTimeStamp();
            initialIntegrate(&param, &atom);
            timer[INTEGRATE_INITIAL] += getTimeStamp() - integrateStart;
        }

        bool reneigh = (param.reneigh_every > 0 && !((n + 1) % param.reneigh_every));
        if (!reneigh && param.displacement_reneigh) {
            double checkStart = getTimeStamp();
            reneigh           = needsReneigh(&atom, &param);
            timer[RENEIGH_CHECK] += getTimeStamp() - checkStart;
        }

        if (!reneigh) {
            if (param.outer_skin > 0.0 && !((n + 1) % param.prune_every)) {
                double prune_start = getTimeStamp();
#ifdef CUDA_TARGET
                pruneNeighborCUDA(&param, &atom, &neighbor);
#else
                pruneNeighbor(&param, &atom, &neighbor);
#endif
                double pruneElapsed = getTimeStamp() - prune_start;
                timer[NEIGH] += pruneElapsed;
                timer[NEIGH_PRUNE] += pruneElapsed;
            }

            timer[FORWARD] += forward(&comm, &atom, &param);
        } else {
#ifdef CUDA_TARGET
            double copyStart = getTimeStamp();
            copyDataFromCUDADevice(&param, &atom);
            timer[DEVICE_COPY] += getTimeStamp() - copyStart;
#endif
            timer[UPDATE] += updateAtoms(&comm, &atom, &param);
            if (param.balance && !((n + 1) % param.balance_every)) {
                timer[BALANCE] += dynamicBalance(&comm,
                    &grid,
                    &atom,
                    &param,
                    timer[FORCE]);
            }

            timer[NEIGH] += reneighbour(&comm, &param, &atom, &neighbor, timer);
            if (param.displacement_reneigh) {
                storeReferencePositions(&atom);
            }
#ifdef CUDA_TARGET
            copyStart = getTimeStamp();
            copyDataToCUDADevice(&param, &atom, &neighbor);
            timer[DEVICE_COPY] += getTimeStamp() - copyStart;
#endif
        }
#if defined(MEM_TRACER) || defined(INDEX_TRACER)
        traceAddresses(&param, &atom, &neighbor, n + 1);
#endif
        timer[FORCE] += computeForce(&param, &atom, &neighbor, &stats);
        timer[REVERSE] += reverse(&comm, &atom, &param);
        {
            double integrateStart = getTimeStamp();
            finalIntegrate(&param, &atom);
            timer[INTEGRATE_FINAL] += getTimeStamp() - integrateStart;
        }

        if (!((n + 1) % param.nstat) && (n + 1) < param.ntimes) {
            double thermoStart = getTimeStamp();
            computeThermo(n + 1, &param, &atom);
            timer[THERMO] += getTimeStamp() - thermoStart;
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
        if (param.verbose) {
            printTimerSubRow("Clusters",
                NEIGH_CLUSTERS,
                sumt,
                mint,
                maxt,
                timer[TOTAL],
                n);
            printTimerSubRow("Pbc/ghost", NEIGH_PBC, sumt, mint, maxt, timer[TOTAL], n);
            printTimerSubRow("Binning", NEIGH_BIN, sumt, mint, maxt, timer[TOTAL], n);
            printTimerSubRow("Building", NEIGH_BUILD, sumt, mint, maxt, timer[TOTAL], n);
            printTimerSubRow("Prune", NEIGH_PRUNE, sumt, mint, maxt, timer[TOTAL], n);
        }
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Rest",
            sumt[REST] / n,
            mint[REST],
            maxt[REST],
            100.0 * sumt[REST] / (n * timer[TOTAL]),
            100.0 * (maxt[REST] - mint[REST]) / (sumt[REST] / n));
        if (param.verbose) {
            printTimerSubRow("Integrate initial",
                INTEGRATE_INITIAL,
                sumt,
                mint,
                maxt,
                timer[TOTAL],
                n);
            printTimerSubRow("Integrate final",
                INTEGRATE_FINAL,
                sumt,
                mint,
                maxt,
                timer[TOTAL],
                n);
            printTimerSubRow("Thermo output", THERMO, sumt, mint, maxt, timer[TOTAL], n);
            printTimerSubRow("Reneigh check",
                RENEIGH_CHECK,
                sumt,
                mint,
                maxt,
                timer[TOTAL],
                n);
#ifdef CUDA_TARGET
            printTimerSubRow("Device copy",
                DEVICE_COPY,
                sumt,
                mint,
                maxt,
                timer[TOTAL],
                n);
#endif
        }
        fprintf(stdout,
            "    %-20s %8.2f   %8.2f   %8.2f    %5.1f%%       %5.1f%%\n",
            "Integration",
            sumt[UPDATE] / n,
            mint[UPDATE],
            maxt[UPDATE],
            100.0 * sumt[UPDATE] / (n * timer[TOTAL]),
            100.0 * (maxt[UPDATE] - mint[UPDATE]) / (sumt[UPDATE] / n));
        if (param.verbose) {
            printTimerSubRow("Load balance", BALANCE, sumt, mint, maxt, timer[TOTAL], n);
        }
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
