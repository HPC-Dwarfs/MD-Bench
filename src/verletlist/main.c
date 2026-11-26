/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>

#include <parameter.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <likwid-marker.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _MPI
#include <mpi.h>
#endif

#include <allocate.h>
#include <atom.h>
#include <balance.h>
#include <cmdline.h>
#include <comm.h>
#include <device.h>
#include <eam.h>
#include <energy.h>
#include <ewald_table.h>
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

extern void copyDataToCUDADevice(Atom*);
extern void copyDataFromCUDADevice(Atom*);

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
    initPbc(atom);
    initStats(stats);
    initNeighbor(neighbor, param);

    if (param->input_file == NULL) {
        if (param->gmxbenchmark) {
            readAtomGmx(atom, param);
        } else {
            createAtom(atom, param);
        }
    } else {
        readAtom(atom, param);
    }
    // Before the atoms are distributed over the MPI ranks, so the ids are global
    initAtomIds(atom);
    // --lj-comb auto and the check of the rule need the atom types
    resolveLJCombRule(param,
        atom->ntypes,
        atom->sqrt_epsilon_per_type,
        atom->sigma3_per_type);
    neighbor->lj_comb_rule = param->lj_comb_rule;
#if COULOMB_IS_EWALD
    param->ewald_beta = computeEwaldCoefficient(param->cutforce, param->ewald_rtol);
#endif
#if COULOMB_RULE == COULOMB_EWALD_TABLE
    initEwaldTable(param);
#endif

#ifndef _MPI
    // Single domain: its box is the simulation box, which the shell methods use to
    // select the ghost atoms of their zones (see ghostZone() in neighbor.c)
    atom->mybox.xprd  = param->xprd;
    atom->mybox.yprd  = param->yprd;
    atom->mybox.zprd  = param->zprd;
    atom->mybox.lo[0] = atom->mybox.lo[1] = atom->mybox.lo[2] = 0.0;
    atom->mybox.hi[0]                                         = param->xprd;
    atom->mybox.hi[1]                                         = param->yprd;
    atom->mybox.hi[2]                                         = param->zprd;
#endif

    setupNeighbor(param);

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

#ifndef CUDA_TARGET
    if (param->resort_every > 0) {
        atom->Nghost = 0;
        sortAtom(atom);
    }
#endif

    // initDevice must precede setupPbc: the CUDA setupPbc reads positions on
    // the device, so the initial H2D upload has to happen first.
    initDevice(param, atom, neighbor);
#if defined(CUDA_TARGET) && !defined(_MPI)
    // Spatial order of the local atoms (--sort-every), before the ghosts are made
    if (param->resort_every != 0) {
        sortAtomsCUDA(atom);
    }
#endif
    setupPbc(atom, param);

#ifdef _MPI
    ghostNeighbor(comm, atom, param);
#ifdef CUDA_TARGET
    copyDataToCUDADevice(atom);
#endif
#else
    updatePbc(atom, param, true);
#endif

    buildNeighbor(atom, neighbor);
    initForce(param);
    timeStop = getTimeStamp();
    return timeStop - timeStart;
}

double reneighbour(
    int n, Parameter* param, Atom* atom, Neighbor* neighbor, Comm* comm, double* timer)
{
    double timeStart, timeStop, stepStart, stepStop;
    timeStart = getTimeStamp();
    LIKWID_MARKER_START("reneighbour");
    // updateAtomsPbc(atom, param, true); function called at updateAtoms
#if defined(CUDA_TARGET) && !defined(_MPI)
    // Spatial order of the local atoms (--sort-every: at every build or every
    // resort_every steps), before the ghosts are made
    if (param->resort_every < 0 ||
        (param->resort_every > 0 && (n + 1) % param->resort_every == 0)) {
        stepStart = getTimeStamp();
        sortAtomsCUDA(atom);
        timer[NEIGH_SORT] += getTimeStamp() - stepStart;
    }
#endif
    stepStart = getTimeStamp();
#ifdef _MPI
    ghostNeighbor(comm, atom, param);
#ifdef CUDA_TARGET
    copyDataToCUDADevice(param, atom);
#endif
#else
    setupPbc(atom, param);
    updatePbc(atom, param, true);
#endif
    stepStop = getTimeStamp();
    timer[NEIGH_PBC] += stepStop - stepStart;

    stepStart = getTimeStamp();
    buildNeighbor(atom, neighbor);
    stepStop = getTimeStamp();
    timer[NEIGH_BUILD] += stepStop - stepStart;

    LIKWID_MARKER_STOP("reneighbour");
    timeStop = getTimeStamp();
    return timeStop - timeStart;
}

void printAtomState(Atom* atom)
{
    printf("Atom counts: Natoms=%d Nlocal=%d Nghost=%d Nmax=%d\n",
        atom->Natoms,
        atom->Nlocal,
        atom->Nghost,
        atom->Nmax);

    int nall = atom->Nlocal + atom->Nghost;
    for (int i = 0; i < nall; i++) {
        printf("%d  %f %f %f\n", i, atom_x(i), atom_y(i), atom_z(i));
    }
}

double updateAtoms(Comm* comm, Atom* atom, Parameter* param)
{
    double timeStart, timeStop;
    timeStart = getTimeStamp();

#ifdef _MPI
#ifdef CUDA_TARGET
    copyDataFromCUDADevice(atom);
#endif
    exchangeComm(comm, atom);
#else
    updateAtomsPbc(atom, param, true);
#endif

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

void writeInput(Parameter* param, Atom* atom)
{
    FILE* fpin = fopen("input.in", "w");
    fprintf(fpin, "0,%f,0,%f,0,%f\n", param->xprd, param->yprd, param->zprd);

    for (int i = 0; i < atom->Nlocal; i++) {
        fprintf(fpin,
            "1,%f,%f,%f,%f,%f,%f\n",
            atom_x(i),
            atom_y(i),
            atom_z(i),
            atom_vx(i),
            atom_vy(i),
            atom_vz(i));
    }

    fclose(fpin);
}

// Writes the forces of the local atoms, ordered by their original index
static void dumpForces(Parameter* param, Atom* atom)
{
    // Collective with MPI, before the non-root processes return
    double elj, ecoul;
    globalEnergies(param, &elj, &ecoul);
#ifdef CUDA_TARGET
    // The forces are on the device
    copyForceFromGPU(atom);
#endif
    MD_FLOAT* f = (MD_FLOAT*)allocate(ALIGNMENT, atom->Natoms * 3 * sizeof(MD_FLOAT));
    for (int i = 0; i < atom->Natoms * 3; i++) {
        f[i] = 0.0;
    }

    for (int i = 0; i < atom->Nlocal; i++) {
        const int id  = atom->id[i];
        f[id * 3 + 0] = atom_fx(i);
        f[id * 3 + 1] = atom_fy(i);
        f[id * 3 + 2] = atom_fz(i);
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
    if (param->compute_energy) {
        fprintf(fp, "# energy_lj %.15g\n", elj);
        fprintf(fp, "# energy_coulomb %.15g\n", ecoul);
    }
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

    long long pairs = 0;
#ifdef CUDA_TARGET
    // GPU builds keep the neighbor counts on the device only
    int* numneigh = (int*)malloc(atom->Nlocal * sizeof(int));
    memcpyFromGPU(numneigh,
        neighbor->d_neighbor.numneigh_inner,
        atom->Nlocal * sizeof(int));
#else
    const int* numneigh = neighbor->numneigh_inner;
#endif
    for (int i = 0; i < atom->Nlocal; i++) {
        pairs += numneigh[i];
    }
#ifdef CUDA_TARGET
    free(numneigh);
#endif

    fprintf(stdout, "\nKERNEL BENCHMARK\n");
    fprintf(stdout, "  Atoms:                 %d\n", atom->Natoms);
    fprintf(stdout, "  Cut-off:               %g\n", rc);
    fprintf(stdout,
        "  Neighbor lists:        %s\n",
        neighbor->half_neigh ? "half" : "full");
    fprintf(stdout, "  Atom pairs in list:    %lld\n", pairs);
    fprintf(stdout, "  Iterations:            %d\n", param->kernel_bench);
    fprintf(stdout, "  Time per iteration:    %.3f us\n", perIter * 1e6);
    fprintf(stdout, "  Useful pairs per iter: %.0f\n", useful);
    fprintf(stdout, "  Computed pairs/us:     %.3f\n", (double)pairs / (perIter * 1e6));
    fprintf(stdout, "  Useful pairs/us:       %.3f\n", useful / (perIter * 1e6));
    if (param->compute_energy) {
        double elj, ecoul;
        globalEnergies(param, &elj, &ecoul);
        fprintf(stdout, "  Energy LJ:             %.15g\n", elj);
        fprintf(stdout, "  Energy Coulomb:        %.15g\n", ecoul);
    }
}

int main(int argc, char** argv)
{
    double timer[NUMTIMER];
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
        LIKWID_MARKER_REGISTER("reneighbour");
        // LIKWID_MARKER_REGISTER("pbc");
    }

    initComm(&argc, &argv, &comm);
    initParameter(&param);
    parseCommandLine(argc, argv, &param, &comm);
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
    traceAddresses(&param, &atom, &neighbor, n + 1); // TODO: trace adress
#endif
    if (param.write_atom_file != NULL) {
        writeAtom(&atom, &param);
    }

    // writeInput(&param, &atom);
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
    timer[NEIGH_SORT]        = 0.0;
    timer[INTEGRATE_INITIAL] = 0.0;
    timer[INTEGRATE_FINAL]   = 0.0;
    timer[THERMO]            = 0.0;
    timer[RENEIGH_CHECK]     = 0.0;
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
        // write_atoms_to_vtk_file(param.vtk_file, &atom, 0);
        printvtk(param.vtk_file, &comm, &atom, &param, 0);
    }
    for (int n = 0; n < param.ntimes; n++) {
        bool resort = param.resort_every > 0 && (n + 1) % param.resort_every == 0;
#ifndef CUDA_TARGET
        if (resort) {
            double sortStart = getTimeStamp();
            DEBUG_MESSAGE("Resorting atoms\n");
            atom.Nghost = 0;
            sortAtom(&atom);
            double sortElapsed = getTimeStamp() - sortStart;
            timer[NEIGH] += sortElapsed;
            timer[NEIGH_SORT] += sortElapsed;
        }
#endif

        // Sorting permutes the atoms and drops the ghosts: the lists have to be rebuilt
        bool reneigh = resort ||
                       (param.reneigh_every > 0 && (n + 1) % param.reneigh_every == 0);
        if (!reneigh && param.displacement_reneigh) {
            double checkStart = getTimeStamp();
            reneigh           = needsReneigh(&atom, &param);
            timer[RENEIGH_CHECK] += getTimeStamp() - checkStart;
        }
        {
            double integrateStart = getTimeStamp();
            initialIntegrate(reneigh, &param, &atom);
            timer[INTEGRATE_INITIAL] += getTimeStamp() - integrateStart;
        }

        if (reneigh) {
            timer[UPDATE] += updateAtoms(&comm, &atom, &param);
            if (param.balance && !((n + 1) % param.balance_every))
                timer[BALANCE] += dynamicBalance(&comm,
                    &grid,
                    &atom,
                    &param,
                    timer[FORCE]);

            timer[NEIGH] += reneighbour(n, &param, &atom, &neighbor, &comm, timer);
            if (param.displacement_reneigh) {
                storeReferencePositions(&atom);
            }
        } else {
            timer[FORWARD] += forward(&comm, &atom, &param);
#ifndef _MPI
            updatePbc(&atom, &param, false);
#endif
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
        }

#if defined(MEM_TRACER) || defined(INDEX_TRACER)
        traceAddresses(&param, &atom, &neighbor, n + 1);
#endif

        timer[FORCE] += computeForce(&param, &atom, &neighbor, &stats);
        timer[REVERSE] += reverse(&comm, &atom, &param);
        {
            double integrateStart = getTimeStamp();
            finalIntegrate(reneigh, &param, &atom);
            timer[INTEGRATE_FINAL] += getTimeStamp() - integrateStart;
        }

        if (!((n + 1) % param.nstat) && (n + 1) < param.ntimes) {
#ifdef CUDA_TARGET
            memcpyFromGPU(atom.x, atom.d_atom.x, atom.Nmax * sizeof(MD_FLOAT) * 3);
#endif
            double thermoStart = getTimeStamp();
            computeThermo(n + 1, &param, &atom);
            timer[THERMO] += getTimeStamp() - thermoStart;
        }

        if (param.vtk_file != NULL) {
            // write_atoms_to_vtk_file(param.vtk_file, &atom, n + 1);
            printvtk(param.vtk_file, &comm, &atom, &param, n + 1);
        }
    }
    barrierComm();
    timer[TOTAL] = getTimeStamp() - timer[TOTAL];
    computeThermo(-1, &param, &atom);
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
            printTimerSubRow("Pbc/ghost", NEIGH_PBC, sumt, mint, maxt, timer[TOTAL], n);
            printTimerSubRow("Building", NEIGH_BUILD, sumt, mint, maxt, timer[TOTAL], n);
            printTimerSubRow("Prune", NEIGH_PRUNE, sumt, mint, maxt, timer[TOTAL], n);
            printTimerSubRow("Sort", NEIGH_SORT, sumt, mint, maxt, timer[TOTAL], n);
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

#ifdef CUDA_TARGET
    if (param.force_field == FF_LJ_TABLE) {
        freeLJTableGPU();
    }
#endif

    endComm(&comm);
    LIKWID_MARKER_CLOSE;
    return EXIT_SUCCESS;
}
