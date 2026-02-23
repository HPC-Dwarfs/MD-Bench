/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

#include <atom.h>
#include <box.h>
#include <parameter.h>

#ifndef __MAP_H_
#define __MAP_H_

#define atom_pos(i) ((dim == 0) ? atom_x((i)) : (dim == 1) ? atom_y((i)) : atom_z((i)))

#ifdef _MPI
    #include <mpi.h>
    #define world       MPI_COMM_WORLD
    static MPI_Datatype type_float = (sizeof(MD_FLOAT) == 4) ? MPI_FLOAT : MPI_DOUBLE;
#endif

typedef struct {
    int balance_every;
    int map_size;
    MD_FLOAT* map;
    //===Param for Staggerd balance
    int nprocs[3];
    int coord[3];
    MD_FLOAT cutneigh[3];
    double Timer;
    //===Param for RCB balance
    MD_FLOAT* buf_send;
    MD_FLOAT* buf_recv;
    int maxsend;
    int maxrecv;
} Grid;

enum { RCB = 1, meanTimeRCB, Staggered };
#ifdef _MPI
void setupGrid(Grid*, Atom*, Parameter*);
void printGrid(Grid*);
typedef MD_FLOAT (*RCB_Method)(Atom*, MPI_Comm, int, double);
void cartisian3d(Grid*, Parameter*, Box*);
void rcbBalance(Grid*, Atom*, Parameter*, RCB_Method, int, double);
void staggeredBalance(Grid*, Atom*, Parameter*, double);
// rcb methods
MD_FLOAT meanBisect(Atom*, MPI_Comm, int, double);
MD_FLOAT meanTimeBisect(Atom*, MPI_Comm, int, double);
#endif
#endif
