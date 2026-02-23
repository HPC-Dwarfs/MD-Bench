#include <atom.h>
#include <comm.h>
#include <grid.h>
#include <neighbor.h>
#include <parameter.h>
#include <timing.h>

#ifdef _MPI
double dynamicBalance(Comm* comm, Grid* grid, Atom* atom, Parameter* param, double time)
{
    double S, E;
    int dims = 3; // TODO: Adjust to do in 3d and 2d
    S        = getTimeStamp();
    if (param->balance == RCB) {
        rcbBalance(grid, atom, param, meanBisect, dims, 0);
        neighComm(comm, param, grid);
    } else if (param->balance == meanTimeRCB) {
        rcbBalance(grid, atom, param, meanTimeBisect, dims, time);
        neighComm(comm, param, grid);
    } else if (param->balance == Staggered) {
        staggeredBalance(grid, atom, param, time);
        neighComm(comm, param, grid);
        exchangeComm(comm, atom);
    } else {
        // Do nothing
    }
    // printGrid(grid);
    // MPI_Barrier(world);
    E = getTimeStamp();
    return E - S;
}

double initialBalance(Parameter* param,
    Atom* atom,
    Neighbor* neighbor,
    Stats* stats,
    Comm* comm,
    Grid* grid)
{
    double E, S, time;
    S = getTimeStamp();
    if (param->balance == meanTimeRCB || param->balance == RCB) {
        rcbBalance(grid, atom, param, meanBisect, 3, 0);
        neighComm(comm, param, grid);
    }
    MPI_Allreduce(&atom->Nlocal, &atom->Natoms, 1, MPI_INT, MPI_SUM, world);
    // printf("Processor:%i, Local atoms:%i, Total atoms:%i\n",
    //     comm->myproc,
    //     atom->Nlocal,
    //     atom->Natoms);
    MPI_Barrier(world);
    E = getTimeStamp();
    return E - S;
}
#else
double dynamicBalance(Comm* comm, Grid* grid, Atom* atom, Parameter* param, double time)
{
    return 0;
}
double initialBalance(Parameter* param,
    Atom* atom,
    Neighbor* neighbor,
    Stats* stats,
    Comm* comm,
    Grid* grid)
{
    return 0;
}
#endif
