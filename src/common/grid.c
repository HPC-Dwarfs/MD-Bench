#include <allocate.h>
#include <grid.h>
#include <math.h>
#include <parameter.h>
#include <stdio.h>
#include <util.h>
#include <string.h>

#ifdef _MPI
MD_FLOAT f_normalization(double* x, MD_FLOAT* fx, MD_FLOAT minx, int nprocs) {
    MD_FLOAT sum = 0;
    for (int n = 0; n < nprocs; n++) {
        fx[n] = MAX(minx, x[n]);
        sum += fx[n];
    }
    for (int n = 0; n < nprocs; n++)
        fx[n] /= sum;
}

void fixedPointIteration(double* x0, int nprocs, MD_FLOAT minx) {
    MD_FLOAT tolerance = 1e-3;
    MD_FLOAT alpha     = 0.5;
    MD_FLOAT* fx       = (MD_FLOAT*)malloc(nprocs * sizeof(MD_FLOAT));
    int maxIterations  = 100;

    for (int i = 0; i < maxIterations; i++) {
        int converged = 1;
        f_normalization(x0, fx, minx, nprocs);

        for (int n = 0; n < nprocs; n++)
            fx[n] = (1 - alpha) * x0[n] + alpha * fx[n];

        for (int n = 0; n < nprocs; n++) {
            if (fabs(fx[n] - x0[n]) >= tolerance) {
                converged = 0;
                break;
            }
        }

        for (int n = 0; n < nprocs; n++)
            x0[n] = fx[n];

        if (converged) {
            for (int n = 0; n < nprocs; n++)
                return;
        }
    }
}

void staggeredBalance(Grid* grid, Atom* atom, Parameter* param, double newTime) {
    int me;
    int* coord  = grid->coord;
    int* nprocs = grid->nprocs;
    // Elapsed time since the last rebalance
    double time = newTime - grid->Timer;

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    grid->Timer = newTime;

    // store the older dimm to compare later for exchange
    MD_FLOAT lo[3], hi[3];
    for (int dim = 0; dim < 3; dim++) {
        lo[dim] = atom->mybox.lo[dim];
        hi[dim] = atom->mybox.hi[dim];
    }

    // Define parameters
    MPI_Comm subComm[3];
    int color[3]    = { 0, 0, 0 };
    int id[3]       = { 0, 0, 0 }; 
    double** load = (double**) malloc(3 * sizeof(double*));

    for(int dim = 0; dim < 3; dim++) {
        load[dim] = (double*)malloc(nprocs[dim] * sizeof(double));
    }

    int maxprocs       = MAX(MAX(nprocs[0], nprocs[1]), nprocs[2]);
    MD_FLOAT* cellSize = (MD_FLOAT*)malloc(maxprocs * sizeof(MD_FLOAT));
    MD_FLOAT* limits   = (MD_FLOAT*)malloc(
        2 * maxprocs * sizeof(MD_FLOAT)); // limits: (x0, x1), (x1, x2)... Repeat values
                                            // in between to perfom MPI_Scatter later
    MD_FLOAT t_sum[3]        = { 0, 0, 0 };
    MD_FLOAT recv_buf[2]     = { 0, 0 }; // Each proc only receives 2 elments per dimension xlo and xhi
    MD_FLOAT balancedLoad[3] = { 0, 0, 0 }; // 1/nprocs
    MD_FLOAT minLoad[3]      = { 0, 0, 0 }; // beta*(1/nprocs)
    MD_FLOAT prd[3]          = { param->xprd, param->yprd, param->zprd };
    MD_FLOAT boundaries[6]   = { 0, 0, 0, 0, 0, 0 }; // xlo,xhi,ylo,yhi,zlo,zhi

    // Create sub-communications along each dimension
    for (int dim = 0; dim < 3; dim++) {
        if (dim == 0) {
            color[0] = (coord[1] == 0 && coord[2] == 0) ? 1 : MPI_UNDEFINED;
            id[0]    = me;
        } else if (dim == 1) {
            color[1] = coord[2] == 0 ? coord[0] : MPI_UNDEFINED;
            id[1]    = (coord[1] == 0 && coord[2] == 0) ? 0 : me;
        } else {
            color[2] = coord[1] * nprocs[0] + coord[0];
            id[2]    = coord[2] == 0 ? 0 : me;
        }

        MPI_Comm_split(world, color[dim], id[dim], &subComm[dim]);
    }

    // Set the minimum load and the balance load
    for (int dim = 0; dim < 3; dim++) {
        balancedLoad[dim] = 1. / nprocs[dim];
        minLoad[dim]      = 0.8 * balancedLoad[dim];
    }
    // set and communicate the workload in reverse order
    for (int dim = 2; dim >= 0; dim--) {
        if (subComm[dim] != MPI_COMM_NULL) {
            MPI_Gather(&time, 1, MPI_DOUBLE, load[dim], 1, MPI_DOUBLE, 0, subComm[dim]);

            if (id[dim] == 0) {
                for (int n = 0; n < nprocs[dim]; n++)
                    t_sum[dim] += load[dim][n];

                for (int n = 0; n < nprocs[dim]; n++)
                    load[dim][n] /= t_sum[dim];
            }
            time = t_sum[dim];
        }

        MPI_Barrier(world);
    }

    // Brodacast the new boundaries along dimensions
    for (int dim = 0; dim < 3; dim++) {
        if (subComm[dim] != MPI_COMM_NULL) {
            MPI_Bcast(boundaries, 6, type_float, 0, subComm[dim]);

            if (id[dim] == 0) {
                MD_FLOAT inv_sum = 0.0;
                MD_FLOAT sum = 0.0;

                fixedPointIteration(load[dim], nprocs[dim], minLoad[dim]);                

                for(int n = 0; n < nprocs[dim]; n++) {
                    inv_sum += (1 / load[dim][n]);
                }

                for(int n = 0; n < nprocs[dim]; n++) {
                    cellSize[n] = (prd[dim] / load[dim][n]) * (1. / inv_sum);
                }

                for(int n = 0; n < nprocs[dim]; n++) {
                    limits[2 * n]     = sum;
                    limits[2 * n + 1] = sum + cellSize[n];
                    sum += cellSize[n];
                }

                limits[2 * nprocs[dim] - 1] = prd[dim];
            }

            MPI_Scatter(limits, 2, type_float, recv_buf, 2, type_float, 0, subComm[dim]);
            boundaries[2 * dim]     = recv_buf[0];
            boundaries[2 * dim + 1] = recv_buf[1];
        }

        MPI_Barrier(world);
    }
 
    atom->mybox.lo[0] = boundaries[0];
    atom->mybox.hi[0] = boundaries[1];
    atom->mybox.lo[1] = boundaries[2];
    atom->mybox.hi[1] = boundaries[3];
    atom->mybox.lo[2] = boundaries[4];
    atom->mybox.hi[2] = boundaries[5];

    MD_FLOAT domain[6] = {
        boundaries[0],
        boundaries[2],
        boundaries[4],
        boundaries[1],
        boundaries[3],
        boundaries[5]
    };

    MPI_Allgather(domain, 6, type_float, grid->map, 6, type_float, world);

    // because cells change dynamically, It is required to increase the
    // neighbouring exchange region
    for (int dim = 0; dim <= 2; dim++) {
        MD_FLOAT dr, dr_max;
        int n             = grid->nprocs[dim];
        MD_FLOAT maxdelta = 0.2 * prd[dim];

        dr                = MAX(
            fabs(lo[dim] - atom->mybox.lo[dim]),
            fabs(hi[dim] - atom->mybox.hi[dim]));

        MPI_Allreduce(&dr, &dr_max, 1, type_float, MPI_MAX, world);
        grid->cutneigh[dim] = param->cutneigh + dr_max;
    }

    for (int dim = 0; dim < 3; dim++) {
        if (subComm[dim] != MPI_COMM_NULL) {
            MPI_Comm_free(&subComm[dim]);
        }

        free(load[dim]);
    }

    free(load);
    free(limits);
}

// RCB Balancing
MD_FLOAT meanTimeBisect(Atom* atom, MPI_Comm subComm, int dim, double time) {
    MD_FLOAT mean = 0, sum = 0, total_sum = 0, weightAtoms = 0, total_weight = 0;

    for (int i = 0; i < atom->Nlocal; i++) {
        sum += atom_pos(i);
    }

    sum *= time;
    weightAtoms = atom->Nlocal * time;

    MPI_Allreduce(&sum, &total_sum, 1, type_float, MPI_SUM, subComm);
    MPI_Allreduce(&weightAtoms, &total_weight, 1, type_float, MPI_SUM, subComm);
    mean = total_sum / total_weight;
    return mean;
}

MD_FLOAT meanBisect(Atom* atom, MPI_Comm subComm, int dim, double time) {    
    int Natoms   = 0;
    MD_FLOAT sum = 0, mean = 0, total_sum = 0;

    for (int i = 0; i < atom->Nlocal; i++) {
        sum += atom_pos(i);
    }

    MPI_Allreduce(&sum, &total_sum, 1, type_float, MPI_SUM, subComm);
    MPI_Allreduce(&atom->Nlocal, &Natoms, 1, MPI_INT, MPI_SUM, subComm);
    mean = total_sum / Natoms;
    return mean;
}

void nextBisectionLevel(Grid* grid,
    Atom* atom,
    RCB_Method method,
    MPI_Comm subComm,
    int dim,
    int* color,
    int ilevel,
    double time) 
{
    int me;
    MPI_Comm_rank(world, &me);
    
    int rank, size;
    int branch = 0, i = 0, m = 0;
    int nsend = 0, nrecv = 0, nrecv2 = 0;
    int values_per_atom = 7;
    MD_FLOAT bisection, pos;
    MPI_Request request[2] = { MPI_REQUEST_NULL, MPI_REQUEST_NULL };
    MPI_Comm_rank(subComm, &rank);
    MPI_Comm_size(subComm, &size);


    int odd       = size % 2;
    int extraProc = odd ? size - 1 : size;
    int half      = (int)(0.5 * size);
    int partner   = (rank < half) ? rank + half : rank - half;
    

    if (odd && rank == extraProc) {
        partner = 0;
    }

    // Apply the bisection
    bisection = method(atom, subComm, dim, time);

    // Define the new boundaries
    if (rank < half) {
        atom->mybox.hi[dim] = bisection;
        branch              = 0;
    } else {
        atom->mybox.lo[dim] = bisection;
        branch              = 1;
    }

    // Define new color for the further communicaton
    *color = (branch << ilevel) | *color;

    // Grow the send buffer
    if (atom->Nlocal >= grid->maxsend) {
        if (grid->buf_send) free(grid->buf_send);
        grid->buf_send = (MD_FLOAT*)malloc(
            atom->Nlocal * values_per_atom * sizeof(MD_FLOAT));
        grid->maxsend = atom->Nlocal;
    }

    // buffer particles to send
    while (i < atom->Nlocal) {
        pos = atom_pos(i);
        if (pos < atom->mybox.lo[dim] || pos >= atom->mybox.hi[dim]) {
            nsend += packExchange(atom, i, &grid->buf_send[nsend]);
            copy(atom, i, atom->Nlocal - 1);
            atom->Nlocal--;
        } else {
            i++;
        }
    }

    // Communicate the number of elements to be sent
    if (rank < extraProc) {
        MPI_Irecv(&nrecv, 1, MPI_INT, partner, 0, subComm, &request[0]);
    }
    if (odd && rank == 0) {
        MPI_Irecv(&nrecv2, 1, MPI_INT, extraProc, 0, subComm, &request[1]);
    }
    MPI_Send(&nsend, 1, MPI_INT, partner, 0, subComm);
    MPI_Waitall(2, request, MPI_STATUS_IGNORE);

    // Grow the recv buffer
    if (nrecv + nrecv2 >= grid->maxrecv) {
        if (grid->buf_recv) free(grid->buf_recv);
        grid->buf_recv = (MD_FLOAT*)malloc(
            (nrecv + nrecv2) * values_per_atom * sizeof(MD_FLOAT));
        grid->maxrecv = nrecv + nrecv2;
    }

    // communicate elements in the buffer
    request[0] = MPI_REQUEST_NULL;
    request[1] = MPI_REQUEST_NULL;

    if (rank < extraProc) {
        MPI_Irecv(grid->buf_recv, nrecv, type_float, partner, 0, subComm, &request[0]);
    }

    if (odd && rank == 0) {
        MPI_Irecv(&grid->buf_recv[nrecv],
            nrecv2,
            type_float,
            extraProc,
            0,
            subComm,
            &request[1]);
    }

    MPI_Send(grid->buf_send, nsend, type_float, partner, 0, subComm);
    MPI_Waitall(2, request, MPI_STATUS_IGNORE);

    // store atoms in atom list
    while (m < nrecv + nrecv2) {
        m += unpackExchange(atom, atom->Nlocal++, &grid->buf_recv[m]);
    }
}

void rcbBalance(Grid* grid, Atom* atom, Parameter* param, RCB_Method method, int ndim, double newTime) {
    int me, nprocs = 0, ilevel = 0, nboxes = 1;
    int color = 0, size = 0;
    int index, prd[3];
    MPI_Comm subComm;
    MPI_Comm_size(world, &nprocs);
    MPI_Comm_rank(world, &me);

    // set the elapsed time since the last dynamic balance
    double time = newTime - grid->Timer;

    prd[0] = atom->mybox.xprd = param->xprd;
    prd[1] = atom->mybox.yprd = param->yprd;
    prd[2] = atom->mybox.zprd = param->zprd;

    // Sort by larger dimension
    int largerDim[3] = { 0, 1, 2 };

    for (int i = 0; i < 2; i++) {
        for (int j = i + 1; j < 3; j++) {
            if (prd[largerDim[j]] > prd[largerDim[i]]) {
                MD_FLOAT tmp = largerDim[j];
                largerDim[j] = largerDim[i];
                largerDim[i] = tmp;
            }
        }
    }

    // Initial Partition
    atom->mybox.lo[0] = 0;
    atom->mybox.hi[0] = atom->mybox.xprd;
    atom->mybox.lo[1] = 0;
    atom->mybox.hi[1] = atom->mybox.yprd;
    atom->mybox.lo[2] = 0;
    atom->mybox.hi[2] = atom->mybox.zprd;

    // Recursion tree
    while (nboxes < nprocs) {
        index = ilevel % ndim;
        MPI_Comm_split(world, color, me, &subComm);
        MPI_Comm_size(subComm, &size);

        if (size > 1) {
            nextBisectionLevel(grid,
                atom,
                method,
                subComm,
                largerDim[index],
                &color,
                ilevel,
                time);
        }

        MPI_Comm_free(&subComm);
        nboxes = pow(2, ++ilevel);
    }

    // Set the new timer grid
    grid->Timer = newTime;

    // Creating the global map
    MD_FLOAT domain[6] = { atom->mybox.lo[0],
        atom->mybox.lo[1],
        atom->mybox.lo[2],
        atom->mybox.hi[0],
        atom->mybox.hi[1],
        atom->mybox.hi[2]
    };

    MPI_Allgather(domain, 6, type_float, grid->map, 6, type_float, world);

    // Define the same cutneighbour in all dimensions for the exchange
    // communication
    for (int dim = 0; dim <= 2; dim++){
        grid->cutneigh[dim] = param->cutneigh;
    }
}

// Regular grid
void init3DCartesian(Grid* grid, Parameter* param, Box* box) {
    int me, nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    MPI_Comm cartesian;
    int numdim     = 3;
    int reorder    = 0;
    int periods[3] = { 1, 1, 1 };
    int mycoord[3] = { 0, 0, 0 };
    int griddim[3] = { 0, 0, 0 };
    MD_FLOAT len[3];
    MD_FLOAT eps = 1e-3; 

    box->xprd = param->xprd;
    box->yprd = param->yprd;
    box->zprd = param->zprd;

    // Creates a cartesian 3d grid
    MPI_Dims_create(nproc, numdim, griddim);
    MPI_Cart_create(world, numdim, griddim, periods, reorder, &cartesian);
    grid->nprocs[0] = griddim[0];
    grid->nprocs[1] = griddim[1];
    grid->nprocs[2] = griddim[2];

    // Coordinates position in the grid
    MPI_Cart_coords(cartesian, me, 3, mycoord);
    grid->coord[0] = mycoord[0];
    grid->coord[1] = mycoord[1];
    grid->coord[2] = mycoord[2];

    // Boundaries of my local box, with origin in (0,0,0).
    len[0] = param->xprd / griddim[0];
    len[1] = param->yprd / griddim[1];
    len[2] = param->zprd / griddim[2];

    box->lo[0] = mycoord[0] * len[0];
    box->hi[0] = (mycoord[0] + 1) * len[0];
    box->lo[1] = mycoord[1] * len[1];
    box->hi[1] = (mycoord[1] + 1) * len[1];
    box->lo[2] = mycoord[2] * len[2];
    box->hi[2] = (mycoord[2] + 1) * len[2];

    if (box->hi[0] + eps > param->xprd) {
        box->hi[0] = param->xprd;
    }

    if (box->hi[1] + eps > param->yprd) {
        box->hi[1] = param->yprd;
    }

    if (box->hi[2] + eps > param->zprd) {
        box->hi[2] = param->zprd;
    }

    MD_FLOAT domain[6] = {
        box->lo[0],
        box->lo[1],
        box->lo[2],
        box->hi[0],
        box->hi[1],
        box->hi[2]
    };

    MPI_Allgather(domain, 6, type_float, grid->map, 6, type_float, world);
    MPI_Comm_free(&cartesian);
}

// Other Functions from the grid
void initGrid(Grid* grid, int nprocs) {
    grid->map_size = 6 * nprocs;
    grid->map      = (MD_FLOAT*)allocate(ALIGNMENT, grid->map_size * sizeof(MD_FLOAT));
    // RCB
    grid->maxsend  = 0;
    grid->maxrecv  = 0;
    grid->buf_send = NULL;
    grid->buf_recv = NULL;
    // Staggered
    grid->Timer = 0.;
}

int readAtomsTempFile(Atom* atom, char* file) {
    FILE *fp = fopen(file, "r");
    if (fp == NULL) {
        perror("Error opening file");
        return -1;
    }

    if(atom->Nmax > 0){
        freeAtom(atom);
        atom->Nmax = 0;
    } 

    MD_FLOAT xlo = atom->mybox.lo[0];
    MD_FLOAT xhi = atom->mybox.hi[0];
    MD_FLOAT ylo = atom->mybox.lo[1];
    MD_FLOAT yhi = atom->mybox.hi[1];
    MD_FLOAT zlo = atom->mybox.lo[2];
    MD_FLOAT zhi = atom->mybox.hi[2];

    MD_FLOAT x, y, z, vx, vy, vz;
    int type;
    int i = 0;

#if PRECISION == 1
    #define SCANF_FORMAT "%f %f %f %f %f %f %d"
#else
    #define SCANF_FORMAT "%lf %lf %lf %lf %lf %lf %d"
#endif

    while (fscanf(fp, SCANF_FORMAT, &x, &y, &z, &vx, &vy, &vz, &type) == 7) {
        if (x >= xlo && x < xhi && y >= ylo && y < yhi && z >= zlo && z < zhi ) {
                if (i == atom->Nmax) {
                    growAtom(atom);
                }
                
                atom_x(i) = x;
                atom_y(i) = y;
                atom_z(i) = z;
                atom_vx(i) = vx;
                atom_vy(i) = vy;
                atom_vz(i) = vz;
                atom->type[i] = type;
                i++;
            }
    }
 
    fclose(fp);
    atom->Nlocal = i;
    return 0;
}

void discardAtomsOutsideSubdomainBox(Atom* atom) {
    MD_FLOAT xlo = atom->mybox.lo[0];
    MD_FLOAT xhi = atom->mybox.hi[0];
    MD_FLOAT ylo = atom->mybox.lo[1];
    MD_FLOAT yhi = atom->mybox.hi[1];
    MD_FLOAT zlo = atom->mybox.lo[2];
    MD_FLOAT zhi = atom->mybox.hi[2];
    int n = 0;

    for(int i = 0; i < atom->Nlocal; i++) {
        MD_FLOAT x = atom_x(i);
        MD_FLOAT y = atom_y(i);
        MD_FLOAT z = atom_z(i);

        if (x >= xlo && x < xhi && y >= ylo && y < yhi && z >= zlo && z < zhi) {
            atom_x(n) = x;
            atom_y(n) = y;
            atom_z(n) = z;
            atom_vx(n) = atom_vx(i);
            atom_vy(n) = atom_vy(i);
            atom_vz(n) = atom_vz(i);
            atom->type[n] = atom->type[i];
            n++;
        }
    }

    atom->Nlocal = n;
}

void setupGrid(Grid* grid, Atom* atom, Parameter* param) {
    int me = 0;
    int nprocs = 1;
    MD_FLOAT xlo, ylo, zlo, xhi, yhi, zhi;

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(world, &nprocs);
    initGrid(grid, nprocs);

    // Set the origin at (0,0,0)
    if(param->input_file != NULL) {
        for (int i = 0; i < atom->Nlocal; i++) {
            atom_x(i) = atom_x(i) - param->xlo;
            atom_y(i) = atom_y(i) - param->ylo;
            atom_z(i) = atom_z(i) - param->zlo;
        }
    }

    // MAP is stored as follows: xlo,ylo,zlo,xhi,yhi,zhi
    init3DCartesian(grid, param, &atom->mybox);
    
    // Define the same cutneighbour in all dimensions for the exchange
    // communication
    for (int dim = 0; dim <= 2; dim++){
        grid->cutneigh[dim] = param->cutneigh;
    }

    MPI_Barrier(world);
    discardAtomsOutsideSubdomainBox(atom);

    //printGrid(grid);
    if (!param->balance) {
        MPI_Allreduce(&atom->Nlocal, &atom->Natoms, 1, MPI_INT, MPI_SUM, world);
        //fprintf(stdout,"Processor:%i, Local atoms:%i, Total atoms:%i\n",
        //    me,
        //    atom->Nlocal,
        //   atom->Natoms);
        //fflush(stdout);
        MPI_Barrier(world);
    }
}

void printGrid(Grid* grid) {
    int me = 0; 
    int nprocs = 1;
    MPI_Comm_size(world, &nprocs);
    MPI_Comm_rank(world, &me);
    MD_FLOAT* map = grid->map;

    if (me == 0) {
        printf("GRID:\n");
        printf("==================================================================="
               "================================\n");
        for (int i = 0; i < nprocs; i++)
            printf("Box:%i\txlo:%.4f\txhi:%.4f\tylo:%.4f\tyhi:%.4f\tzlo:%.4f\tzhi:%.4f\n",
                i,
                map[6 * i],
                map[6 * i + 3],
                map[6 * i + 1],
                map[6 * i + 4],
                map[6 * i + 2],
                map[6 * i + 5]);

        printf("\n\n");
        // printf("Box processor:%i\n xlo:%.4f\txhi:%.4f\n ylo:%.4f\tyhi:%.4f\n
        // zlo:%.4f\tzhi:%.4f\n",
        // i,map[6*i],map[6*i+3],map[6*i+1],map[6*i+4],map[6*i+2],map[6*i+5]);
    }

    MPI_Barrier(world);
}
#endif
