/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
//---
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
//---
#include <likwid-marker.h>

#include <device.h>

extern "C" {
#include <allocate.h>
#include <atom.h>
#include <comm.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <timing.h>
#include <util.h>
}

__global__ void computeForceLJCudaFullNeigh(DeviceAtom a,
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    int Nlocal,
    int neigh_maxneighs,
    int* neigh_neighbors,
    int* neigh_numneigh,
    int ntypes)
{

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom    = &a;
    const int numneighs = neigh_numneigh[i];

    MD_FLOAT xtmp = atom_x(i);
    MD_FLOAT ytmp = atom_y(i);
    MD_FLOAT ztmp = atom_z(i);

    MD_FLOAT fix = 0;
    MD_FLOAT fiy = 0;
    MD_FLOAT fiz = 0;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
    const int type_i = atom->type[i];
#endif

    for (int k = 0; k < numneighs; k++) {
        int j         = neighs(neigh_neighbors, i, k, Nlocal, neigh_maxneighs);
        MD_FLOAT delx = xtmp - atom_x(j);
        MD_FLOAT dely = ytmp - atom_y(j);
        MD_FLOAT delz = ztmp - atom_z(j);
        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
        const int type_j          = atom->type[j];
        const int type_ij         = type_i * ntypes + type_j;
        const MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
        const MD_FLOAT sigma6     = atom->sigma6[type_ij];
        const MD_FLOAT epsilon    = atom->epsilon[type_ij];
#endif

        if (rsq < cutforcesq) {
            MD_FLOAT sr2   = (MD_FLOAT)1.0 / rsq;
            MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
            MD_FLOAT force = (MD_FLOAT)48.0 * sr6 * (sr6 - (MD_FLOAT)0.5) * sr2 * epsilon;
            fix += delx * force;
            fiy += dely * force;
            fiz += delz * force;
        }
    }

    atom_fx(i) = fix;
    atom_fy(i) = fiy;
    atom_fz(i) = fiz;
}

__global__ void computeForceLJCudaHalfNeigh(DeviceAtom a,
    MD_FLOAT cutforcesq,
    MD_FLOAT sigma6,
    MD_FLOAT epsilon,
    int Nlocal,
    int neigh_maxneighs,
    int* neigh_neighbors,
    int* neigh_numneigh,
    int ntypes)
{

    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom    = &a;
    const int numneighs = neigh_numneigh[i];

    MD_FLOAT xtmp = atom_x(i);
    MD_FLOAT ytmp = atom_y(i);
    MD_FLOAT ztmp = atom_z(i);

    MD_FLOAT fix = 0;
    MD_FLOAT fiy = 0;
    MD_FLOAT fiz = 0;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
    const int type_i = atom->type[i];
#endif

    for (int k = 0; k < numneighs; k++) {
        int j         = neighs(neigh_neighbors, i, k, Nlocal, neigh_maxneighs);
        MD_FLOAT delx = xtmp - atom_x(j);
        MD_FLOAT dely = ytmp - atom_y(j);
        MD_FLOAT delz = ztmp - atom_z(j);
        MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#if LJ_COMB_RULE != LJ_COMB_SINGLE
        const int type_j          = atom->type[j];
        const int type_ij         = type_i * ntypes + type_j;
        const MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
        const MD_FLOAT sigma6     = atom->sigma6[type_ij];
        const MD_FLOAT epsilon    = atom->epsilon[type_ij];
#endif

        if (rsq < cutforcesq) {
            MD_FLOAT sr2   = (MD_FLOAT)1.0 / rsq;
            MD_FLOAT sr6   = sr2 * sr2 * sr2 * sigma6;
            MD_FLOAT force = (MD_FLOAT)48.0 * sr6 * (sr6 - (MD_FLOAT)0.5) * sr2 * epsilon;
            MD_FLOAT partial_force_x = delx * force;
            MD_FLOAT partial_force_y = dely * force;
            MD_FLOAT partial_force_z = delz * force;

            atomicAdd(&atom_fx(j), -partial_force_x);
            atomicAdd(&atom_fy(j), -partial_force_y);
            atomicAdd(&atom_fz(j), -partial_force_z);

            fix += partial_force_x;
            fiy += partial_force_y;
            fiz += partial_force_z;
        }
    }

    atomicAdd(&atom_fx(i), fix);
    atomicAdd(&atom_fy(i), fiy);
    atomicAdd(&atom_fz(i), fiz);
}

__global__ void kernel_initial_integrate(
    MD_FLOAT dtforce, MD_FLOAT dt, int Nlocal, DeviceAtom a)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom = &a;

    atom_vx(i) += dtforce * atom_fx(i);
    atom_vy(i) += dtforce * atom_fy(i);
    atom_vz(i) += dtforce * atom_fz(i);
    atom_x(i) = atom_x(i) + dt * atom_vx(i);
    atom_y(i) = atom_y(i) + dt * atom_vy(i);
    atom_z(i) = atom_z(i) + dt * atom_vz(i);
}

__global__ void kernel_final_integrate(MD_FLOAT dtforce, int Nlocal, DeviceAtom a)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= Nlocal) {
        return;
    }

    DeviceAtom* atom = &a;

    atom_vx(i) += dtforce * atom_fx(i);
    atom_vy(i) += dtforce * atom_fy(i);
    atom_vz(i) += dtforce * atom_fz(i);
}

extern "C" {

void finalIntegrateCUDA(bool reneigh, Parameter* param, Atom* atom)
{
    DEBUG_MESSAGE("finalIntegrateCUDA begin\n");

    const int Nlocal                = atom->Nlocal;
    const int num_threads_per_block = get_cuda_num_threads();
    const int num_blocks            = ceil((float)Nlocal / (float)num_threads_per_block);

    kernel_final_integrate<<<num_blocks, num_threads_per_block>>>(param->dtforce,
        Nlocal,
        atom->d_atom);
    cuda_assert("kernel_final_integrate", cudaPeekAtLastError());
    cuda_assert("kernel_final_integrate", cudaDeviceSynchronize());

    if (reneigh) {
        memcpyFromGPU(atom->vx, atom->d_atom.vx, sizeof(MD_FLOAT) * atom->Nlocal * 3);
    }

    DEBUG_MESSAGE("finalIntegrateCUDA end\n");
}

void initialIntegrateCUDA(bool reneigh, Parameter* param, Atom* atom)
{
    DEBUG_MESSAGE("initialIntegrateCUDA begin\n");

    const int Nlocal                = atom->Nlocal;
    const int num_threads_per_block = get_cuda_num_threads();
    const int num_blocks            = ceil((float)Nlocal / (float)num_threads_per_block);

    kernel_initial_integrate<<<num_blocks, num_threads_per_block>>>(param->dtforce,
        param->dt,
        Nlocal,
        atom->d_atom);
    cuda_assert("kernel_initial_integrate", cudaPeekAtLastError());
    cuda_assert("kernel_initial_integrate", cudaDeviceSynchronize());

    if (reneigh) {
        memcpyFromGPU(atom->vx, atom->d_atom.vx, sizeof(MD_FLOAT) * atom->Nlocal * 3);
    }

    DEBUG_MESSAGE("initialIntegrateCUDA end\n");
}

double computeForceLJCUDA(Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    DEBUG_MESSAGE("computeForceLJCUDA begin\n");

    const int num_threads_per_block = get_cuda_num_threads();
    int Nlocal                      = atom->Nlocal;
    int Nmax                        = atom->Nmax;
    MD_FLOAT cutforcesq             = param->cutforce * param->cutforce;
    MD_FLOAT sigma6                 = param->sigma6;
    MD_FLOAT epsilon                = param->epsilon;

    /*
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    size_t free, total;
    for(int i = 0; i < nDevices; ++i) {
        cudaMemGetInfo( &free, &total );
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("DEVICE %d/%d NAME: %s\r\n with %ld MB/%ld MB memory used", i + 1,
    nDevices, prop.name, free / 1024 / 1024, total / 1024 / 1024);
    }
    */

    // HINT: Run with cuda-memcheck ./MDBench-NVCC in case of error
    // memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nlocal * 3);

    cudaProfilerStart();
    const int num_blocks = ceil((float)Nlocal / (float)num_threads_per_block);
    double S             = getTimeStamp();
    LIKWID_MARKER_START("force");

    if (neighbor->half_neigh) {
#ifdef ATOM_POSITION_AOS
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax * 3);
#else
        memsetGPU(atom->d_atom.fx, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fy, 0, sizeof(MD_FLOAT) * Nmax);
        memsetGPU(atom->d_atom.fz, 0, sizeof(MD_FLOAT) * Nmax);
#endif
        computeForceLJCudaHalfNeigh<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
            cutforcesq,
            sigma6,
            epsilon,
            Nlocal,
            neighbor->maxneighs,
            neighbor->d_neighbor.neighbors,
            neighbor->d_neighbor.numneigh,
            atom->ntypes);
    } else {
        computeForceLJCudaFullNeigh<<<num_blocks, num_threads_per_block>>>(atom->d_atom,
            cutforcesq,
            sigma6,
            epsilon,
            Nlocal,
            neighbor->maxneighs,
            neighbor->d_neighbor.neighbors,
            neighbor->d_neighbor.numneigh,
            atom->ntypes);
    }

    cuda_assert("computeForceLJCuda", cudaPeekAtLastError());
    cuda_assert("computeForceLJCuda", cudaDeviceSynchronize());
    cudaProfilerStop();

    LIKWID_MARKER_STOP("force");
    double E = getTimeStamp();
    DEBUG_MESSAGE("computeForceLJCUDA end\n");
    return E - S;
}
}

extern "C" void copyGhostFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->x, atom->d_atom.x, atom->Nlocal * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyGhostToGPU(Atom* atom)
{
    memcpyToGPU(&atom->d_atom.x[atom->Nlocal * 3],
        &atom->x[atom->Nlocal * 3],
        atom->Nghost * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyForceFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->fx, atom->d_atom.fx, atom->Nmax * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyForceToGPU(Atom* atom)
{
    memcpyToGPU(atom->d_atom.fx, atom->fx, atom->Nmax * sizeof(MD_FLOAT) * 3);
}

extern "C" void copyDataFromCUDADevice(Parameter* param, Atom* atom)
{
    memcpyFromGPU(atom->x, atom->d_atom.x, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyFromGPU(atom->vx, atom->d_atom.vx, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyFromGPU(atom->type, atom->d_atom.type, atom->Nmax * sizeof(int));
}

extern "C" void copyDataToCUDADevice(Parameter* param, Atom* atom)
{
    memcpyToGPU(atom->d_atom.x, atom->x, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyToGPU(atom->d_atom.vx, atom->vx, atom->Nmax * sizeof(MD_FLOAT) * 3);
    memcpyToGPU(atom->d_atom.type, atom->type, atom->Nmax * sizeof(int));
}

/*
__global__ void pack_forward_Cuda(DeviceAtom a,
    int n,
    int* cuda_list,
    MD_FLOAT* buf,
    int PBCx,
    int PBCy,
    int PBCz,
    MD_FLOAT xprd,
    MD_FLOAT yprd,
    MD_FLOAT zprd)
{
    DeviceAtom* atom    = &a;
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    int j        = cuda_list[i];
    buf_x(i) = atom_x(j) + PBCx * xprd;
    buf_y(i) = atom_y(j) + PBCy * yprd;
    buf_z(i) = atom_z(j) + PBCz * zprd;
}

__global__ void unpack_forward_Cuda(DeviceAtom a,
    int n,
    int first,
    MD_FLOAT* buf)
{
    DeviceAtom* atom    = &a;
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    atom_x((first + i)) = buf_x(i);
    atom_y((first + i)) = buf_y(i);
    atom_z((first + i)) = buf_z(i);
}

extern "C" void forwardCommCUDA(Comm* comm, Atom* atom, int iswap)
{
    int nrqst = 0, offset = 0, nsend = 0, nrecv = 0;
    int pbc[3];
    int size    = comm->forwardSize;
    int maxrqst = comm->numneigh;
    int* cuda_sendlist;
    int max_list_size;
    cuda_buf_send =  (MD_FLOAT*)allocateGPU(comm->maxsend * sizeof(MD_FLOAT));
    cuda_buf_recv =  (MD_FLOAT*)allocateGPU(comm->maxrecv * sizeof(MD_FLOAT));

    //use a single buffer and takes the highes list size to move list of cluster to send
    for (int ineigh = 0; ineigh < comm->numneigh; ineigh++){
        max_list_size = comm->maxsendlist[ineigh];
    }
    //allocate the memory for the unique buffer
    cuda_sendlist = (int*)allocateGPU(max_list_size * sizeof(int));

    for (int ineigh = comm->sendfrom[iswap]; ineigh < comm->sendtill[iswap]; ineigh++) {
        offset  = comm->off_atom_send[ineigh];
        pbc[_x] = comm->pbc_x[ineigh];
        pbc[_y] = comm->pbc_y[ineigh];
        pbc[_z] = comm->pbc_z[ineigh];
        //copy lists into the buffer
        memcpyToGPU(cuda_sendlist, comm->sendlist[ineigh], comm->atom_send[ineigh] *
sizeof(int));

        const int num_threads_per_block = get_cuda_num_threads();
        const int num_blocks = ceil((float) comm->atom_send[ineigh] /
(float)num_threads_per_block);

        pack_forward_Cuda<<<num_blocks, num_threads_per_block>>>(
                                            atom->d_atom,                   //MD_FLOAT*
-->need to be in tye device comm->atom_send[ineigh],        //int cuda_sendlist, //int*
-->need to be in tye device &cuda_buf_send[offset*size],    //MD_FLOAT*  -->need to be in
tye device pbc[_x],                        //int pbc[_y],                        //int
                                            pbc[_z],                        //int
                                            atom->mybox.xprd, //MD_FLOAT atom->mybox.yprd,
//MD_FLOAT atom->mybox.zprd);                    //MD_FLOAT cudaDeviceSynchronize();
        }

#ifdef _MPI
    MPI_Request requests[maxrqst];
    // Receives elements
    if (comm->othersend[iswap])
        for (int ineigh = comm->recvfrom[iswap]; ineigh < comm->recvtill[iswap];ineigh++)
{ offset = comm->off_atom_recv[ineigh] * size; nrecv  = comm->atom_recv[ineigh] * size;
            MPI_Irecv(&cuda_buf_recv[offset],
                    nrecv,
                    type,
                    comm->nrecv[ineigh],
                    0,
                    world,
                    &requests[nrqst++]);
        }

    // Send elements
    if (comm->othersend[iswap])
        for (int ineigh = comm->sendfrom[iswap]; ineigh < comm->sendtill[iswap]; ineigh++)
{ offset = comm->off_atom_send[ineigh] * size; nsend  = comm->atom_send[ineigh] * size;
            MPI_Send(&cuda_buf_send[offset], nsend, type, comm->nsend[ineigh], 0, world);
        }

    if (comm->othersend[iswap]) MPI_Waitall(nrqst, requests, MPI_STATUS_IGNORE);
#endif

    // unpack buffer
    for (int ineigh = comm->recvfrom[iswap]; ineigh < comm->recvtill[iswap]; ineigh++) {
        offset = comm->off_atom_recv[ineigh];
        MD_FLOAT *buf = (comm->othersend[iswap]) ? cuda_buf_recv : cuda_buf_send;

        const int num_threads_per_block = get_cuda_num_threads();
        const int num_blocks = ceil((float) comm->atom_send[ineigh] /
(float)num_threads_per_block);

        unpack_forward_Cuda<<<num_blocks, num_threads_per_block>>>(
                                        atom->d_atom,                       //MD_FLOAT*
--> need to be in the device comm->atom_recv[ineigh],            //int
                                        comm->firstrecv[iswap] + offset,    //int
                                        &buf[offset * size]);               //MD_FLOAT*
--> need to be in the devic cudaDeviceSynchronize();
        }
    cuda_assert("cudaDeviceFree", cudaFree(cuda_sendlist));
    cuda_assert("cudaDeviceFree", cudaFree(cuda_buf_recv));
    cuda_assert("cudaDeviceFree", cudaFree(cuda_buf_send));
}
*/
