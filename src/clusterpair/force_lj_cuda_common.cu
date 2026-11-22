/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */

#include <device.h>

extern "C" {
#include <stdio.h>
//---
#include <likwid-marker.h>
//---
#include <allocate.h>
#include <atom.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>
}

#include <force_lj_cuda_common.h>

extern "C" {
MD_FLOAT* cuda_cl_x;
MD_FLOAT* cuda_cl_v;
MD_FLOAT* cuda_cl_f;
int* cuda_neighbors;
unsigned int* cuda_neighbors_imask;
int* cuda_numneigh;
int* cuda_numneigh_inner;
int* cuda_natoms;
int* natoms;
int* ngatoms;
int* cuda_border_map;
int* cuda_jclusters_natoms;
int *cuda_PBCx, *cuda_PBCy, *cuda_PBCz;

// LJ data of the combination rules: allocated for every rule (--lj-comb is a runtime
// option), copied only for the active one
int* cuda_cl_t;
MD_FLOAT* cuda_cl_sqrt_epsilon;
MD_FLOAT* cuda_cl_sigma3;
// Charges (only with electrostatics), molecules (copied when the system has molecules)
// and atom ids (negative for padding slots)
MD_FLOAT* cuda_cl_q;
int* cuda_cl_mol;
int* cuda_cl_id;

MD_FLOAT* cuda_cutforcesq;
MD_FLOAT* cuda_sigma6;
MD_FLOAT* cuda_epsilon;

MD_FLOAT4* cuda_cl_xq;
AtomParam* cuda_cl_par;
MD_FLOAT4* cuda_nbfp;
}

// Packs cuda_cl_xq and, with pack_par, cuda_cl_par for the combination rule rule
__global__ void packAtomsCuda(const MD_FLOAT* cuda_cl_x,
    const MD_FLOAT* cuda_cl_q,
    const MD_FLOAT* cuda_cl_sqrt_epsilon,
    const MD_FLOAT* cuda_cl_sigma3,
    const int* cuda_cl_mol,
    const int* cuda_cl_id,
    const int* cuda_cl_t,
    MD_FLOAT4* xq,
    AtomParam* par,
    int nslots,
    int pack_par,
    int rule,
    int excl)
{
    const int a = blockIdx.x * blockDim.x + threadIdx.x;
    if (a >= nslots) {
        return;
    }

    const MD_FLOAT* x = &cuda_cl_x[CJ_VECTOR_BASE_INDEX(a / CLUSTER_N)];
    const int i       = a % CLUSTER_N;
    MD_FLOAT4 v;
    v.x   = x[CL_X_INDEX(i)];
    v.y   = x[CL_Y_INDEX(i)];
    v.z   = x[CL_Z_INDEX(i)];
    v.w   = chargeCuda(cuda_cl_q, a);
    xq[a] = v;

    if (pack_par) {
        AtomParam p = {};
        if (rule == LJ_COMB_GEOM) {
            p.sqrt_epsilon = cuda_cl_sqrt_epsilon[a];
            p.sigma3       = cuda_cl_sigma3[a];
        } else if (rule == LJ_COMB_FULL) {
            p.type = cuda_cl_t[a];
        }
        p.mol  = cuda_cl_id[a] < 0 ? -1 : (excl ? cuda_cl_mol[a] : 0);
        par[a] = p;
    }
}

// Packs all local and ghost atom slots
extern "C" void launchPackAtomsCuda(Parameter* param, Atom* atom, int pack_par)
{
    const int nslots = (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_N;
    if (nslots == 0) {
        return;
    }
    packAtomsCuda<<<(nslots + 255) / 256, 256>>>(cuda_cl_x,
        cuda_cl_q,
        cuda_cl_sqrt_epsilon,
        cuda_cl_sigma3,
        cuda_cl_mol,
        cuda_cl_id,
        cuda_cl_t,
        cuda_cl_xq,
        cuda_cl_par,
        nslots,
        pack_par,
        param->lj_comb_rule,
        atom->molsize > 1);
    cuda_assert("packAtomsCuda", cudaPeekAtLastError());
}

// packs the arguments for the CUDA force kernel
NbForceArgs makeNbForceArgs(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    NbForceArgs nb;
    nb.xq              = cuda_cl_xq;
    nb.par             = cuda_cl_par;
    nb.f               = cuda_cl_f;
    nb.nclusters_local = atom->Nclusters_local;
    nb.numneigh        = cuda_numneigh_inner;
    nb.cj              = cuda_neighbors;
    nb.cj_mask         = NULL;
    nb.maxneighs       = neighbor->maxneighs;
    nb.cutforcesq      = param->cutforce * param->cutforce;
    nb.sigma6          = param->sigma6;
    nb.epsilon         = param->epsilon;
    nb.nbfp            = cuda_nbfp;
    nb.ntypes          = atom->ntypes;
    nb.cs              = coulombScalarInit(param);
    nb.energies        = NULL;
    return nb;
}

extern "C" void initDevice(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    cuda_cl_x = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * ATOM_DIM * sizeof(MD_FLOAT));
    cuda_cl_v = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
    cuda_cl_f = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
    cuda_cl_t = (int*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
    cuda_cl_sqrt_epsilon = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
    cuda_cl_sigma3 = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
    // cutforcesq, sigma6 and epsilon per type pair: separate arrays for the reference
    // kernels, packed {cutforcesq, sigma6, epsilon} for the optimized ones
    const int ntypes2 = atom->ntypes * atom->ntypes;
    cuda_cutforcesq   = (MD_FLOAT*)allocateGPU(ntypes2 * sizeof(MD_FLOAT));
    cuda_sigma6       = (MD_FLOAT*)allocateGPU(ntypes2 * sizeof(MD_FLOAT));
    cuda_epsilon      = (MD_FLOAT*)allocateGPU(ntypes2 * sizeof(MD_FLOAT));
    cuda_nbfp         = (MD_FLOAT4*)allocateGPU(ntypes2 * sizeof(MD_FLOAT4));
#if COULOMB_RULE != COULOMB_NONE
    cuda_cl_q = (MD_FLOAT*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
#endif
    cuda_cl_mol = (int*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
    cuda_cl_id = (int*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
    cuda_cl_xq = (MD_FLOAT4*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT4));
    cuda_cl_par = (AtomParam*)allocateGPU(
        atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(AtomParam));
    cuda_natoms = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_jclusters_natoms = (int*)allocateGPU(
        atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_border_map = (int*)allocateGPU(
        atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_PBCx     = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_PBCy     = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_PBCz     = (int*)allocateGPU(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    cuda_numneigh = (int*)allocateGPU(atom->Nclusters_max * sizeof(int));
    cuda_numneigh_inner = (int*)allocateGPU(atom->Nclusters_max * sizeof(int));
    cuda_neighbors      = (int*)allocateGPU(
        atom->Nclusters_max * neighbor->maxneighs * sizeof(int));
    cuda_neighbors_imask = (unsigned int*)allocateGPU(
        atom->Nclusters_max * neighbor->maxneighs * sizeof(unsigned int));
    // Pinned staging buffers, refilled on every list upload
    natoms = (int*)allocateHostPinned(atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
    ngatoms = (int*)allocateHostPinned(
        atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
}

extern "C" void copyDataToCUDADevice(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    memcpyToGPU(cuda_cl_x,
        atom->cl_x,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost + 1) * CLUSTER_M * ATOM_DIM * sizeof(MD_FLOAT));
    memcpyToGPU(cuda_cl_v,
        atom->cl_v,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));

    if (param->lj_comb_rule == LJ_COMB_FULL) {
        memcpyToGPU(cuda_cl_t,
            atom->cl_t,
            (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
                sizeof(int));

        // cutforcesq, sigma6 and epsilon per type pair
        const int ntypes2 = atom->ntypes * atom->ntypes;
        memcpyToGPU(cuda_cutforcesq, atom->cutforcesq, ntypes2 * sizeof(MD_FLOAT));
        memcpyToGPU(cuda_sigma6, atom->sigma6, ntypes2 * sizeof(MD_FLOAT));
        memcpyToGPU(cuda_epsilon, atom->epsilon, ntypes2 * sizeof(MD_FLOAT));
        MD_FLOAT4* nbfp = (MD_FLOAT4*)malloc(ntypes2 * sizeof(MD_FLOAT4));
        for (int t = 0; t < ntypes2; t++) {
            nbfp[t].x = atom->cutforcesq[t];
            nbfp[t].y = atom->sigma6[t];
            nbfp[t].z = atom->epsilon[t];
            nbfp[t].w = (MD_FLOAT)0.0;
        }
        memcpyToGPU(cuda_nbfp, nbfp, ntypes2 * sizeof(MD_FLOAT4));
        free(nbfp);
    }
    if (param->lj_comb_rule == LJ_COMB_GEOM) {
        memcpyToGPU(cuda_cl_sqrt_epsilon,
            atom->cl_sqrt_epsilon,
            (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
                sizeof(MD_FLOAT));
        memcpyToGPU(cuda_cl_sigma3,
            atom->cl_sigma3,
            (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
                sizeof(MD_FLOAT));
    }
#if COULOMB_RULE != COULOMB_NONE
    memcpyToGPU(cuda_cl_q,
        atom->cl_q,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
            sizeof(MD_FLOAT));
#endif
    if (atom->molsize > 1) {
        memcpyToGPU(cuda_cl_mol,
            atom->cl_mol,
            (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
                sizeof(int));
    }
    memcpyToGPU(cuda_cl_id,
        atom->cl_id,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M *
            sizeof(int));
    launchPackAtomsCuda(param, atom, 1);

    for (int ci = 0; ci < atom->Nclusters_local * SCLUSTER_SIZE; ci++) {
        natoms[ci] = atom->iclusters[ci].natoms;
    }

    // int ncj  = get_ncj_from_nci(atom->Nclusters_local);
    int ncj = CJ0_FROM_CI(atom->Nclusters_local);
    for (int cg = 0; cg < atom->Nclusters_ghost; cg++) {
        const int cj = ncj + cg;
        ngatoms[cg]  = atom->jclusters[cj].natoms;
    }

    memcpyToGPU(cuda_natoms, natoms, atom->Nclusters_local * SCLUSTER_SIZE * sizeof(int));
    memcpyToGPU(cuda_jclusters_natoms, ngatoms, atom->Nclusters_ghost * sizeof(int));
#ifndef _MPI
    memcpyToGPU(cuda_border_map, atom->border_map, atom->Nclusters_ghost * sizeof(int));
    memcpyToGPU(cuda_PBCx, atom->PBCx, atom->Nclusters_ghost * sizeof(int));
    memcpyToGPU(cuda_PBCy, atom->PBCy, atom->Nclusters_ghost * sizeof(int));
    memcpyToGPU(cuda_PBCz, atom->PBCz, atom->Nclusters_ghost * sizeof(int));
#endif

    memcpyToGPU(cuda_numneigh, neighbor->numneigh, atom->Nclusters_local * sizeof(int));
    memcpyToGPU(cuda_numneigh_inner,
        neighbor->numneigh_inner,
        atom->Nclusters_local * sizeof(int));
    memcpyToGPU(cuda_neighbors,
        neighbor->neighbors,
        atom->Nclusters_local * neighbor->maxneighs * sizeof(int));

    if (param->super_clustering) {
        memcpyToGPU(cuda_neighbors_imask,
            neighbor->neighbors_imask,
            atom->Nclusters_local * neighbor->maxneighs * sizeof(unsigned int));
        launchPackNeighborsSup(atom, neighbor);
    } else {
        memsetGPU(cuda_neighbors_imask,
            0xff,
            atom->Nclusters_local * neighbor->maxneighs * sizeof(unsigned int));
    }
}

extern "C" void copyDataFromCUDADevice(Parameter* param, Atom* atom)
{
    memcpyFromGPU(atom->cl_x,
        cuda_cl_x,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost + 1) * CLUSTER_M * ATOM_DIM * sizeof(MD_FLOAT));
    memcpyFromGPU(atom->cl_v,
        cuda_cl_v,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 *
            sizeof(MD_FLOAT));
}

extern "C" void cudaDeviceFree(Parameter* param)
{
    GPUfree(cuda_cl_x);
    GPUfree(cuda_cl_v);
    GPUfree(cuda_cl_f);
    GPUfree(cuda_cl_t);
    GPUfree(cuda_cl_sqrt_epsilon);
    GPUfree(cuda_cl_sigma3);
    GPUfree(cuda_cutforcesq);
    GPUfree(cuda_sigma6);
    GPUfree(cuda_epsilon);
    GPUfree(cuda_nbfp);
#if COULOMB_RULE != COULOMB_NONE
    GPUfree(cuda_cl_q);
#endif
    GPUfree(cuda_cl_mol);
    GPUfree(cuda_cl_id);
    GPUfree(cuda_cl_xq);
    GPUfree(cuda_cl_par);
    GPUfree(cuda_numneigh);
    GPUfree(cuda_numneigh_inner);
    GPUfree(cuda_neighbors);
    GPUfree(cuda_neighbors_imask);
    GPUfree(cuda_natoms);
    GPUfree(cuda_border_map);
    GPUfree(cuda_jclusters_natoms);
    GPUfree(cuda_PBCx);
    GPUfree(cuda_PBCy);
    GPUfree(cuda_PBCz);

    freeHostPinned(natoms);
    freeHostPinned(ngatoms);
    freeSupData();
}

extern "C" void copyGhostFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->cl_x,
        cuda_cl_x,
        atom->Nclusters_local * CLUSTER_M * SCLUSTER_SIZE * ATOM_DIM * sizeof(MD_FLOAT));
}

extern "C" void copyGhostToGPU(Atom* atom)
{
    int ncj               = atom->Nclusters_local * SCLUSTER_SIZE;
    int cj_vec_base       = CJ_VECTOR_BASE_INDEX(ncj);
    MD_FLOAT* host_cl_x   = &atom->cl_x[cj_vec_base];
    MD_FLOAT* device_cl_x = &cuda_cl_x[cj_vec_base];
    memcpyToGPU(device_cl_x,
        host_cl_x,
        atom->Nclusters_ghost * CLUSTER_N * ATOM_DIM * sizeof(MD_FLOAT));
}

extern "C" void copyForceFromGPU(Atom* atom)
{
    memcpyFromGPU(atom->cl_f,
        cuda_cl_f,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 * sizeof(MD_FLOAT));
}

extern "C" void copyForceToGPU(Atom* atom)
{
    memcpyToGPU(cuda_cl_f,
        atom->cl_f,
        (atom->Nclusters_local * SCLUSTER_SIZE + atom->Nclusters_ghost) * CLUSTER_M * 3 * sizeof(MD_FLOAT));
}

extern "C" void growClustersCUDA(Atom* atom)
{
    if (cuda_cl_x) {
        cuda_cl_x             = (MD_FLOAT*)reallocateGPU(cuda_cl_x,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * ATOM_DIM *
                sizeof(MD_FLOAT));
        cuda_cl_v             = (MD_FLOAT*)reallocateGPU(cuda_cl_v,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
        cuda_cl_f             = (MD_FLOAT*)reallocateGPU(cuda_cl_f,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * 3 * sizeof(MD_FLOAT));
        cuda_natoms           = (int*)reallocateGPU(cuda_natoms,
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
        cuda_jclusters_natoms = (int*)reallocateGPU(cuda_jclusters_natoms,
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
        cuda_numneigh         = (int*)reallocateGPU(cuda_numneigh,
            atom->Nclusters_max * sizeof(int));
        cuda_numneigh_inner   = (int*)reallocateGPU(cuda_numneigh_inner,
            atom->Nclusters_max * sizeof(int));

        freeHostPinned(natoms);
        freeHostPinned(ngatoms);

        natoms = (int*)allocateHostPinned(
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
        ngatoms = (int*)allocateHostPinned(
            atom->Nclusters_max * SCLUSTER_SIZE * sizeof(int));
        cuda_cl_t = (int*)reallocateGPU(cuda_cl_t,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
        cuda_cl_sqrt_epsilon = (MD_FLOAT*)reallocateGPU(cuda_cl_sqrt_epsilon,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
        cuda_cl_sigma3       = (MD_FLOAT*)reallocateGPU(cuda_cl_sigma3,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
#if COULOMB_RULE != COULOMB_NONE
        cuda_cl_q = (MD_FLOAT*)reallocateGPU(cuda_cl_q,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT));
#endif
        cuda_cl_mol = (int*)reallocateGPU(cuda_cl_mol,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
        cuda_cl_id  = (int*)reallocateGPU(cuda_cl_id,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(int));
        cuda_cl_xq  = (MD_FLOAT4*)reallocateGPU(cuda_cl_xq,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(MD_FLOAT4));
        cuda_cl_par = (AtomParam*)reallocateGPU(cuda_cl_par,
            atom->Nclusters_max * CLUSTER_M * SCLUSTER_SIZE * sizeof(AtomParam));
    }
}

extern void growNeighborCUDA(Atom* atom, Neighbor* neighbor)
{
    cuda_neighbors       = (int*)reallocateGPU(cuda_neighbors,
        atom->Nclusters_max * neighbor->maxneighs * sizeof(int));
    cuda_neighbors_imask = (unsigned int*)reallocateGPU(cuda_neighbors_imask,
        atom->Nclusters_max * neighbor->maxneighs * sizeof(unsigned int));
}

/* update coordinates of ghost atoms */
/* uses mapping created in setupPbc */
extern "C" void updatePbcCUDA(Atom* atom, Parameter* param)
{
    if (param->super_clustering) {
        updatePbcCUDASup(atom, param);
    } else {
        updatePbcCUDASimple(atom, param);
    }
}

extern "C" void reverseGhostForcesCUDA(Atom* atom, Parameter* param)
{
    if (param->super_clustering) {
        reverseGhostForcesCUDASup(atom, param);
    } else {
        reverseGhostForcesCUDASimple(atom, param);
    }
}

extern "C" void pruneNeighborCUDA(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    if (param->super_clustering) {
        pruneNeighborCUDASup(param, atom, neighbor);
    } else {
        pruneNeighborCUDASimple(param, atom, neighbor);
    }
}

// Resident blocks per SM of a force kernel (one wave = N_SM of these)
static int residentBlocksPerSM(const void* kernel)
{
    static const void* kernels[8];
    static int blocks[8];
    static int nkernels = 0;
    for (int i = 0; i < nkernels; i++) {
        if (kernels[i] == kernel) return blocks[i];
    }

    int n = 0;
    cuda_assert("cudaOccupancyMaxActiveBlocksPerMultiprocessor",
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&n, kernel, CLUSTER_N * CLUSTER_M, 0));
    if (n < 1) n = 1;
    if (nkernels < 8) {
        kernels[nkernels] = kernel;
        blocks[nkernels]  = n;
        nkernels++;
    }
    return n;
}

// j-slices S = K * sqrt(N_SM * n_tot) / N_c, at most 10 waves of resident blocks
extern "C" int chooseJSlicesK(
    const Atom* atom, const Neighbor* neighbor, double k_slice, const void* kernel)
{
    static int nsm = 0;
    if (nsm == 0) {
        int dev = 0;
        cuda_assert("cudaGetDevice", cudaGetDevice(&dev));
        cuda_assert("cudaDeviceGetAttribute",
            cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev));
    }
    const int blocks_per_sm = residentBlocksPerSM(kernel);

    const int nc = atom->Nclusters_local;
    if (nc <= 0) return 1;

    long n_tot = 0;
    int  n_max = 1;
    for (int ci = 0; ci < nc; ci++) {
        const int n = neighbor->numneigh[ci];
        n_tot += n;
        if (n > n_max) n_max = n;
    }

    int S = (int)lround(k_slice * sqrt((double)nsm) * sqrt((double)n_tot) / (double)nc);
    const long max_blocks = 10L * nsm * blocks_per_sm; // 10 waves
    if (S > max_blocks / nc) S = (int)(max_blocks / nc);
    if (S < 1)     S = 1;
    if (S > n_max) S = n_max;
    if (S > 65535) S = 65535;
    return S;
}

/*
extern "C" void forwardCommCUDA(Comm* comm, Atom* atom, int iswap)
{
    int nrqst = 0, offset = 0, nsend = 0, nrecv = 0;
    int pbc[3];
    int size    = comm->forwardSize;
    int maxrqst = comm->numneigh;
    int* cuda_sendlist;
    int max_list_size;
    cuda_buf_send =  (MD_FLOAT*)allocateGPU(comm->max_send * sizeof(MD_FLOAT));
    cuda_buf_recv =  (MD_FLOAT*)allocateGPU(comm->max_recv * sizeof(MD_FLOAT));

    //use a single buffer and takes the highes list size to move list of cluster to send
    for (int ineigh = 0; ineigh < comm->numneigh; ineigh){
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

        const int threads_num = 64;
        dim3 block_size       = dim3(threads_num, 1, 1);
        dim3 grid_size = dim3(( comm->atom_send[ineigh] + threads_num - 1) / threads_num,
1, 1); packForwardCuda<<<grid_size, block_size>>>( cuda_cl_x, //MD_FLOAT*  -->need to be
in tye device cuda_jclusters_natoms,       //int*       -->need to be in tye device
                                            comm->atom_send[ineigh],     //int
                                            cuda_sendlist,               //int* -->need to
be in tye device cuda_buf_send[offset*size],  //MD_FLOAT*  -->need to be in tye device
                                            pbc[_x],                     //int
                                            pbc[_y],                     //int
                                            pbc[_z],                     //int
                                            atom->xprd,                  //MD_FLOAT
                                            atom->yrpd,                  //MD_FLOAT
                                            atom->zprd);                 //MD_FLOAT
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

    //unpack buffer
    for (int ineigh = comm->recvfrom[iswap]; ineigh < comm->recvtill[iswap]; ineigh++) {
        offset = comm->off_atom_recv[ineigh];
        MD_FLOAT *buf = (comm->othersend[iswap]) ? cuda_buf_recv : cuda_buf_send;
        const int threads_num = 64;
        dim3 block_size       = dim3(threads_num, 1, 1);
        dim3 grid_size = dim3(( comm->atom_recv[ineigh] + threads_num - 1) / threads_num,
1, 1);

        unpackForwardCuda<<<grid_size, block_size>>>(
                                        cuda_cl_x,                          //MD_FLOAT*
--> need to be in the device cuda_jclusters_natoms,              //int*      --> need to
be in the device comm->atom_recv[ineigh],            //int comm->firstrecv[iswap] +
offset,    //int &buf[offset * size]);               //MD_FLOAT* --> need to be in the
devic
    }
    GPUfree(cuda_sendlist);
    GPUfree(cuda_buf_recv);
    GPUfree(cuda_buf_send);
}

__global__ void packForwardCuda(MD_FLOAT* cuda_cl_x, int* cuda_jclusters_natoms, int nc,
int* cuda_list, MD_FLOAT* cuda_buf, int PBCx, int PBCY, int PBCz, MD_FLOAT xprd, MD_FLOAT
yrpd, MD_FLOAT zprd)
{
    unsigned int j = blockDim.x * blockIdx.x + threadIdx.x;
    if (j >= nc) return;

    int cj          = cuda_list[j];
    int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
    MD_FLOAT* cj_x  = &cuda_cl_x[cj_vec_base];
    int displ       = j * CLUSTER_N;

    for (int cjj = 0; cjj < cuda_jclusters_natoms[cj]; cjj++) {
        cuda_buf[3 * (displ + cjj) + 0] = cj_x[CL_X_INDEX(cjj)] +
                                    PBCx * xprd;
        cuda_buf[3 * (displ + cjj) + 1] = cj_x[CL_Y_INDEX(cjj)] +
                                    PBCy * yprd;
        cuda_buf[3 * (displ + cjj) + 2] = cj_x[CL_Z_INDEX(cjj)] +
                                    PBCz * zprd;
    }

    for (int cjj = cuda_jclusters_natoms[cj]; cjj < CLUSTER_N; cjj++) {
        cuda_buf[3 * (displ + cjj) + 0] = -1; // x
        cuda_buf[3 * (displ + cjj) + 1] = -1; // y
        cuda_buf[3 * (displ + cjj) + 2] = -1; // z
    }

}
__global__ void unpackForwardCuda(MD_FLOAT* cuda_cl_x, int* cuda_jclusters_natoms, int nc,
int c0, MD_FLOAT* cuda_buf)
{
    unsigned int j = blockDim.x * blockIdx.x + threadIdx.x;
    if (j >= nc) return;
    int cj = c0 + j;
    int cj_vec_base = CJ_VECTOR_BASE_INDEX(cj);
    MD_FLOAT* cj_x  = &cuda_cl_x[cj_vec_base];
    int displ       = j * CLUSTER_N;

    for (int cjj = 0; cjj < cuda_jclusters_natoms[cj]; cjj++) {
        if (cj_x[CL_X_INDEX(cjj)] < INF)
            cj_x[CL_X_INDEX(cjj)] = cuda_buf[3 * (displ + cjj) + 0];
        if (cj_x[CL_Y_INDEX(cjj)] < INF)
            cj_x[CL_Y_INDEX(cjj)] = cuda_buf[3 * (displ + cjj) + 1];
        if (cj_x[CL_Z_INDEX(cjj)] < INF)
            cj_x[CL_Z_INDEX(cjj)] = cuda_buf[3 * (displ + cjj) + 2];
    }
}
*/
