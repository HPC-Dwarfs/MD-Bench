#include <atom.h>
#include <box.h>
#include <force.h>
#include <grid.h>
#include <parameter.h>

#ifndef COMM_H
#define COMM_H

#ifdef CLUSTER_PAIR
#define FORWARD_SIZE  (3 * CLUSTER_N)
#define REVERSE_SIZE  (3 * CLUSTER_N)
#define GHOST_SIZE    (4 * CLUSTER_N + 10)
#define EXCHANGE_SIZE 7
#define JFAC          MAX(1, CLUSTER_N / CLUSTER_M)
#ifdef CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS
#define LOCAL         (atom->Nclusters_local * SCLUSTER_SIZE)
#else
#define LOCAL         (atom->Nclusters_local / JFAC)
#endif
#define GHOST         atom->Nclusters_ghost

#define IsInRegionToSend(cj)                                                             \
       ((atom->jclusters[(cj)].bbmaxx >= xlo && atom->jclusters[(cj)].bbminx < xhi) &&   \
        (atom->jclusters[(cj)].bbmaxy >= ylo && atom->jclusters[(cj)].bbminy < yhi) &&   \
        (atom->jclusters[(cj)].bbmaxz >= zlo && atom->jclusters[(cj)].bbminz < zhi))
#else
#define FORWARD_SIZE  3
#define REVERSE_SIZE  3
#define GHOST_SIZE    4
#define EXCHANGE_SIZE 7
#define LOCAL         atom->Nlocal
#define GHOST         atom->Nghost

#define IsInRegionToSend(i)                                                              \
       ((atom_x((i)) >= xlo && atom_x((i)) < xhi) &&                                     \
        (atom_y((i)) >= ylo && atom_y((i)) < yhi) &&                                     \
        (atom_z((i)) >= zlo && atom_z((i)) < zhi))
#endif

typedef struct {
    int myproc;  // my proc ID
    int numproc; // # of processors
    int numneigh;    // # of all my neighs along all swaps
    int maxneigh;    // Buffer size for my neighs
    int sendfrom[6]; // return the lowest neigh index to send in each swap
    int sendtill[6]; // return the highest neigh index to send in each swao
    int recvfrom[6]; // return the lowest neigh index to recv in each swap
    int recvtill[6]; // return the highest neigh index to recv in each swap
    int* nsend;      // neigh whose I want to send
    int* nrecv;      // neigh whose I want to recv

    int* pbc_x; // if pbc in x
    int* pbc_y; // if pbc in y
    int* pbc_z; // if pbc in z

    int *atom_send, *atom_recv; // # of atoms to send/recv for each of my neighs
    int* off_atom_send;         // atom offset to send, inside of a swap
    int* off_atom_recv;         // atom offset to recv, inside of a swap

    int* nexch;       // procs to exchange
    int numneighexch; // # of neighbours to exchange
    int maxneighexch; // max buff size to store neighbours

    int numswap;      // # of swaps to perform, it is 6
    int swapdim[6];   // dimension of the swap (_x, _y or _z)
    int swapdir[6];   // direction of the swap 0 or 1
    int swap[3][2];   // given a dim and dir, knows the swap
    int othersend[6]; // Determine if a proc interact with more procs in a given
                      // swap

    int firstrecv[6]; // where to put 1st recv atom in each swap
    int** sendlist;   // list of atoms to send in each swap
    int* maxsendlist; // max # of atoms send in each list-swap

    int maxsend;        // max elements in buff sender
    int maxrecv;        // max elements in buff receiver
    MD_FLOAT* buf_send; // sender buffer for all comm
    MD_FLOAT* buf_recv; // receicer buffer for all comm

    int forwardSize;  // # of paramaters per atom in forward comm.
    int reverseSize;  // # of parameters per atom in reverse
    int exchangeSize; // # of parameters per atom in exchange
    int ghostSize;    // # of parameters per atom in ghost list

    int iterAtom; // last atom to iterate in each swap.
    Box* boxes;   // Boundaries to  be sent to other procs as ghost.
} Comm;

void initComm(int*, char***, Comm*); // Init MPI
void endComm(Comm*);                 // End MPI
void barrierComm();                  // Barrier comm

#ifdef _MPI
void setupComm(Comm*, Parameter*, Grid*); // Creates a 3d grid or rcb grid
void neighComm(Comm*,
    Parameter*,
    Grid*); // Find neighbours within cut-off and defines ghost regions
void forwardComm(Comm*, Atom*, int); // Send info in one direction
void reverseComm(Comm*, Atom*, int); // Return info after forward
                                     // communication
void exchangeComm(Comm*, Atom*);     // Exchange info between procs
void ghostComm(Comm*, Parameter*, Atom*, int);   // Build the ghost neighbours to send during next forwards
void growSend(Comm*, int);      // Grows the size of the buffer sender
void growRecv(Comm*, int);      // Grows the size of the buffer receiver
void growList(Comm*, int, int); // Grows the size of the list to send
#endif
#endif
