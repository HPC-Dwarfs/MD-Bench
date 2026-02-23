/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <atom.h>
#include <neighbor.h>
#include <parameter.h>
#include <util.h>

#define SMALL  1.0e-6
#define FACTOR 0.999

#ifdef CUDA_TARGET
BuildNeighborFunction buildNeighbor = buildNeighborCUDA;
#else
BuildNeighborFunction buildNeighbor = buildNeighborCPU;
#endif

MD_FLOAT xprd, yprd, zprd;
MD_FLOAT bininvx, bininvy, bininvz;
int mbinxlo, mbinylo, mbinzlo; // not used
int nbinx, nbiny, nbinz;
int mbinx, mbiny, mbinz; // n bins in x, y, z
int* bincount;
int* bins;
int mbins;         // total number of bins
int atoms_per_bin; // max atoms per bin
MD_FLOAT cutneigh;
MD_FLOAT cutneighsq; // neighbor cutoff squared
int nmax;
int nstencil; // # of bins in stencil
int* stencil; // stencil list of bin offsets
MD_FLOAT binsizex, binsizey, binsizez;
static int coord2bin(MD_FLOAT, MD_FLOAT, MD_FLOAT);
static MD_FLOAT bindist(int, int, int);

// MPI version
int pad_x, pad_y, pad_z; // version MPI
int me;                  // rank
int method;              // method
int half_stencil;        // If half stencil exist
int shellMethod;         // If shell method exist
static int ghostZone(Atom*, int);
static int eightZone(Atom*, int);
static int halfZone(Atom*, int);
static void neighborGhost(Atom*, Neighbor*);
static inline int skipNeigh(Atom* atom, int i, int j);

/* exported subroutines */
void initNeighbor(Neighbor* neighbor, Parameter* param)
{
    MD_FLOAT neighscale = 5.0 / 6.0;
    xprd                = param->nx * param->lattice;
    yprd                = param->ny * param->lattice;
    zprd                = param->nz * param->lattice;
    cutneigh            = param->cutneigh;
    nbinx               = MAX(1, neighscale * param->nx);
    nbiny               = MAX(1, neighscale * param->ny);
    nbinz               = MAX(1, neighscale * param->nz);
    nmax                = 0;
    atoms_per_bin       = 8;
    stencil             = NULL;
    bins                = NULL;
    bincount            = NULL;
    neighbor->maxneighs = 100;
    neighbor->numneigh  = NULL;
    neighbor->neighbors = NULL;
    //========== MPI =============
    method = param->method;
    if (method == halfShell || method == eightShell) {
        param->half_neigh = 1;
        shellMethod       = 1;
    }

    if (method == halfStencil) {
        param->half_neigh = 0;
        half_stencil      = 1;
    }
    neighbor->half_neigh = param->half_neigh;
    
    me = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#endif
    neighbor->Nshell        = 0;
    neighbor->numNeighShell = NULL;
    neighbor->neighshell    = NULL;
    neighbor->listshell     = NULL;
}

void setupNeighbor(Parameter* param)
{
    MD_FLOAT coord;
    int mbinxhi, mbinyhi, mbinzhi;
    int nextx, nexty, nextz;

    if (param->input_file != NULL) {
        xprd = param->xprd;
        yprd = param->yprd;
        zprd = param->zprd;
    }

    // TODO: update lo and hi for standard case and use them here instead
    MD_FLOAT xlo = 0.0;
    MD_FLOAT xhi = xprd;
    MD_FLOAT ylo = 0.0;
    MD_FLOAT yhi = yprd;
    MD_FLOAT zlo = 0.0;
    MD_FLOAT zhi = zprd;

    cutneighsq = cutneigh * cutneigh;

    if (param->input_file != NULL) {
        binsizex = cutneigh * 0.5;
        binsizey = cutneigh * 0.5;
        binsizez = cutneigh * 0.5;
        nbinx    = (int)((param->xhi - param->xlo) / binsizex);
        nbiny    = (int)((param->yhi - param->ylo) / binsizey);
        nbinz    = (int)((param->zhi - param->zlo) / binsizez);
        if (nbinx == 0) {
            nbinx = 1;
        }
        if (nbiny == 0) {
            nbiny = 1;
        }
        if (nbinz == 0) {
            nbinz = 1;
        }
        bininvx = nbinx / (param->xhi - param->xlo);
        bininvy = nbiny / (param->yhi - param->ylo);
        bininvz = nbinz / (param->zhi - param->zlo);
    } else {
        binsizex = xprd / nbinx;
        binsizey = yprd / nbiny;
        binsizez = zprd / nbinz;
        bininvx  = 1.0 / binsizex;
        bininvy  = 1.0 / binsizey;
        bininvz  = 1.0 / binsizez;
    }
   
    pad_x = (int)(cutneigh * bininvx);
    while (pad_x * binsizex < FACTOR * cutneigh)
        pad_x++;
    pad_y = (int)(cutneigh * bininvy);
    while (pad_y * binsizey < FACTOR * cutneigh)
        pad_y++;
    pad_z = (int)(cutneigh * bininvz);
    while (pad_z * binsizez < FACTOR * cutneigh)
        pad_z++;

    coord   = xlo - cutneigh - SMALL * xprd;
    mbinxlo = (int)(coord * bininvx);
    if (coord < 0.0) {
        mbinxlo = mbinxlo - 1;
    }
    coord   = xhi + cutneigh + SMALL * xprd;
    mbinxhi = (int)(coord * bininvx);

    coord   = ylo - cutneigh - SMALL * yprd;
    mbinylo = (int)(coord * bininvy);
    if (coord < 0.0) {
        mbinylo = mbinylo - 1;
    }
    coord   = yhi + cutneigh + SMALL * yprd;
    mbinyhi = (int)(coord * bininvy);

    coord   = zlo - cutneigh - SMALL * zprd;
    mbinzlo = (int)(coord * bininvz);
    if (coord < 0.0) {
        mbinzlo = mbinzlo - 1;
    }
    coord   = zhi + cutneigh + SMALL * zprd;
    mbinzhi = (int)(coord * bininvz);

    mbinxlo = mbinxlo - 1;
    mbinxhi = mbinxhi + 1;
    mbinylo = mbinylo - 1;
    mbinyhi = mbinyhi + 1;

    mbinzlo = mbinzlo - 1;
    mbinzhi = mbinzhi + 1;

/*
    mbinxlo = mbinxlo - 1;
    mbinxhi = mbinxhi + 1;
    mbinx   = mbinxhi - mbinxlo + 1;

    mbinylo = mbinylo - 1;
    mbinyhi = mbinyhi + 1;
    mbiny   = mbinyhi - mbinylo + 1;

    mbinzlo = mbinzlo - 1;
    mbinzhi = mbinzhi + 1;
    mbinz   = mbinzhi - mbinzlo + 1;
*/
    nextx = (int)(cutneigh * bininvx);
    while (nextx * binsizex < FACTOR * cutneigh) { 
        nextx++;
        pad_x++;
    }

    nexty = (int)(cutneigh * bininvy);
    while (nexty * binsizey < FACTOR * cutneigh) {
        nexty++;
        pad_y++;
    }

    nextz = (int)(cutneigh * bininvz);
    while (nextz * binsizez < FACTOR * cutneigh) {
        nextz++;
        pad_z++;
    }

    mbinx = MAX(1, nbinx + 4 * pad_x);
    mbiny = MAX(1, nbiny + 4 * pad_y);
    mbinz = MAX(1, nbinz + 4 * pad_z);

    if (stencil) {
        free(stencil);
    }
    stencil = (int*)malloc(
        (2 * nextz + 1) * (2 * nexty + 1) * (2 * nextx + 1) * sizeof(int));
    nstencil   = 0;
    int kstart = -nextz;
    int ibin   = 0;
    for (int k = kstart; k <= nextz; k++) {
        for (int j = -nexty; j <= nexty; j++) {
            for (int i = -nextx; i <= nextx; i++) {
                if (bindist(i, j, k) < cutneighsq) {
                    int jbin = k * mbiny * mbinx + j * mbinx + i;
                    if (ibin > jbin && half_stencil) continue;
                    stencil[nstencil++] = jbin;
                    // stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
                }
            }
        }
    }

    mbins = mbinx * mbiny * mbinz;
    if (bincount) {
        free(bincount);
    }
    bincount = (int*)malloc(mbins * sizeof(int));
    if (bins) {
        free(bins);
    }
    bins = (int*)malloc(mbins * atoms_per_bin * sizeof(int));
}

void buildNeighborCPU(Atom* atom, Neighbor* neighbor) {
    DEBUG_MESSAGE("buildNeighborCPU begin\n");
    int nall = atom->Nlocal + atom->Nghost;

    /* extend atom arrays if necessary */
    if (nall > nmax) {
        nmax = nall;
        if (neighbor->numneigh) free(neighbor->numneigh);
        if (neighbor->neighbors) free(neighbor->neighbors);
        neighbor->numneigh  = (int*)malloc(nmax * sizeof(int));
        neighbor->neighbors = (int*)malloc(nmax * neighbor->maxneighs * sizeof(int*));
    }

    /* bin local & ghost atoms */
    binatoms(atom);
    int resize = 1;

    /* loop over each atom, storing neighbors */
    while (resize) {
        int new_maxneighs = neighbor->maxneighs;
        resize            = 0;

        for (int i = 0; i < atom->Nlocal; i++) {
            int n         = 0;
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            int ibin      = coord2bin(xtmp, ytmp, ztmp);
#ifndef ONE_ATOM_TYPE
            int type_i = atom->type[i];
#endif
            for (int k = 0; k < nstencil; k++) {
                int jbin     = ibin + stencil[k];
                int* loc_bin = &bins[jbin * atoms_per_bin];

                for (int m = 0; m < bincount[jbin]; m++) {
                    int j = loc_bin[m];
                    if (i == j) continue;
                    if (neighbor->half_neigh && j < i) continue;
                    if (half_stencil && ibin == jbin && skipNeigh(atom, i, j)) continue;
                    // if ((j == i) || (neighbor->half_neigh && (j < i))) continue;

                    MD_FLOAT delx = xtmp - atom_x(j);
                    MD_FLOAT dely = ytmp - atom_y(j);
                    MD_FLOAT delz = ztmp - atom_z(j);
                    MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#ifndef ONE_ATOM_TYPE
                    int type_j = atom->type[j];
                    const MD_FLOAT cutoff =
                        atom->cutneighsq[type_i * atom->ntypes + type_j];
#else
                    const MD_FLOAT cutoff = cutneighsq;
#endif
                    if (rsq <= cutoff) {
                        neighs(neighbor->neighbors, i, n, atom->Nlocal, neighbor->maxneighs) = j;
                        n++;
                    }
                }
            }

            neighbor->numneigh[i] = n;
            if (n >= neighbor->maxneighs) {
                resize = 1;

                if (n >= new_maxneighs) {
                    new_maxneighs = n;
                }
            }
        }

        if (resize) {
            printf("RESIZE %d, PROC %d\n", neighbor->maxneighs, me);
            neighbor->maxneighs = new_maxneighs * 1.2;
            free(neighbor->neighbors);
            neighbor->neighbors = (int*)malloc(
                atom->Nmax * neighbor->maxneighs * sizeof(int));
        }
    }

    if (method == eightShell) {
        neighborGhost(atom, neighbor);
    }

    DEBUG_MESSAGE("buildNeighborCPU end\n");
}

/* internal subroutines */
MD_FLOAT bindist(int i, int j, int k) {
    MD_FLOAT delx, dely, delz;

    if (i > 0) {
        delx = (i - 1) * binsizex;
    } else if (i == 0) {
        delx = 0.0;
    } else {
        delx = (i + 1) * binsizex;
    }

    if (j > 0) {
        dely = (j - 1) * binsizey;
    } else if (j == 0) {
        dely = 0.0;
    } else {
        dely = (j + 1) * binsizey;
    }

    if (k > 0) {
        delz = (k - 1) * binsizez;
    } else if (k == 0) {
        delz = 0.0;
    } else {
        delz = (k + 1) * binsizez;
    }

    return (delx * delx + dely * dely + delz * delz);
}

int coord2bin(MD_FLOAT xin, MD_FLOAT yin, MD_FLOAT zin) {
    int ix, iy, iz;
    MD_FLOAT eps = 1e-9;
    MD_FLOAT xlo = 0.0;
    MD_FLOAT ylo = 0.0;
    MD_FLOAT zlo = 0.0;
    xlo          = fabs(xlo - pad_x * binsizex) + eps;
    ylo          = fabs(ylo - pad_y * binsizey) + eps;
    zlo          = fabs(zlo - pad_z * binsizez) + eps;
    ix           = (int)((xin + xlo) * bininvx);
    iy           = (int)((yin + ylo) * bininvy);
    iz           = (int)((zin + zlo) * bininvz);

    return (iz * mbiny * mbinx + iy * mbinx + ix);

    /*
    if (xin >= xprd) {
        ix = (int)((xin - xprd) * bininvx) + nbinx - mbinxlo;
    } else if (xin >= 0.0) {
        ix = (int)(xin * bininvx) - mbinxlo;
    } else {
        ix = (int)(xin * bininvx) - mbinxlo - 1;
    }

    if (yin >= yprd) {
        iy = (int)((yin - yprd) * bininvy) + nbiny - mbinylo;
    } else if (yin >= 0.0) {
        iy = (int)(yin * bininvy) - mbinylo;
    } else {
        iy = (int)(yin * bininvy) - mbinylo - 1;
    }

    if (zin >= zprd) {
        iz = (int)((zin - zprd) * bininvz) + nbinz - mbinzlo;
    } else if (zin >= 0.0) {
        iz = (int)(zin * bininvz) - mbinzlo;
    } else {
        iz = (int)(zin * bininvz) - mbinzlo - 1;
    }

    return (iz * mbiny * mbinx + iy * mbinx + ix + 1);
    */
}

void binatoms(Atom* atom) {
    DEBUG_MESSAGE("binatoms begin\n");

    int nall   = atom->Nlocal + atom->Nghost;
    int resize = 1;

    while (resize > 0) {
        resize = 0;

        for (int i = 0; i < mbins; i++) {
            bincount[i] = 0;
        }

        for (int i = 0; i < nall; i++) {
            int ibin = coord2bin(atom_x(i), atom_y(i), atom_z(i));
            if (shellMethod && !ghostZone(atom, i)) continue;
            if (bincount[ibin] < atoms_per_bin) {
                int ac                          = bincount[ibin]++;
                bins[ibin * atoms_per_bin + ac] = i;
            } else {
                resize = 1;
            }
        }

        if (resize) {
            free(bins);
            atoms_per_bin *= 2;
            bins = (int*)malloc(mbins * atoms_per_bin * sizeof(int));
        }
    }

    DEBUG_MESSAGE("binatoms end\n");
}

void sortAtom(Atom* atom) {
    DEBUG_MESSAGE("sortAtom begin\n");

    binatoms(atom);
    int Nmax    = atom->Nmax;
    int* binpos = bincount;

    for (int i = 1; i < mbins; i++) {
        binpos[i] += binpos[i - 1];
    }

#ifdef ATOM_POSITION_AOS
    MD_FLOAT* new_x  = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT) * 3);
    MD_FLOAT* new_vx = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT) * 3);
#else
    MD_FLOAT* new_x  = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_y  = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_z  = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vx = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vy = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT));
    MD_FLOAT* new_vz = (MD_FLOAT*)malloc(Nmax * sizeof(MD_FLOAT));
#endif
    MD_FLOAT* old_x  = atom->x;
    MD_FLOAT* old_y  = atom->y;
    MD_FLOAT* old_z  = atom->z;
    MD_FLOAT* old_vx = atom->vx;
    MD_FLOAT* old_vy = atom->vy;
    MD_FLOAT* old_vz = atom->vz;

    for (int mybin = 0; mybin < mbins; mybin++) {
        int start = mybin > 0 ? binpos[mybin - 1] : 0;
        int count = binpos[mybin] - start;
        for (int k = 0; k < count; k++) {
            int new_i = start + k;
            int old_i = bins[mybin * atoms_per_bin + k];
#ifdef ATOM_POSITION_AOS
            new_x[new_i * 3 + 0]  = old_x[old_i * 3 + 0];
            new_x[new_i * 3 + 1]  = old_x[old_i * 3 + 1];
            new_x[new_i * 3 + 2]  = old_x[old_i * 3 + 2];
            new_vx[new_i * 3 + 0] = old_vx[old_i * 3 + 0];
            new_vx[new_i * 3 + 1] = old_vx[old_i * 3 + 1];
            new_vx[new_i * 3 + 2] = old_vx[old_i * 3 + 2];
#else
            new_x[new_i]  = old_x[old_i];
            new_y[new_i]  = old_y[old_i];
            new_z[new_i]  = old_z[old_i];
            new_vx[new_i] = old_vx[old_i];
            new_vy[new_i] = old_vy[old_i];
            new_vz[new_i] = old_vz[old_i];
#endif
        }
    }

    free(atom->x);
    free(atom->vx);
    atom->x  = new_x;
    atom->vx = new_vx;
#ifndef ATOM_POSITION_AOS
    free(atom->y);
    free(atom->z);
    free(atom->vy);
    free(atom->vz);
    atom->y  = new_y;
    atom->z  = new_z;
    atom->vy = new_vy;
    atom->vz = new_vz;
#endif

    DEBUG_MESSAGE("sortAtom end\n");
}

/* internal subroutines
Added with MPI*/

static int ghostZone(Atom* atom, int i) {
    if (i < atom->Nlocal) return 1;
    else if (method == halfShell)
        return halfZone(atom, i);
    else if (method == eightShell)
        return eightZone(atom, i);
    else
        return 0;
}

static int eightZone(Atom* atom, int i) {
    // Mapping: 0->0, 1->1, 2->2, 3->6, 4->3, 5->5, 6->4, 7->7
    int zoneMapping[] = { 0, 1, 2, 6, 3, 5, 4, 7 };
    MD_FLOAT* hi      = atom->mybox.hi;
    int zone          = 0;

    if (BigOrEqual(atom_x(i), hi[0])) {
        zone += 1;
    }

    if (BigOrEqual(atom_y(i), hi[1])) {
        zone += 2;
    }

    if (BigOrEqual(atom_z(i), hi[2])) {
        zone += 4;
    }

    return zoneMapping[zone];
}

static int halfZone(Atom* atom, int i) {
    MD_FLOAT* hi = atom->mybox.hi;
    MD_FLOAT* lo = atom->mybox.lo;

    if (atom_x(i) < lo[0] && atom_y(i) < hi[1] && atom_z(i) < hi[2]) {
        return 0;
    } else if (atom_y(i) < lo[1] && atom_z(i) < hi[2]) {
        return 0;
    } else if (atom_z(i) < lo[2]) {
        return 0;
    } else {
        return 1;
    }
}

static void neighborGhost(Atom* atom, Neighbor* neighbor) {
    DEBUG_MESSAGE("neighborGhost begin\n");

    int Nshell = 0;
    int Nlocal = atom->Nlocal;
    int Nghost = atom->Nghost;
    if (neighbor->listshell) free(neighbor->listshell);
    neighbor->listshell = (int*)malloc(Nghost * sizeof(int));
    int* listzone       = (int*)malloc(8 * Nghost * sizeof(int));
    int countAtoms[8]   = { 0, 0, 0, 0, 0, 0, 0, 0 };

    // Selecting ghost atoms for interaction
    for (int i = Nlocal; i < Nlocal + Nghost; i++) {
        int izone = ghostZone(atom, i);
        int* list = &listzone[Nghost * izone];
        int n     = countAtoms[izone];
        list[n]   = i;
        countAtoms[izone]++;
    }

    for (int zone = 1; zone <= 3; zone++) {
        int* list = &listzone[Nghost * zone];
        for (int n = 0; n < countAtoms[zone]; n++)
            neighbor->listshell[Nshell++] = list[n];
    }

    neighbor->Nshell = Nshell;
    if (neighbor->numNeighShell) free(neighbor->numNeighShell);
    if (neighbor->neighshell) free(neighbor->neighshell);
    neighbor->neighshell    = (int*)malloc(Nshell * neighbor->maxneighs * sizeof(int));
    neighbor->numNeighShell = (int*)malloc(Nshell * sizeof(int));
    int resize              = 1;

    while (resize) {
        resize = 0;
        for (int i = 0; i < Nshell; i++) {
            int n           = 0;
            int iatom       = neighbor->listshell[i];
            int izone       = ghostZone(atom, iatom);
            MD_FLOAT xtmp   = atom_x(iatom);
            MD_FLOAT ytmp   = atom_y(iatom);
            MD_FLOAT ztmp   = atom_z(iatom);
            int ibin        = coord2bin(xtmp, ytmp, ztmp);

#ifdef EXPLICIT_TYPES
            int type_i = atom->type[iatom];
#endif

            for (int k = 0; k < nstencil; k++) {
                int jbin     = ibin + stencil[k];
                int* loc_bin = &bins[jbin * atoms_per_bin];
                for (int m = 0; m < bincount[jbin]; m++) {
                    int jatom = loc_bin[m];

                    int jzone = ghostZone(atom, jatom);

                    if (jzone <= izone) continue;
                    if (izone == 1 && (jzone == 5 || jzone == 6 || jzone == 7)) continue;
                    if (izone == 2 && (jzone == 4 || jzone == 6 || jzone == 7)) continue;
                    if (izone == 3 && (jzone == 4 || jzone == 5 || jzone == 7)) continue;

                    MD_FLOAT delx = xtmp - atom_x(jatom);
                    MD_FLOAT dely = ytmp - atom_y(jatom);
                    MD_FLOAT delz = ztmp - atom_z(jatom);
                    MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;

#ifdef EXPLICIT_TYPES
                    int type_j = atom->type[jatom];
                    const MD_FLOAT cutoff =
                        atom->cutneighsq[type_i * atom->ntypes + type_j];
#else
                    const MD_FLOAT cutoff = cutneighsq;
#endif
                    if (rsq <= cutoff) {
                        neighs(neighbor->neighshell, i, n, Nshell, neighbor->maxneighs) = jatom;
                        n++;
                    }
                }
            }

            neighbor->numNeighShell[i] = n;
            if (n >= neighbor->maxneighs) {
                resize              = 1;
                neighbor->maxneighs = n * 1.2;
                break;
            }
        }

        if (resize) {
            free(neighbor->neighshell);
            neighbor->neighshell = (int*)malloc(
                Nshell * neighbor->maxneighs * sizeof(int));
        }
    }

    free(listzone);
    DEBUG_MESSAGE("neighborGhost end\n");
}

static inline int skipNeigh(Atom* atom, int i, int j) {
    if (i > j && j < atom->Nlocal) {
        return 1;
    } else if (atom_z(i) > atom_z(j) && j >= atom->Nlocal) {
        return 1;
    } else if (Equal(atom_z(i), atom_z(j)) && atom_y(i) > atom_y(j) &&
               j >= atom->Nlocal) {
        return 1;
    } else if (Equal(atom_z(i), atom_z(j)) && Equal(atom_y(i), atom_y(j)) &&
               atom_x(i) > atom_x(j) && j >= atom->Nlocal) {
        return 1;
    } else {
        return 0;
    }
}
