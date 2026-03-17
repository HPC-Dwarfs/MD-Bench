/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <likwid-marker.h>
#include <math.h>

#include <allocate.h>
#include <atom.h>
#include <eam.h>
#include <neighbor.h>
#include <parameter.h>
#include <stats.h>
#include <timing.h>
#include <util.h>

double computeForceEam(Parameter* param, Atom* atom, Neighbor* neighbor, Stats* stats)
{
    if (eam.nmax < atom->Nmax) {
        eam.nmax = atom->Nmax;
        if (eam.fp != NULL) {
            free(eam.fp);
        }
        eam.fp = (MD_FLOAT*)allocate(ALIGNMENT, atom->Nmax * sizeof(MD_FLOAT));
    }

    int Nlocal = atom->Nlocal;
    int ntypes            = atom->ntypes;
    MD_FLOAT* fp          = eam.fp;
    MD_FLOAT* rhor_spline = eam.rhor_spline;
    MD_FLOAT* frho_spline = eam.frho_spline;
    MD_FLOAT* z2r_spline  = eam.z2r_spline;
    MD_FLOAT rdr          = eam.rdr;
    int nr                = eam.nr;
    int nr_tot            = eam.nr_tot;
    MD_FLOAT rdrho        = eam.rdrho;
    int nrho              = eam.nrho;
    int nrho_tot          = eam.nrho_tot;
    double timeStart      = getTimeStamp();

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for
        for (int i = 0; i < Nlocal; i++) {
            int numneighs = neighbor->numneigh[i];
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            MD_FLOAT rhoi = 0;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            const int type_i = atom->type[i];
#endif
#pragma ivdep
            for (int k = 0; k < numneighs; k++) {
                int j         = neighs(neighbor->neighbors, i, k, Nlocal, neighbor->maxneighs);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                const int type_j          = atom->type[j];
                const int type_ij         = type_i * ntypes + type_j;
                const MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
#else
                const MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif
                if (rsq < cutforcesq) {
                    MD_FLOAT p = sqrt(rsq) * rdr + 1.0;
                    int m      = (int)(p);
                    m          = m < nr - 1 ? m : nr - 1;
                    p -= m;
                    p = p < 1.0 ? p : 1.0;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                    rhoi += ((rhor_spline[type_ij * nr_tot + m * 7 + 3] * p +
                                 rhor_spline[type_ij * nr_tot + m * 7 + 4]) *
                                    p +
                                rhor_spline[type_ij * nr_tot + m * 7 + 5]) *
                                p +
                            rhor_spline[type_ij * nr_tot + m * 7 + 6];
#else
                    rhoi += ((rhor_spline[m * 7 + 3] * p + rhor_spline[m * 7 + 4]) * p +
                                rhor_spline[m * 7 + 5]) *
                                p +
                            rhor_spline[m * 7 + 6];
#endif
                }
            }

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            const int type_ii = type_i * type_i;
#endif
            MD_FLOAT p = 1.0 * rhoi * rdrho + 1.0;
            int m      = (int)(p);
            m          = MAX(1, MIN(m, nrho - 1));
            p -= m;
            p = MIN(p, 1.0);
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            fp[i] = (frho_spline[type_ii * nrho_tot + m * 7 + 0] * p +
                        frho_spline[type_ii * nrho_tot + m * 7 + 1]) *
                        p +
                    frho_spline[type_ii * nrho_tot + m * 7 + 2];
#else
            fp[i] = (frho_spline[m * 7 + 0] * p + frho_spline[m * 7 + 1]) * p +
                    frho_spline[m * 7 + 2];
#endif
        }

        LIKWID_MARKER_STOP("force");
    }

    // We still need to update fp for PBC atoms
    for (int i = 0; i < atom->Nghost; i++) {
        fp[Nlocal + i] = fp[atom->border_map[i]];
    }

#pragma omp parallel
    {
        LIKWID_MARKER_START("force");

#pragma omp for
        for (int i = 0; i < Nlocal; i++) {
            int numneighs = neighbor->numneigh[i];
            MD_FLOAT xtmp = atom_x(i);
            MD_FLOAT ytmp = atom_y(i);
            MD_FLOAT ztmp = atom_z(i);
            MD_FLOAT fix  = 0;
            MD_FLOAT fiy  = 0;
            MD_FLOAT fiz  = 0;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
            const int type_i = atom->type[i];
#endif

#pragma ivdep
            for (int k = 0; k < numneighs; k++) {
                int j         = neighs(neighbor->neighbors, i, k, Nlocal, neighbor->maxneighs);
                MD_FLOAT delx = xtmp - atom_x(j);
                MD_FLOAT dely = ytmp - atom_y(j);
                MD_FLOAT delz = ztmp - atom_z(j);
                MD_FLOAT rsq  = delx * delx + dely * dely + delz * delz;
#if LJ_COMB_RULE != LJ_COMB_SINGLE
                const int type_j          = atom->type[j];
                const int type_ij         = type_i * ntypes + type_j;
                const MD_FLOAT cutforcesq = atom->cutforcesq[type_ij];
#else
                const MD_FLOAT cutforcesq = param->cutforce * param->cutforce;
#endif

                if (rsq < cutforcesq) {
                    MD_FLOAT r = sqrt(rsq);
                    MD_FLOAT p = r * rdr + 1.0;
                    int m      = (int)(p);
                    m          = m < nr - 1 ? m : nr - 1;
                    p -= m;
                    p = p < 1.0 ? p : 1.0;

                    // rhoip = derivative of (density at atom j due to atom i)
                    // rhojp = derivative of (density at atom i due to atom j)
                    // phi = pair potential energy
                    // phip = phi'
                    // z2 = phi * r
                    // z2p = (phi * r)' = (phi' r) + phi
                    // psip needs both fp[i] and fp[j] terms since r_ij appears in two
                    //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
                    //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

#if LJ_COMB_RULE != LJ_COMB_SINGLE
                    MD_FLOAT rhoip = (rhor_spline[type_ij * nr_tot + m * 7 + 0] * p +
                                         rhor_spline[type_ij * nr_tot + m * 7 + 1]) *
                                         p +
                                     rhor_spline[type_ij * nr_tot + m * 7 + 2];

                    MD_FLOAT z2p = (z2r_spline[type_ij * nr_tot + m * 7 + 0] * p +
                                       z2r_spline[type_ij * nr_tot + m * 7 + 1]) *
                                       p +
                                   z2r_spline[type_ij * nr_tot + m * 7 + 2];

                    MD_FLOAT z2 = ((z2r_spline[type_ij * nr_tot + m * 7 + 3] * p +
                                       z2r_spline[type_ij * nr_tot + m * 7 + 4]) *
                                          p +
                                      z2r_spline[type_ij * nr_tot + m * 7 + 5]) *
                                      p +
                                  z2r_spline[type_ij * nr_tot + m * 7 + 6];
#else
                    MD_FLOAT rhoip = (rhor_spline[m * 7 + 0] * p +
                                         rhor_spline[m * 7 + 1]) *
                                         p +
                                     rhor_spline[m * 7 + 2];
                    MD_FLOAT z2p = (z2r_spline[m * 7 + 0] * p + z2r_spline[m * 7 + 1]) *
                                       p +
                                   z2r_spline[m * 7 + 2];
                    MD_FLOAT z2 = ((z2r_spline[m * 7 + 3] * p + z2r_spline[m * 7 + 4]) *
                                          p +
                                      z2r_spline[m * 7 + 5]) *
                                      p +
                                  z2r_spline[m * 7 + 6];
#endif

                    MD_FLOAT recip = 1.0 / r;
                    MD_FLOAT phi   = z2 * recip;
                    MD_FLOAT phip  = z2p * recip - phi * recip;
                    MD_FLOAT psip  = fp[i] * rhoip + fp[j] * rhoip + phip;
                    MD_FLOAT fpair = -psip * recip;

                    fix += delx * fpair;
                    fiy += dely * fpair;
                    fiz += delz * fpair;
                    // fpair *= 0.5;
                }
            }

            atom_fx(i) = fix;
            atom_fy(i) = fiy;
            atom_fz(i) = fiz;
            addStat(stats->total_force_neighs, numneighs);
            addStat(stats->total_force_iters,
                (numneighs + VECTOR_WIDTH - 1) / VECTOR_WIDTH);
        }

        LIKWID_MARKER_STOP("force");
    }

    double timeStop = getTimeStamp();
    return timeStop - timeStart;
}
