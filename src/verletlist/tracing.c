/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <neighbor.h>
#include <parameter.h>
#include <tracing.h>

void traceAddresses(Parameter* param, Atom* atom, Neighbor* neighbor, int timestep)
{
    MEM_TRACER_INIT;
    INDEX_TRACER_INIT;
    int Nlocal = atom->Nlocal;

    INDEX_TRACE_NATOMS(Nlocal, atom->Nghost, neighbor->maxneighs);
    for (int i = 0; i < Nlocal; i++) {
        int numneighs = neighbor->numneigh[i];
        MEM_TRACE(atom_x(i), 'R');
        MEM_TRACE(atom_y(i), 'R');
        MEM_TRACE(atom_z(i), 'R');
        INDEX_TRACE_ATOM(i);

#if LJ_COMB_RULE != LJ_COMB_SINGLE
        MEM_TRACE(atom->type[i], 'R');
#endif

        DIST_TRACE_SORT(neighbor->neighbors, i, numneighs, Nlocal, neighbor->maxneighs);
        INDEX_TRACE(neighbor->neighbors, i, numneighs, Nlocal, neighbor->maxneighs);
        DIST_TRACE(neighbor->neighbors, i, numneighs, Nlocal, neighbor->maxneighs);

        for (int k = 0; k < numneighs; k++) {
            MEM_TRACE(neighs(neighbor->neighbors, i, k, Nlocal, neighbor->maxneighs), 'R');
            MEM_TRACE(atom_x(j), 'R');
            MEM_TRACE(atom_y(j), 'R');
            MEM_TRACE(atom_z(j), 'R');

#if LJ_COMB_RULE != LJ_COMB_SINGLE
            MEM_TRACE(atom->type[j], 'R');
#endif
        }

        MEM_TRACE(atom_fx(i), 'R');
        MEM_TRACE(atom_fx(i), 'W');
        MEM_TRACE(atom_fy(i), 'R');
        MEM_TRACE(atom_fy(i), 'W');
        MEM_TRACE(atom_fz(i), 'R');
        MEM_TRACE(atom_fz(i), 'W');
    }

    INDEX_TRACER_END;
    MEM_TRACER_END;
}
