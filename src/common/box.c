/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <box.h>
#include <parameter.h>
#include <stdio.h>
#include <util.h>

int overlapBox(int dim,
    int dir,
    const Box* mybox,
    const Box* other,
    Box* cut,
    MD_FLOAT xprd,
    MD_FLOAT cutneigh)
{
    int pbc = -100;
    MD_FLOAT min[3], max[3];
    int same = (mybox->id == other->id) ? 1 : 0;

    // projections
    min[0] = MAX(mybox->lo[0], other->lo[0]);
    max[0] = MIN(mybox->hi[0], other->hi[0]);
    min[1] = MAX(mybox->lo[1], other->lo[1]);
    max[1] = MIN(mybox->hi[1], other->hi[1]);
    min[2] = MAX(mybox->lo[2], other->lo[2]);
    max[2] = MIN(mybox->hi[2], other->hi[2]);

    // Intersection no periodic case
    if (!same) {
        if (dir == 0) max[dim] = MIN(mybox->hi[dim], other->hi[dim] + cutneigh);
        if (dir == 1) min[dim] = MAX(mybox->lo[dim], other->lo[dim] - cutneigh);
        if ((min[0] < max[0]) && (min[1] < max[1]) && (min[2] < max[2])) pbc = 0;
    }

    // Intersection periodic case
    if (pbc < 0) {
        if (dir == 0) {
            min[dim] = MAX(mybox->lo[dim], other->lo[dim] - xprd);
            max[dim] = MIN(mybox->hi[dim], other->hi[dim] - xprd + cutneigh);

        } else {
            min[dim] = MAX(mybox->lo[dim], other->lo[dim] + xprd - cutneigh);
            max[dim] = MIN(mybox->hi[dim], other->hi[dim] + xprd);
        }
        if ((min[0] < max[0]) && (min[1] < max[1]) && (min[2] < max[2]))
            pbc = (dir == 0) ? 1 : -1;
    }

    // storing the cuts
    cut->lo[0] = min[0];
    cut->hi[0] = max[0];
    cut->lo[1] = min[1];
    cut->hi[1] = max[1];
    cut->lo[2] = min[2];
    cut->hi[2] = max[2];

    return pbc;
}

int overlapFullBox(
    Parameter* param, MD_FLOAT* cutneigh, const Box* mybox, const Box* other)
{
    MD_FLOAT min[3], max[3];
    MD_FLOAT xprd = param->xprd;
    MD_FLOAT yprd = param->yprd;
    MD_FLOAT zprd = param->zprd;

    for (int k = -1; k < 2; k++) {
        for (int j = -1; j < 2; j++) {
            for (int i = -1; i < 2; i++) {
                min[0] = MAX(mybox->lo[0], other->lo[0] - cutneigh[0] + i * xprd);
                min[1] = MAX(mybox->lo[1], other->lo[1] - cutneigh[1] + j * yprd);
                min[2] = MAX(mybox->lo[2], other->lo[2] - cutneigh[2] + k * zprd);
                max[0] = MIN(mybox->hi[0], other->hi[0] + cutneigh[0] + i * xprd);
                max[1] = MIN(mybox->hi[1], other->hi[1] + cutneigh[1] + j * yprd);
                max[2] = MIN(mybox->hi[2], other->hi[2] + cutneigh[2] + k * zprd);
                if ((min[0] < max[0]) && (min[1] < max[1]) && (min[2] < max[2])) return 1;
            }
        }
    }

    return 0;
}

void expandBox(int iswap, const Box* me, const Box* other, Box* cut, MD_FLOAT cutneigh)
{
    if (iswap == 2 || iswap == 3) {
        if (me->lo[0] <= other->lo[0]) cut->lo[0] -= cutneigh;
        if (me->hi[0] >= other->hi[0]) cut->hi[0] += cutneigh;
    }

    if (iswap == 4 || iswap == 5) {
        if (me->lo[0] <= other->lo[0]) cut->lo[0] -= cutneigh;
        if (me->hi[0] >= other->hi[0]) cut->hi[0] += cutneigh;
        if (me->lo[1] <= other->lo[1]) cut->lo[1] -= cutneigh;
        if (me->hi[1] >= other->hi[1]) cut->hi[1] += cutneigh;
    }
}
