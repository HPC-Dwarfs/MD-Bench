/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#ifndef __TIMERS_H_
#define __TIMERS_H_

// typedef enum { TOTAL = 0, NEIGH, FORCE, NUMTIMER } timertype;

typedef enum {
    TOTAL = 0,
    FORCE,
    NEIGH,
    FORWARD,
    REVERSE,
    UPDATE,
    BALANCE,
    SETUP,
    REST,
    NUMTIMER
} timertype;

#endif
