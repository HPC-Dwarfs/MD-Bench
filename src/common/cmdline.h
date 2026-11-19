/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#ifndef __CMDLINE_H_
#define __CMDLINE_H_

#include <comm.h>
#include <parameter.h>

// Command-line options of both schemes, validates them and derives cutneigh
void parseCommandLine(int argc, char** argv, Parameter* param, Comm* comm);

#endif
