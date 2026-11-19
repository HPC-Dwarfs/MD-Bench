/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
// Tabulated Ewald correction erf(beta r)/r (COULOMB_RULE=ewald_table), as GROMACS' tables
#ifndef __EWALD_TABLE_H_
#define __EWALD_TABLE_H_

#include <parameter.h>

// Builds the tables for param->ewald_beta and param->cutforce
void initEwaldTable(Parameter* param);

#endif
