/*
 * Copyright (C)  NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of MD-Bench.
 * Use of this source code is governed by a LGPL-3.0
 * license that can be found in the LICENSE file.
 */
#include <atom.h>
#include <comm.h>
#include <parameter.h>

#ifndef __VTK_H_
#define __VTK_H_
extern void write_data_to_vtk_file(const char* filename, Atom* atom, Parameter* param, int timestep);
extern int write_local_atoms_to_vtk_file(const char* filename, Atom* atom, Parameter* param, int timestep);
extern int write_ghost_atoms_to_vtk_file(const char* filename, Atom* atom, int timestep);
extern int write_local_cluster_edges_to_vtk_file(
    const char* filename, Atom* atom, int timestep);
extern int write_ghost_cluster_edges_to_vtk_file(
    const char* filename, Atom* atom, int timestep);
extern void printvtk(
    const char* filename, Comm* comm, Atom* atom, Parameter* param, int timestep);
extern int write_super_clusters_to_vtk_file(
    const char* filename, Atom* atom, int timestep);
#endif
