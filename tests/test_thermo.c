#include "test_runner.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <atom.h>
#include <force.h>
#include <parameter.h>
#include <thermo.h>

/*
 * Unit tests for thermodynamic routines (adjustThermo).
 *
 * adjustThermo performs two operations:
 *   1. Zeroes center-of-mass velocity
 *   2. Rescales all velocities so the kinetic temperature matches param->temp
 */

/* Helper: set up a minimal atom with per-atom velocity arrays. */
static void setup_atom_velocities(Atom* atom, int natoms, MD_FLOAT* vx, MD_FLOAT* vy,
    MD_FLOAT* vz)
{
    memset(atom, 0, sizeof(Atom));
    atom->Natoms = natoms;
    atom->Nlocal = natoms;
    atom->Nmax   = natoms;

    atom->vx = (MD_FLOAT*)malloc(natoms * sizeof(MD_FLOAT));
    atom->vy = (MD_FLOAT*)malloc(natoms * sizeof(MD_FLOAT));
    atom->vz = (MD_FLOAT*)malloc(natoms * sizeof(MD_FLOAT));

    for (int i = 0; i < natoms; i++) {
        atom->vx[i] = vx[i];
        atom->vy[i] = vy[i];
        atom->vz[i] = vz[i];
    }

    /* Allocate position arrays (used by atom_x macros but not by adjustThermo).
     * Needed in case some include-path macro references them. */
#ifdef ATOM_POSITION_AOS
    atom->x = (MD_FLOAT*)calloc(natoms * 3, sizeof(MD_FLOAT));
#else
    atom->x = (MD_FLOAT*)calloc(natoms, sizeof(MD_FLOAT));
    atom->y = (MD_FLOAT*)calloc(natoms, sizeof(MD_FLOAT));
    atom->z = (MD_FLOAT*)calloc(natoms, sizeof(MD_FLOAT));
#endif
}

static void free_atom_velocities(Atom* atom)
{
    free(atom->vx);
    free(atom->vy);
    free(atom->vz);
    free(atom->x);
#ifndef ATOM_POSITION_AOS
    free(atom->y);
    free(atom->z);
#endif
}

/* Test 1: adjustThermo zeroes center-of-mass velocity.
 * After adjustThermo, the mean velocity in each dimension should be ~0.
 */
static int test_adjustThermo_zeroes_com(void)
{
    const int N = 4;
    MD_FLOAT vx[4] = { 1.0, 2.0, 3.0, 4.0 };    /* COM_x = 2.5 */
    MD_FLOAT vy[4] = { 0.5, 1.5, -0.5, 2.5 };    /* COM_y = 1.0 */
    MD_FLOAT vz[4] = { -1.0, 0.0, 1.0, 2.0 };    /* COM_z = 0.5 */

    Atom atom;
    setup_atom_velocities(&atom, N, vx, vy, vz);

    Parameter param;
    initParameter(&param);
    param.force_field = FF_LJ;
    param.temp        = 1.0;
    param.mass        = 1.0;
    param.xprd        = 10.0;
    param.yprd        = 10.0;
    param.zprd        = 10.0;

    setupThermo(&param, N);
    adjustThermo(&param, &atom);

    /* Check COM velocity is zero */
    MD_FLOAT com_x = 0.0, com_y = 0.0, com_z = 0.0;
    for (int i = 0; i < N; i++) {
        com_x += atom.vx[i];
        com_y += atom.vy[i];
        com_z += atom.vz[i];
    }
    com_x /= N;
    com_y /= N;
    com_z /= N;

    ASSERT_NEAR(com_x, 0.0, 1e-10, "COM vx zeroed");
    ASSERT_NEAR(com_y, 0.0, 1e-10, "COM vy zeroed");
    ASSERT_NEAR(com_z, 0.0, 1e-10, "COM vz zeroed");

    free_atom_velocities(&atom);
    return 0;
}

/* Test 2: adjustThermo rescales to target temperature.
 * After adjustThermo, the kinetic temperature should equal param->temp.
 *
 * For LJ units: T = sum(m*v^2) / (3N-3)
 * (with mvv2e=1.0 and mass=1.0)
 */
static int test_adjustThermo_rescales_temp(void)
{
    const int N = 4;
    MD_FLOAT vx[4] = { 1.0, 2.0, 3.0, 4.0 };
    MD_FLOAT vy[4] = { 0.5, 1.5, -0.5, 2.5 };
    MD_FLOAT vz[4] = { -1.0, 0.0, 1.0, 2.0 };

    Atom atom;
    setup_atom_velocities(&atom, N, vx, vy, vz);

    Parameter param;
    initParameter(&param);
    param.force_field = FF_LJ;
    param.temp        = 2.0;
    param.mass        = 1.0;
    param.xprd        = 10.0;
    param.yprd        = 10.0;
    param.zprd        = 10.0;

    setupThermo(&param, N);
    adjustThermo(&param, &atom);

    /* Compute kinetic temperature: T = sum(m*v^2) * t_scale
     * For LJ: t_scale = 1.0 / (3*N - 3) */
    double t_scale = 1.0 / (3.0 * N - 3.0);
    double ke      = 0.0;
    for (int i = 0; i < N; i++) {
        ke += (atom.vx[i] * atom.vx[i] + atom.vy[i] * atom.vy[i] +
                  atom.vz[i] * atom.vz[i]) *
              param.mass;
    }
    double T = ke * t_scale;

    ASSERT_NEAR(T, param.temp, 1e-10, "temperature rescaled to target");

    free_atom_velocities(&atom);
    return 0;
}

/* Test 3: adjustThermo preserves velocity ratios after rescaling.
 * The relative structure of velocities (after COM removal) should be
 * unchanged by the uniform rescaling factor.
 */
static int test_adjustThermo_preserves_structure(void)
{
    const int N = 3;
    MD_FLOAT vx[3] = { 3.0, 6.0, 9.0 };
    MD_FLOAT vy[3] = { 0.0, 0.0, 0.0 };
    MD_FLOAT vz[3] = { 0.0, 0.0, 0.0 };

    Atom atom;
    setup_atom_velocities(&atom, N, vx, vy, vz);

    Parameter param;
    initParameter(&param);
    param.force_field = FF_LJ;
    param.temp        = 1.0;
    param.mass        = 1.0;
    param.xprd        = 10.0;
    param.yprd        = 10.0;
    param.zprd        = 10.0;

    setupThermo(&param, N);
    adjustThermo(&param, &atom);

    /* After COM removal: vx = {-3, 0, 3} (COM was 6.0).
     * Ratios should be preserved: vx[2] = -vx[0], vx[1] = 0. */
    ASSERT_NEAR(atom.vx[0] + atom.vx[2], 0.0, 1e-10,
        "antisymmetric velocities preserved");
    ASSERT_NEAR(atom.vx[1], 0.0, 1e-10, "zero velocity preserved");

    free_atom_velocities(&atom);
    return 0;
}

int run_thermo_tests(void)
{
    int rc = 0;

    tr_log("  thermo: adjustThermo zeroes COM velocity");
    rc = test_adjustThermo_zeroes_com();
    if (rc) return rc;

    tr_log("  thermo: adjustThermo rescales to target temperature");
    rc = test_adjustThermo_rescales_temp();
    if (rc) return rc;

    tr_log("  thermo: adjustThermo preserves velocity structure");
    rc = test_adjustThermo_preserves_structure();
    if (rc) return rc;

    return 0;
}
