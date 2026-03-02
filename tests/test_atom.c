#include "test_runner.h"

#include <atom.h>

/* Forward declaration for helper that is only defined in atom.c. */
int typeStr2int(const char* type);

/* Basic unit tests for selected atom/clusterpair helpers. These are designed
 * to avoid pulling in the full MD-Bench dependency graph, so they only touch
 * functions that do not require MPI or neighbor structures.
 */

static int test_typeStr2int_valid_and_error(void)
{
    int t = typeStr2int("Ar");
    ASSERT_INT_EQ(t, 0, "typeStr2int(Ar) == 0");

    /* NOTE: Invalid types cause the program to exit in the production code,
     * so we do not test that path here to keep the unit tests simple.
     */
    return 0;
}

static int test_pbc_wraps_positions(void)
{
    Atom atom_storage;
    Atom* atom = &atom_storage;
    initAtom(atom);

    atom->Nlocal     = 2;
    atom->Nmax       = 2;
    atom->mybox.xprd = 10.0;
    atom->mybox.yprd = 10.0;
    atom->mybox.zprd = 10.0;

    /* Allocate minimal position/velocity/type arrays to satisfy macros. */
#ifdef ATOM_POSITION_AOS
    atom->x = (MD_FLOAT*)malloc(atom->Nmax * 3 * sizeof(MD_FLOAT));
    ASSERT_TRUE(atom->x != NULL, "alloc x (AOS)");
#else
    atom->x = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    atom->y = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    atom->z = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    ASSERT_TRUE(atom->x != NULL && atom->y != NULL && atom->z != NULL, "alloc x,y,z (SOA)");
#endif
    atom->vx   = (MD_FLOAT*)calloc(atom->Nmax, sizeof(MD_FLOAT));
    atom->vy   = (MD_FLOAT*)calloc(atom->Nmax, sizeof(MD_FLOAT));
    atom->vz   = (MD_FLOAT*)calloc(atom->Nmax, sizeof(MD_FLOAT));
    atom->type = (int*)calloc(atom->Nmax, sizeof(int));
    ASSERT_TRUE(atom->vx && atom->vy && atom->vz && atom->type, "alloc velocities/types");

    /* Put one atom slightly below 0 and one slightly above box size. */
    atom_x(0) = -0.1;
    atom_y(0) = -0.2;
    atom_z(0) = -0.3;

    atom_x(1) = 10.1;
    atom_y(1) = 10.2;
    atom_z(1) = 10.3;

    pbc(atom);

    ASSERT_NEAR(atom_x(0), 9.9, 1e-12, "pbc x < 0");
    ASSERT_NEAR(atom_y(0), 9.8, 1e-12, "pbc y < 0");
    ASSERT_NEAR(atom_z(0), 9.7, 1e-12, "pbc z < 0");

    ASSERT_NEAR(atom_x(1), 0.1, 1e-12, "pbc x >= box");
    ASSERT_NEAR(atom_y(1), 0.2, 1e-12, "pbc y >= box");
    ASSERT_NEAR(atom_z(1), 0.3, 1e-12, "pbc z >= box");

    free(atom->x);
#ifndef ATOM_POSITION_AOS
    free(atom->y);
    free(atom->z);
#endif
    free(atom->vx);
    free(atom->vy);
    free(atom->vz);
    free(atom->type);

    return 0;
}

static int test_pack_unpack_exchange_roundtrip(void)
{
    Atom atom_storage;
    Atom* atom = &atom_storage;
    initAtom(atom);
    atom->Nlocal = 1;
    atom->Nmax   = 1;

#ifdef ATOM_POSITION_AOS
    atom->x = (MD_FLOAT*)malloc(atom->Nmax * 3 * sizeof(MD_FLOAT));
    ASSERT_TRUE(atom->x != NULL, "alloc x (AOS)");
#else
    atom->x = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    atom->y = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    atom->z = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    ASSERT_TRUE(atom->x != NULL && atom->y != NULL && atom->z != NULL, "alloc x,y,z (SOA)");
#endif
    atom->vx   = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    atom->vy   = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    atom->vz   = (MD_FLOAT*)malloc(atom->Nmax * sizeof(MD_FLOAT));
    atom->type = (int*)malloc(atom->Nmax * sizeof(int));

    ASSERT_TRUE(atom->vx && atom->vy && atom->vz && atom->type, "alloc velocities/types");

    /* Original values. */
    atom_x(0)    = 1.0;
    atom_y(0)    = 2.0;
    atom_z(0)    = 3.0;
    atom_vx(0)   = 0.1;
    atom_vy(0)   = 0.2;
    atom_vz(0)   = 0.3;
    atom->type[0] = 7;

    MD_FLOAT buf[7];
    int m1 = packExchange(atom, 0, buf);

    /* Clear atom and unpack. */
    atom_x(0)    = 0.0;
    atom_y(0)    = 0.0;
    atom_z(0)    = 0.0;
    atom_vx(0)   = 0.0;
    atom_vy(0)   = 0.0;
    atom_vz(0)   = 0.0;
    atom->type[0] = 0;

    int m2 = unpackExchange(atom, 0, buf);

    ASSERT_INT_EQ(m1, 7, "packExchange element count");
    ASSERT_INT_EQ(m2, 7, "unpackExchange element count");

    ASSERT_NEAR(atom_x(0), 1.0, 1e-12, "x roundtrip");
    ASSERT_NEAR(atom_y(0), 2.0, 1e-12, "y roundtrip");
    ASSERT_NEAR(atom_z(0), 3.0, 1e-12, "z roundtrip");
    ASSERT_NEAR(atom_vx(0), 0.1, 1e-12, "vx roundtrip");
    ASSERT_NEAR(atom_vy(0), 0.2, 1e-12, "vy roundtrip");
    ASSERT_NEAR(atom_vz(0), 0.3, 1e-12, "vz roundtrip");
    ASSERT_INT_EQ(atom->type[0], 7, "type roundtrip");

    free(atom->x);
#ifndef ATOM_POSITION_AOS
    free(atom->y);
    free(atom->z);
#endif
    free(atom->vx);
    free(atom->vy);
    free(atom->vz);
    free(atom->type);

    return 0;
}

int run_atom_tests(void)
{
    int rc = 0;

    tr_log("  atom: typeStr2int");
    rc = test_typeStr2int_valid_and_error();
    if (rc)
        return rc;

    tr_log("  atom: pbc wraps positions");
    rc = test_pbc_wraps_positions();
    if (rc)
        return rc;

    tr_log("  atom: pack/unpack exchange roundtrip");
    rc = test_pack_unpack_exchange_roundtrip();
    if (rc)
        return rc;

    return 0;
}

