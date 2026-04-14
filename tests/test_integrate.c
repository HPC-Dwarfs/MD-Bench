#include "test_runner.h"

#include <stdlib.h>
#include <string.h>

#include <atom.h>
#include <force.h>
#include <integrate.h>
#include <parameter.h>

/*
 * Unit tests for the Velocity Verlet integrator (clusterpair variant).
 *
 * The integrator splits into two half-steps:
 *   initialIntegrateCPU: v += dtforce * f,  x += dt * v
 *   finalIntegrateCPU:   v += dtforce * f
 *
 * Together they implement a full Velocity Verlet time step
 * (dtforce = dt/2).
 */

/* Helper: allocate cluster arrays for nci i-clusters and zero them. */
static int alloc_clusters(Atom* atom, int nci, int atoms_per_cluster)
{
    /* Each cluster vector block needs 3 * max(CLUSTER_M, CLUSTER_N) floats. */
#if CLUSTER_M >= CLUSTER_N
    int cl_vec_size = 3 * CLUSTER_M;
#else
    int cl_vec_size = 3 * CLUSTER_N;
#endif
    int total = nci * cl_vec_size;

    atom->Nclusters_local = nci;
    atom->iclusters        = (Cluster*)calloc(nci, sizeof(Cluster));
    ASSERT_TRUE(atom->iclusters != NULL, "alloc iclusters");

    for (int ci = 0; ci < nci; ci++) {
        atom->iclusters[ci].natoms = atoms_per_cluster;
    }

    atom->cl_x = (MD_FLOAT*)calloc(total, sizeof(MD_FLOAT));
    atom->cl_v = (MD_FLOAT*)calloc(total, sizeof(MD_FLOAT));
    atom->cl_f = (MD_FLOAT*)calloc(total, sizeof(MD_FLOAT));
    ASSERT_TRUE(atom->cl_x && atom->cl_v && atom->cl_f, "alloc cl arrays");
    return 0;
}

static void free_clusters(Atom* atom)
{
    free(atom->iclusters);
    free(atom->cl_x);
    free(atom->cl_v);
    free(atom->cl_f);
}

/* Test 1: free particle (F=0).
 * With zero force, velocity is unchanged and position advances by dt*v.
 */
static int test_integrate_free_particle(void)
{
    Parameter param;
    initParameter(&param);
    param.dt      = 0.005;
    param.dtforce = 0.5 * param.dt;

    Atom atom;
    memset(&atom, 0, sizeof(Atom));
    alloc_clusters(&atom, 1, 1);

    int base = CI_VECTOR3_BASE_INDEX(0);
    MD_FLOAT vx0 = 1.0, vy0 = 2.0, vz0 = 3.0;
    MD_FLOAT x0 = 10.0, y0 = 20.0, z0 = 30.0;

    atom.cl_v[base + CL_X_INDEX_3D(0)] = vx0;
    atom.cl_v[base + CL_Y_INDEX_3D(0)] = vy0;
    atom.cl_v[base + CL_Z_INDEX_3D(0)] = vz0;
    atom.cl_x[base + CL_X_INDEX_3D(0)] = x0;
    atom.cl_x[base + CL_Y_INDEX_3D(0)] = y0;
    atom.cl_x[base + CL_Z_INDEX_3D(0)] = z0;
    /* cl_f already zeroed by calloc */

    initialIntegrateCPU(&param, &atom);

    /* Velocity unchanged (dtforce * 0 = 0) */
    ASSERT_NEAR(atom.cl_v[base + CL_X_INDEX_3D(0)], vx0, 1e-12,
        "free particle: vx unchanged");
    ASSERT_NEAR(atom.cl_v[base + CL_Y_INDEX_3D(0)], vy0, 1e-12,
        "free particle: vy unchanged");
    ASSERT_NEAR(atom.cl_v[base + CL_Z_INDEX_3D(0)], vz0, 1e-12,
        "free particle: vz unchanged");

    /* Position: x_new = x0 + dt * v */
    ASSERT_NEAR(atom.cl_x[base + CL_X_INDEX_3D(0)], x0 + param.dt * vx0, 1e-12,
        "free particle: x advanced");
    ASSERT_NEAR(atom.cl_x[base + CL_Y_INDEX_3D(0)], y0 + param.dt * vy0, 1e-12,
        "free particle: y advanced");
    ASSERT_NEAR(atom.cl_x[base + CL_Z_INDEX_3D(0)], z0 + param.dt * vz0, 1e-12,
        "free particle: z advanced");

    /* finalIntegrate with F=0 should also leave velocity unchanged */
    finalIntegrateCPU(&param, &atom);
    ASSERT_NEAR(atom.cl_v[base + CL_X_INDEX_3D(0)], vx0, 1e-12,
        "free particle final: vx unchanged");

    free_clusters(&atom);
    return 0;
}

/* Test 2: constant force.
 * After a full Velocity Verlet step (initial + final) with constant force:
 *   v(dt) = v0 + dt * a        (a = f/m, here m=1 is implicit)
 *   x(dt) = x0 + dt * v0 + 0.5 * dt^2 * a
 */
static int test_integrate_constant_force(void)
{
    Parameter param;
    initParameter(&param);
    param.dt      = 0.005;
    param.dtforce = 0.5 * param.dt;

    Atom atom;
    memset(&atom, 0, sizeof(Atom));
    alloc_clusters(&atom, 1, 1);

    int base   = CI_VECTOR3_BASE_INDEX(0);
    MD_FLOAT x0 = 1.0, y0 = 2.0, z0 = 3.0;
    MD_FLOAT vx0 = 0.5, vy0 = -0.3, vz0 = 0.7;
    MD_FLOAT fx = 2.0, fy = -1.0, fz = 0.5;

    atom.cl_x[base + CL_X_INDEX_3D(0)] = x0;
    atom.cl_x[base + CL_Y_INDEX_3D(0)] = y0;
    atom.cl_x[base + CL_Z_INDEX_3D(0)] = z0;
    atom.cl_v[base + CL_X_INDEX_3D(0)] = vx0;
    atom.cl_v[base + CL_Y_INDEX_3D(0)] = vy0;
    atom.cl_v[base + CL_Z_INDEX_3D(0)] = vz0;
    atom.cl_f[base + CL_X_INDEX_3D(0)] = fx;
    atom.cl_f[base + CL_Y_INDEX_3D(0)] = fy;
    atom.cl_f[base + CL_Z_INDEX_3D(0)] = fz;

    /* initialIntegrate: v += dtforce*f, then x += dt*v_new */
    initialIntegrateCPU(&param, &atom);

    MD_FLOAT vx_half = vx0 + param.dtforce * fx;
    MD_FLOAT vy_half = vy0 + param.dtforce * fy;
    MD_FLOAT vz_half = vz0 + param.dtforce * fz;

    ASSERT_NEAR(atom.cl_v[base + CL_X_INDEX_3D(0)], vx_half, 1e-12,
        "const force: vx after initial");
    ASSERT_NEAR(atom.cl_x[base + CL_X_INDEX_3D(0)], x0 + param.dt * vx_half, 1e-12,
        "const force: x after initial");
    ASSERT_NEAR(atom.cl_x[base + CL_Y_INDEX_3D(0)], y0 + param.dt * vy_half, 1e-12,
        "const force: y after initial");

    /* finalIntegrate: v += dtforce*f  (force unchanged for constant F) */
    finalIntegrateCPU(&param, &atom);

    MD_FLOAT vx_full = vx_half + param.dtforce * fx;
    MD_FLOAT vy_full = vy_half + param.dtforce * fy;
    MD_FLOAT vz_full = vz_half + param.dtforce * fz;

    ASSERT_NEAR(atom.cl_v[base + CL_X_INDEX_3D(0)], vx_full, 1e-12,
        "const force: vx after full step");
    ASSERT_NEAR(atom.cl_v[base + CL_Y_INDEX_3D(0)], vy_full, 1e-12,
        "const force: vy after full step");
    ASSERT_NEAR(atom.cl_v[base + CL_Z_INDEX_3D(0)], vz_full, 1e-12,
        "const force: vz after full step");

    /* v(dt) should equal v0 + dt*a */
    ASSERT_NEAR(vx_full, vx0 + param.dt * fx, 1e-12,
        "const force: v = v0 + dt*a (x)");
    ASSERT_NEAR(vy_full, vy0 + param.dt * fy, 1e-12,
        "const force: v = v0 + dt*a (y)");
    ASSERT_NEAR(vz_full, vz0 + param.dt * fz, 1e-12,
        "const force: v = v0 + dt*a (z)");

    /* x(dt) should equal x0 + dt*v0 + 0.5*dt^2*a */
    double dt2a_half = 0.5 * param.dt * param.dt;
    ASSERT_NEAR(atom.cl_x[base + CL_X_INDEX_3D(0)],
        x0 + param.dt * vx0 + dt2a_half * fx, 1e-12,
        "const force: x = x0 + dt*v0 + 0.5*dt^2*a (x)");
    ASSERT_NEAR(atom.cl_x[base + CL_Y_INDEX_3D(0)],
        y0 + param.dt * vy0 + dt2a_half * fy, 1e-12,
        "const force: x = x0 + dt*v0 + 0.5*dt^2*a (y)");
    ASSERT_NEAR(atom.cl_x[base + CL_Z_INDEX_3D(0)],
        z0 + param.dt * vz0 + dt2a_half * fz, 1e-12,
        "const force: x = x0 + dt*v0 + 0.5*dt^2*a (z)");

    free_clusters(&atom);
    return 0;
}

/* Test 3: time reversibility with constant force.
 * Apply forward step, negate velocity, apply another step.
 * Position should return to the original value and velocity
 * should be the negative of the original.
 */
static int test_integrate_time_reversible(void)
{
    Parameter param;
    initParameter(&param);
    param.dt      = 0.005;
    param.dtforce = 0.5 * param.dt;

    Atom atom;
    memset(&atom, 0, sizeof(Atom));
    alloc_clusters(&atom, 1, 1);

    int base   = CI_VECTOR3_BASE_INDEX(0);
    MD_FLOAT x0 = 5.0, y0 = 6.0, z0 = 7.0;
    MD_FLOAT vx0 = 1.2, vy0 = -0.8, vz0 = 0.4;
    MD_FLOAT fx = 3.0, fy = -2.0, fz = 1.5;

    atom.cl_x[base + CL_X_INDEX_3D(0)] = x0;
    atom.cl_x[base + CL_Y_INDEX_3D(0)] = y0;
    atom.cl_x[base + CL_Z_INDEX_3D(0)] = z0;
    atom.cl_v[base + CL_X_INDEX_3D(0)] = vx0;
    atom.cl_v[base + CL_Y_INDEX_3D(0)] = vy0;
    atom.cl_v[base + CL_Z_INDEX_3D(0)] = vz0;
    atom.cl_f[base + CL_X_INDEX_3D(0)] = fx;
    atom.cl_f[base + CL_Y_INDEX_3D(0)] = fy;
    atom.cl_f[base + CL_Z_INDEX_3D(0)] = fz;

    /* Forward step */
    initialIntegrateCPU(&param, &atom);
    finalIntegrateCPU(&param, &atom);

    /* Negate velocity */
    atom.cl_v[base + CL_X_INDEX_3D(0)] = -atom.cl_v[base + CL_X_INDEX_3D(0)];
    atom.cl_v[base + CL_Y_INDEX_3D(0)] = -atom.cl_v[base + CL_Y_INDEX_3D(0)];
    atom.cl_v[base + CL_Z_INDEX_3D(0)] = -atom.cl_v[base + CL_Z_INDEX_3D(0)];

    /* Reverse step (same force since it's constant) */
    initialIntegrateCPU(&param, &atom);
    finalIntegrateCPU(&param, &atom);

    /* Position should return to x0 */
    ASSERT_NEAR(atom.cl_x[base + CL_X_INDEX_3D(0)], x0, 1e-10,
        "reversible: x returned");
    ASSERT_NEAR(atom.cl_x[base + CL_Y_INDEX_3D(0)], y0, 1e-10,
        "reversible: y returned");
    ASSERT_NEAR(atom.cl_x[base + CL_Z_INDEX_3D(0)], z0, 1e-10,
        "reversible: z returned");

    /* Velocity should be -v0 */
    ASSERT_NEAR(atom.cl_v[base + CL_X_INDEX_3D(0)], -vx0, 1e-10,
        "reversible: vx negated");
    ASSERT_NEAR(atom.cl_v[base + CL_Y_INDEX_3D(0)], -vy0, 1e-10,
        "reversible: vy negated");
    ASSERT_NEAR(atom.cl_v[base + CL_Z_INDEX_3D(0)], -vz0, 1e-10,
        "reversible: vz negated");

    free_clusters(&atom);
    return 0;
}

int run_integrate_tests(void)
{
    int rc = 0;

    tr_log("  integrate: free particle (F=0)");
    rc = test_integrate_free_particle();
    if (rc) return rc;

    tr_log("  integrate: constant force full step");
    rc = test_integrate_constant_force();
    if (rc) return rc;

    tr_log("  integrate: time reversibility");
    rc = test_integrate_time_reversible();
    if (rc) return rc;

    return 0;
}
