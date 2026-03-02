#include "test_runner.h"

#include <math.h>

/* Mirror the Lennard-Jones force calculation used in computeForceLJRef. */
static void lj_force(
    double epsilon, double sigma6, double cutforcesq, double dx, double dy, double dz, double* fx, double* fy, double* fz)
{
    double rsq = dx * dx + dy * dy + dz * dz;
    if (rsq >= cutforcesq) {
        *fx = 0.0;
        *fy = 0.0;
        *fz = 0.0;
        return;
    }

    double sr2   = 1.0 / rsq;
    double sr6   = sr2 * sr2 * sr2 * sigma6;
    double force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;

    *fx = dx * force;
    *fy = dy * force;
    *fz = dz * force;
}

static int test_lj_zero_force_at_minimum(void)
{
    const double epsilon    = 1.0;
    const double sigma      = 1.0;
    const double sigma2     = sigma * sigma;
    const double sigma6     = sigma2 * sigma2 * sigma2;
    const double cutforcesq = 100.0; /* large enough that cutoff does not apply */

    /* At r = 2^(1/6) * sigma, the LJ force should be zero. */
    const double r = pow(2.0, 1.0 / 6.0) * sigma;
    double fx, fy, fz;
    lj_force(epsilon, sigma6, cutforcesq, r, 0.0, 0.0, &fx, &fy, &fz);

    ASSERT_NEAR(fx, 0.0, 1e-10, "LJ force x-component at minimum");
    ASSERT_NEAR(fy, 0.0, 1e-10, "LJ force y-component at minimum");
    ASSERT_NEAR(fz, 0.0, 1e-10, "LJ force z-component at minimum");
    return 0;
}

static int test_lj_newtons_third_law(void)
{
    const double epsilon    = 1.0;
    const double sigma      = 1.0;
    const double sigma2     = sigma * sigma;
    const double sigma6     = sigma2 * sigma2 * sigma2;
    const double cutforcesq = 100.0;

    const double r = 1.3 * sigma;
    double fi_x, fi_y, fi_z;
    double fj_x, fj_y, fj_z;

    /* Force on i due to j at +r */
    lj_force(epsilon, sigma6, cutforcesq, r, 0.0, 0.0, &fi_x, &fi_y, &fi_z);
    /* Force on j due to i at -r */
    lj_force(epsilon, sigma6, cutforcesq, -r, 0.0, 0.0, &fj_x, &fj_y, &fj_z);

    ASSERT_NEAR(fi_x + fj_x, 0.0, 1e-12, "Newton 3rd law (x)");
    ASSERT_NEAR(fi_y + fj_y, 0.0, 1e-12, "Newton 3rd law (y)");
    ASSERT_NEAR(fi_z + fj_z, 0.0, 1e-12, "Newton 3rd law (z)");
    return 0;
}

static int test_lj_cutoff_gating(void)
{
    const double epsilon    = 1.0;
    const double sigma      = 1.0;
    const double sigma2     = sigma * sigma;
    const double sigma6     = sigma2 * sigma2 * sigma2;
    const double rcut       = 2.5 * sigma;
    const double cutforcesq = rcut * rcut;

    double fx_in, fy_in, fz_in;
    double fx_out, fy_out, fz_out;

    /* Slightly inside the cutoff */
    lj_force(epsilon, sigma6, cutforcesq, 0.99 * rcut, 0.0, 0.0, &fx_in, &fy_in, &fz_in);
    /* Slightly outside the cutoff */
    lj_force(epsilon, sigma6, cutforcesq, 1.01 * rcut, 0.0, 0.0, &fx_out, &fy_out, &fz_out);

    ASSERT_TRUE(fabs(fx_in) > 0.0, "non-zero force inside cutoff");
    ASSERT_NEAR(fx_out, 0.0, 1e-15, "zero force outside cutoff (x)");
    ASSERT_NEAR(fy_out, 0.0, 1e-15, "zero force outside cutoff (y)");
    ASSERT_NEAR(fz_out, 0.0, 1e-15, "zero force outside cutoff (z)");
    return 0;
}

int run_force_tests(void)
{
    int rc = 0;

    tr_log("  force: LJ zero at minimum");
    rc = test_lj_zero_force_at_minimum();
    if (rc)
        return rc;

    tr_log("  force: Newton's third law");
    rc = test_lj_newtons_third_law();
    if (rc)
        return rc;

    tr_log("  force: cutoff gating");
    rc = test_lj_cutoff_gating();
    if (rc)
        return rc;

    return 0;
}

