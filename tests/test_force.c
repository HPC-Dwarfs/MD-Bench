#include "test_runner.h"

#include <math.h>
#include <string.h>

#include <ljtable.h>
#include <parameter.h>

/* Mirror the Lennard-Jones force calculation used in computeForceLJRef. */
static void lj_force(double epsilon,
    double sigma6,
    double cutforcesq,
    double dx,
    double dy,
    double dz,
    double* fx,
    double* fy,
    double* fz)
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
    lj_force(epsilon,
        sigma6,
        cutforcesq,
        1.01 * rcut,
        0.0,
        0.0,
        &fx_out,
        &fy_out,
        &fz_out);

    ASSERT_TRUE(fabs(fx_in) > 0.0, "non-zero force inside cutoff");
    ASSERT_NEAR(fx_out, 0.0, 1e-15, "zero force outside cutoff (x)");
    ASSERT_NEAR(fy_out, 0.0, 1e-15, "zero force outside cutoff (y)");
    ASSERT_NEAR(fz_out, 0.0, 1e-15, "zero force outside cutoff (z)");
    return 0;
}

static int test_lj_geom_combination_formula(void)
{
    const double cutforcesq = 100.0;
    const double r          = 1.3;

    /* Type A: eps=1.0, sig=1.0.  Type B: eps=4.0, sig=2.0. */
    const double eps_A = 1.0, sig_A = 1.0;
    const double eps_B = 4.0, sig_B = 2.0;

    /* Geometric combination */
    const double eps_AB    = sqrt(eps_A) * sqrt(eps_B); /* = 2.0 */
    const double sigma3_A  = sig_A * sig_A * sig_A;
    const double sigma3_B  = sig_B * sig_B * sig_B;
    const double sigma6_AB = sigma3_A * sigma3_B; /* = 1.0 * 8.0 = 8.0 */

    double fx_geom, fy_geom, fz_geom;
    lj_force(eps_AB, sigma6_AB, cutforcesq, r, 0.0, 0.0, &fx_geom, &fy_geom, &fz_geom);

    /* Independently recompute sigma6 via sigma^2 = cbrt(sigma^6). */
    const double sig2_AB         = cbrt(sigma6_AB);
    const double sigma6_explicit = sig2_AB * sig2_AB * sig2_AB;
    double fx_exp, fy_exp, fz_exp;
    lj_force(eps_AB, sigma6_explicit, cutforcesq, r, 0.0, 0.0, &fx_exp, &fy_exp, &fz_exp);

    ASSERT_NEAR(fx_geom, fx_exp, 1e-10, "geometric force x matches explicit");
    ASSERT_NEAR(fy_geom, fy_exp, 1e-10, "geometric force y matches explicit");
    ASSERT_NEAR(fz_geom, fz_exp, 1e-10, "geometric force z matches explicit");

    /* Single-type equivalence: with uniform types, geometric == single */
    const double eps1 = 1.0, sig1 = 1.0;
    const double sigma6_single = sig1 * sig1 * sig1 * sig1 * sig1 * sig1;
    const double sigma6_geom1  = (sig1 * sig1 * sig1) * (sig1 * sig1 * sig1);
    const double eps_geom1     = sqrt(eps1) * sqrt(eps1);

    double fx_s, fy_s, fz_s, fx_g, fy_g, fz_g;
    lj_force(eps1, sigma6_single, cutforcesq, r, 0.0, 0.0, &fx_s, &fy_s, &fz_s);
    lj_force(eps_geom1, sigma6_geom1, cutforcesq, r, 0.0, 0.0, &fx_g, &fy_g, &fz_g);

    ASSERT_NEAR(fx_g, fx_s, 1e-12, "single-type: geometric == single (x)");
    ASSERT_NEAR(fy_g, fy_s, 1e-12, "single-type: geometric == single (y)");
    ASSERT_NEAR(fz_g, fz_s, 1e-12, "single-type: geometric == single (z)");

    return 0;
}

/* ----------------------------------------------------------------------------
 * Tabulated/spline-interpolated LJ forces (src/common/ljtable.{h,c}).
 * These build the real coefficient table via initLJTable and check that the
 * cubic-spline lookup reproduces the analytic LJ force-over-distance.
 * -------------------------------------------------------------------------- */

/* Analytic LJ force-over-distance (fpair), matching the kernels:
 * force = 48*sr6*(sr6-0.5)*sr2*epsilon, with sr2=1/r^2, sr6=sigma6/r^6. */
static double analytic_fpair(double r, double epsilon, double sigma6)
{
    double rsq = r * r;
    double sr2 = 1.0 / rsq;
    double sr6 = sr2 * sr2 * sr2 * sigma6;
    return 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;
}

/* Replicate the kernel table lookup. Works for both index modes: the kernels
 * use LJ_TABLE_COORD(rsq), i.e. sqrt(rsq)=r for r-indexed and rsq for rsq. */
static double table_fpair(double r, double epsilon, double sigma6)
{
    double coord;
#ifdef LJ_TABLE_RSQ
    coord = r * r;
#else
    coord = r;
#endif
    double u = coord * (double)ljtable.inv_h;
    int m    = (int)u;
    int nm1  = ljtable.n - 1;
    if (m > nm1) m = nm1;
    double eps         = u - (double)m;
    const MD_FLOAT* cc = &ljtable.coeff[m * LJ_TABLE_STRIDE];
    double hrep        = cc[0] + eps * (cc[1] + eps * (cc[2] + eps * cc[3]));
    double gdisp       = cc[4] + eps * (cc[5] + eps * (cc[6] + eps * cc[7]));
    return epsilon * sigma6 * (sigma6 * hrep + gdisp);
}

static void build_test_table(double cut, int points)
{
    Parameter param;
    memset(&param, 0, sizeof(param));
    param.cutforce        = (MD_FLOAT)cut;
    param.lj_table_points = points;
    initLJTable(&param);
}

/* At knot positions the spline value equals the tabulated function exactly
 * (A0 coefficient), so the table must reproduce the analytic force to rounding. */
static int test_lj_table_knot_exactness(void)
{
    const double cut = 2.5;
    build_test_table(cut, 1000);

    const double epsilon = 1.0;
    const double sigma6  = 1.0;

    for (int m = 400; m <= 960; m += 40) {
#ifdef LJ_TABLE_RSQ
        double r = sqrt(m * (double)ljtable.h);
#else
        double r = m * (double)ljtable.h;
#endif
        if (r < 0.9) continue; /* skip the steep, physically-unreached inner region */
        double ftab = table_fpair(r, epsilon, sigma6);
        double fana = analytic_fpair(r, epsilon, sigma6);
        double rel  = fabs(ftab - fana) / (fabs(fana) + 1e-300);
        ASSERT_TRUE(rel < 1e-4, "tabulated LJ matches analytic at knot");
    }
    return 0;
}

/* Between knots the cubic interpolation must stay close to analytic. */
static int test_lj_table_interp_accuracy(void)
{
    const double cut = 2.5;
    build_test_table(cut, 1000);

    const double epsilon = 1.0;
    const double sigma6  = 1.0;

    double max_rel = 0.0;
    for (double r = 0.95; r < cut - 1e-6; r += 0.013) {
        double ftab = table_fpair(r, epsilon, sigma6);
        double fana = analytic_fpair(r, epsilon, sigma6);
        double rel  = fabs(ftab - fana) / (fabs(fana) + 1e-12);
        if (rel > max_rel) max_rel = rel;
    }
    ASSERT_TRUE(max_rel < 1e-3, "tabulated LJ interpolation error within 1e-3");
    return 0;
}

/* The decomposition (eps*sigma^12)*Hrep + (eps*sigma^6)*Gdisp must reproduce
 * the analytic force for non-trivial per-pair epsilon/sigma too. */
static int test_lj_table_prefactor_decomposition(void)
{
    const double cut = 3.0;
    build_test_table(cut, 2000);

    /* Type pair with eps=2.0, sigma=1.2 -> sigma6 = 1.2^6. */
    const double epsilon = 2.0;
    const double sigma   = 1.2;
    const double sigma6  = pow(sigma, 6.0);

    for (double r = 1.0; r < cut - 1e-6; r += 0.05) {
        double ftab = table_fpair(r, epsilon, sigma6);
        double fana = analytic_fpair(r, epsilon, sigma6);
        double rel  = fabs(ftab - fana) / (fabs(fana) + 1e-12);
        ASSERT_TRUE(rel < 2e-3, "tabulated LJ prefactor decomposition matches analytic");
    }
    return 0;
}

/* Near the LJ minimum r = 2^(1/6) the force-over-distance crosses zero. */
static int test_lj_table_force_at_minimum(void)
{
    const double cut = 2.5;
    build_test_table(cut, 1000);

    const double r    = pow(2.0, 1.0 / 6.0);
    double ftab       = table_fpair(r, 1.0, 1.0);
    ASSERT_NEAR(ftab, 0.0, 1e-3, "tabulated LJ force near zero at minimum");
    return 0;
}

int run_force_tests(void)
{
    int rc = 0;

    tr_log("  force: LJ zero at minimum");
    rc = test_lj_zero_force_at_minimum();
    if (rc) return rc;

    tr_log("  force: Newton's third law");
    rc = test_lj_newtons_third_law();
    if (rc) return rc;

    tr_log("  force: cutoff gating");
    rc = test_lj_cutoff_gating();
    if (rc) return rc;

    tr_log("  force: geometric combination formula");
    rc = test_lj_geom_combination_formula();
    if (rc) return rc;

    tr_log("  force: tabulated LJ knot exactness");
    rc = test_lj_table_knot_exactness();
    if (rc) return rc;

    tr_log("  force: tabulated LJ interpolation accuracy");
    rc = test_lj_table_interp_accuracy();
    if (rc) return rc;

    tr_log("  force: tabulated LJ prefactor decomposition");
    rc = test_lj_table_prefactor_decomposition();
    if (rc) return rc;

    tr_log("  force: tabulated LJ force at minimum");
    rc = test_lj_table_force_at_minimum();
    if (rc) return rc;

    return 0;
}
