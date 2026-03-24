#include "test_runner.h"

#include <parameter.h>
#include <force.h>

/* Simple helper to write a small temporary parameter file. */
static const char* write_temp_param_file(void)
{
    const char* fname = "test_params.conf";
    FILE* fp          = fopen(fname, "w");
    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    /* Only a subset of parameters needed to verify parsing logic. */
    fprintf(fp, "epsilon 0.8\n");
    fprintf(fp, "sigma 1.1\n");
    fprintf(fp, "nx 10\n");
    fprintf(fp, "ny 20\n");
    fprintf(fp, "nz 30\n");
    fprintf(fp, "dt 0.001\n");
    fprintf(fp, "cutforce 2.0\n");
    fprintf(fp, "skin 0.2\n");
    fprintf(fp, "reneigh_every 5\n");
    fprintf(fp, "balance_every 3\n");

    fclose(fp);
    return fname;
}

static int test_initParameter_defaults(void)
{
    Parameter p;
    initParameter(&p);

    ASSERT_TRUE(p.input_file == NULL, "input_file should default to NULL");
    ASSERT_TRUE(p.vtk_file == NULL, "vtk_file should default to NULL");
    ASSERT_TRUE(p.eam_file == NULL, "eam_file should default to NULL");
    ASSERT_TRUE(p.write_atom_file == NULL, "write_atom_file should default to NULL");

    ASSERT_INT_EQ(p.force_field, FF_LJ, "default force_field is LJ");
    ASSERT_NEAR(p.epsilon, 1.0, 1e-12, "default epsilon");
    ASSERT_NEAR(p.sigma, 1.0, 1e-12, "default sigma");
    ASSERT_NEAR(p.rho, 0.8442, 1e-12, "default rho");

    ASSERT_INT_EQ(p.ntimes, 200, "default ntimes");
    ASSERT_NEAR(p.dt, 0.005, 1e-12, "default dt");
    ASSERT_INT_EQ(p.nx, 32, "default nx");
    ASSERT_INT_EQ(p.ny, 32, "default ny");
    ASSERT_INT_EQ(p.nz, 32, "default nz");

    ASSERT_NEAR(p.cutforce, 2.5, 1e-12, "default cutforce");
    ASSERT_NEAR(p.skin, 0.3, 1e-12, "default skin");
    ASSERT_NEAR(p.cutneigh, p.cutforce + p.skin, 1e-12, "cutneigh = cutforce + skin");

    ASSERT_NEAR(p.dtforce, 0.5 * p.dt, 1e-12, "dtforce = 0.5 * dt");

    ASSERT_INT_EQ(p.reneigh_every, 20, "default reneigh_every");
    ASSERT_INT_EQ(p.balance_every, p.reneigh_every, "default balance_every = reneigh_every");

    return 0;
}

static int test_readParameter_overrides_and_derived(void)
{
    Parameter p;
    initParameter(&p);

    const char* fname = write_temp_param_file();
    readParameter(&p, fname);

    /* Direct values from file */
    ASSERT_NEAR(p.epsilon, 0.8, 1e-12, "epsilon from file");
    ASSERT_NEAR(p.sigma, 1.1, 1e-12, "sigma from file");
    ASSERT_INT_EQ(p.nx, 10, "nx from file");
    ASSERT_INT_EQ(p.ny, 20, "ny from file");
    ASSERT_INT_EQ(p.nz, 30, "nz from file");
    ASSERT_NEAR(p.dt, 0.001, 1e-12, "dt from file");
    ASSERT_NEAR(p.cutforce, 2.0, 1e-12, "cutforce from file");
    ASSERT_NEAR(p.skin, 0.2, 1e-12, "skin from file");

    /* Derived values updated in readParameter */
    ASSERT_NEAR(p.dtforce, 0.5 * p.dt, 1e-12, "dtforce updated");

    /* sigma6 recomputed as sigma^6 (within tight tolerance) */
    double s2      = p.sigma * p.sigma;
    double sigma6e = s2 * s2 * s2;
    ASSERT_NEAR(p.sigma6, sigma6e, 1e-12, "sigma6 recomputed");

    /* balance_every scaled by reneigh_every as in parameter.c */
    ASSERT_INT_EQ(p.balance_every, 3 * p.reneigh_every, "balance_every scaled");

    return 0;
}

static int test_computePerTypeLJParameters(void)
{
    Parameter p;
    initParameter(&p);

    MD_FLOAT eps_arr[2] = { 1.0, 4.0 };
    MD_FLOAT sig_arr[2] = { 1.0, 2.0 };
    p.ntypes           = 2;
    p.epsilon_per_type = eps_arr;
    p.sigma_per_type   = sig_arr;

    MD_FLOAT sqrt_eps[2], sigma3[2];
    computePerTypeLJParameters(2, &p, sqrt_eps, sigma3);

    ASSERT_NEAR(sqrt_eps[0], 1.0, 1e-12, "sqrt_epsilon[0] == sqrt(1.0)");
    ASSERT_NEAR(sqrt_eps[1], 2.0, 1e-12, "sqrt_epsilon[1] == sqrt(4.0)");
    ASSERT_NEAR(sigma3[0], 1.0, 1e-12, "sigma3[0] == 1.0^3");
    ASSERT_NEAR(sigma3[1], 8.0, 1e-12, "sigma3[1] == 2.0^3");

    /* Fallback: no per-type arrays → use global epsilon/sigma */
    p.epsilon_per_type = NULL;
    p.sigma_per_type   = NULL;
    p.epsilon          = 9.0;
    p.sigma            = 3.0;
    computePerTypeLJParameters(2, &p, sqrt_eps, sigma3);

    ASSERT_NEAR(sqrt_eps[0], 3.0, 1e-12, "fallback sqrt_epsilon[0] == sqrt(9.0)");
    ASSERT_NEAR(sqrt_eps[1], 3.0, 1e-12, "fallback sqrt_epsilon[1] == sqrt(9.0)");
    ASSERT_NEAR(sigma3[0], 27.0, 1e-12, "fallback sigma3[0] == 3.0^3");
    ASSERT_NEAR(sigma3[1], 27.0, 1e-12, "fallback sigma3[1] == 3.0^3");

    return 0;
}

static int test_computeTypePairLJParameters(void)
{
    MD_FLOAT sqrt_eps[2] = { 1.0, 2.0 };
    MD_FLOAT sigma3[2]   = { 1.0, 8.0 };
    MD_FLOAT epsilon[4], sigma6[4];

    computeTypePairLJParameters(2, sqrt_eps, sigma3, epsilon, sigma6);

    /* epsilon[i*2+j] = sqrt_eps[i] * sqrt_eps[j] */
    ASSERT_NEAR(epsilon[0 * 2 + 0], 1.0, 1e-12, "epsilon[0,0] = 1*1");
    ASSERT_NEAR(epsilon[0 * 2 + 1], 2.0, 1e-12, "epsilon[0,1] = 1*2");
    ASSERT_NEAR(epsilon[1 * 2 + 0], 2.0, 1e-12, "epsilon[1,0] = 2*1");
    ASSERT_NEAR(epsilon[1 * 2 + 1], 4.0, 1e-12, "epsilon[1,1] = 2*2");

    /* sigma6[i*2+j] = sigma3[i] * sigma3[j] */
    ASSERT_NEAR(sigma6[0 * 2 + 0], 1.0, 1e-12, "sigma6[0,0] = 1*1");
    ASSERT_NEAR(sigma6[0 * 2 + 1], 8.0, 1e-12, "sigma6[0,1] = 1*8");
    ASSERT_NEAR(sigma6[1 * 2 + 0], 8.0, 1e-12, "sigma6[1,0] = 8*1");
    ASSERT_NEAR(sigma6[1 * 2 + 1], 64.0, 1e-12, "sigma6[1,1] = 8*8");

    return 0;
}

static const char* write_temp_multitype_param_file(void)
{
    const char* fname = "test_multitype_params.conf";
    FILE* fp          = fopen(fname, "w");
    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "epsilon 9.0\n");
    fprintf(fp, "sigma 3.0\n");
    fprintf(fp, "ntypes 2\n");
    fprintf(fp, "epsilon_type_0 1.0\n");
    fprintf(fp, "epsilon_type_1 4.0\n");
    fprintf(fp, "sigma_type_0 1.0\n");
    fprintf(fp, "sigma_type_1 2.0\n");

    fclose(fp);
    return fname;
}

static int test_readParameter_multitype_lj(void)
{
    Parameter p;
    initParameter(&p);

    const char* fname = write_temp_multitype_param_file();
    readParameter(&p, fname);

    ASSERT_INT_EQ(p.ntypes, 2, "ntypes == 2");
    ASSERT_TRUE(p.epsilon_per_type != NULL, "epsilon_per_type allocated");
    ASSERT_TRUE(p.sigma_per_type != NULL, "sigma_per_type allocated");

    ASSERT_NEAR(p.epsilon_per_type[0], 1.0, 1e-12, "epsilon_type_0");
    ASSERT_NEAR(p.epsilon_per_type[1], 4.0, 1e-12, "epsilon_type_1");
    ASSERT_NEAR(p.sigma_per_type[0], 1.0, 1e-12, "sigma_type_0");
    ASSERT_NEAR(p.sigma_per_type[1], 2.0, 1e-12, "sigma_type_1");

    /* Global epsilon/sigma remain as set in the file */
    ASSERT_NEAR(p.epsilon, 9.0, 1e-12, "global epsilon unchanged");

    return 0;
}

int run_parameter_tests(void)
{
    int rc = 0;

    tr_log("  parameter: initParameter defaults");
    rc = test_initParameter_defaults();
    if (rc)
        return rc;

    tr_log("  parameter: readParameter overrides and derived values");
    rc = test_readParameter_overrides_and_derived();
    if (rc)
        return rc;

    tr_log("  parameter: computePerTypeLJParameters math and fallback");
    rc = test_computePerTypeLJParameters();
    if (rc)
        return rc;

    tr_log("  parameter: computeTypePairLJParameters NxN matrix");
    rc = test_computeTypePairLJParameters();
    if (rc)
        return rc;

    tr_log("  parameter: readParameter multi-type LJ parsing");
    rc = test_readParameter_multitype_lj();
    if (rc)
        return rc;

    return 0;
}

