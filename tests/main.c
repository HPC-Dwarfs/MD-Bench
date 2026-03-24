#include "test_runner.h"

int run_parameter_tests(void);
int run_atom_tests(void);
int run_force_tests(void);
int run_neighbor_tests(void);

int main(void)
{
    int rc = 0;

    tr_log("Running MD-Bench unit tests...");

    rc = run_parameter_tests();
    if (rc) {
        tr_log("parameter tests FAILED");
        return rc;
    }
    tr_log("parameter tests OK");

    rc = run_atom_tests();
    if (rc) {
        tr_log("atom tests FAILED");
        return rc;
    }
    tr_log("atom tests OK");

    rc = run_force_tests();
    if (rc) {
        tr_log("force tests FAILED");
        return rc;
    }
    tr_log("force tests OK");

    rc = run_neighbor_tests();
    if (rc) {
        tr_log("neighbor tests FAILED");
        return rc;
    }
    tr_log("neighbor tests OK");

    tr_log("All tests OK");
    return 0;
}
