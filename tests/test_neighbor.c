#include "test_runner.h"

#include <atom.h>
#include <force.h>
#include <neighbor.h>
#include <parameter.h>
#include <pbc.h>
#include <util.h>

/* Local copy of the bounding box distance used in neighbor.c. */
static MD_FLOAT getBoundingBoxDistanceSq_test(Atom* atom, int ci, int cj)
{
    MD_FLOAT dl  = atom->iclusters[ci].bbminx - atom->jclusters[cj].bbmaxx;
    MD_FLOAT dh  = atom->jclusters[cj].bbminx - atom->iclusters[ci].bbmaxx;
    MD_FLOAT dm  = MAX(dl, dh);
    MD_FLOAT dm0 = MAX(dm, 0.0);
    MD_FLOAT d2  = dm0 * dm0;

    dl  = atom->iclusters[ci].bbminy - atom->jclusters[cj].bbmaxy;
    dh  = atom->jclusters[cj].bbminy - atom->iclusters[ci].bbmaxy;
    dm  = MAX(dl, dh);
    dm0 = MAX(dm, 0.0);
    d2 += dm0 * dm0;

    dl  = atom->iclusters[ci].bbminz - atom->jclusters[cj].bbmaxz;
    dh  = atom->jclusters[cj].bbminz - atom->iclusters[ci].bbmaxz;
    dm  = MAX(dl, dh);
    dm0 = MAX(dm, 0.0);
    d2 += dm0 * dm0;
    return d2;
}

/* Build a small, deterministic system and its neighbor lists. */
static void build_small_system(Parameter* param, Atom* atom, Neighbor* neighbor)
{
    initParameter(param);
    /* Keep the default LJ solid but reduce system size. */
    param->nx         = 4;
    param->ny         = 4;
    param->nz         = 4;
    param->ntimes     = 0;
    param->half_neigh = 0;

    initAtom(atom);
    /* Neighbor setup does not require full force initialization. */
    initNeighbor(neighbor, param);

    /* Atom positions from FCC lattice. */
    createAtom(atom, param);

    /* Build neighbor infrastructure following clusterpair/main.c::setup. */
    setupNeighbor(param, atom);
    buildClusters(atom);
    defineJClusters(param, atom);
    setupPbc(atom, param);
    binJClusters(param, atom);
    buildNeighbor(atom, neighbor);
}

static int test_neighbor_vs_bruteforce_bounding_boxes(void)
{
    Parameter param;
    Atom atom;
    Neighbor neighbor;

    build_small_system(&param, &atom, &neighbor);

    const MD_FLOAT cutneigh   = param.cutneigh;
    const MD_FLOAT cutneighsq = cutneigh * cutneigh;
    const int nci       = atom.Nclusters_local;
    const int ncj_total = atom.ncj + atom.Nclusters_ghost;
    const int nbM       = nci;

    /* If no clusters or neighbor arrays are present (e.g., degenerate setup),
       skip this check without failing the whole test suite. */
    if (nci == 0 || neighbor.numneigh == NULL || neighbor.neighbors == NULL) {
        return 0;
    }

    /* For each i-cluster, ensure all bbox-close j-clusters are in the list,
       and no far j-clusters are erroneously present. */
    for (int ci = 0; ci < nci; ++ci) {
        int numneigh = neighbor.numneigh[ci];
        ASSERT_TRUE(numneigh >= 0, "numneigh non-negative");

        /* Mark neighbors present in the list. */
        int* present = (int*)calloc(ncj_total, sizeof(int));
        ASSERT_TRUE(present != NULL, "alloc present[]");

        for (int k = 0; k < numneigh; ++k) {
            int cj = neighs(neighbor.neighbors, ci, k, nbM, &neighbor);
            ASSERT_TRUE(cj >= 0 && cj < ncj_total, "neighbor index in range");
            present[cj] = 1;
        }

        /* Check completeness: every bbox-close cj must appear. */
        int missing = 0;
        for (int cj = 0; cj < ncj_total; ++cj) {
            MD_FLOAT d2 = getBoundingBoxDistanceSq_test(&atom, ci, cj);
            if (d2 < cutneighsq) {
                if (!present[cj]) {
                    missing = 1;
                    break;
                }
            }
        }

        /* Check soundness: any listed neighbor should be reasonably close. */
        int spurious          = 0;
        const MD_FLOAT margin = (MD_FLOAT)1.01; /* tiny numerical slack */
        for (int cj = 0; cj < ncj_total && !spurious; ++cj) {
            if (present[cj]) {
                MD_FLOAT d2 = getBoundingBoxDistanceSq_test(&atom, ci, cj);
                if (d2 > margin * cutneighsq) {
                    spurious = 1;
                }
            }
        }

        free(present);

        ASSERT_TRUE(!missing, "all bbox-close clusters appear in neighbor list");
        ASSERT_TRUE(!spurious, "no far clusters appear in neighbor list");
    }

    return 0;
}

static int test_store_and_check_displacement(void)
{
    Parameter param;
    Atom atom;
    Neighbor neighbor;

    build_small_system(&param, &atom, &neighbor);

    if (atom.Nclusters_local == 0 || atom.iclusters[0].natoms == 0) {
        return 0;
    }

    // After storing references, no displacement has occurred yet.
    storeReferencePositions(&atom);
    ASSERT_TRUE(!needsReneigh(&atom, &param),
        "needsReneigh false immediately after storeReferencePositions");

    // Move one atom by less than skin/2 — threshold = (skin/2)^2.
    // skin defaults to 0.3, so threshold = 0.0225; move by 0.14 → sq=0.0196, no trigger.
    int base                        = CI_VECTOR3_BASE_INDEX(0);
    atom.cl_x[base + CL_X_INDEX_3D(0)] += (MD_FLOAT)0.14;
    ASSERT_TRUE(!needsReneigh(&atom, &param),
        "needsReneigh false when displacement < skin/2");

    // Move past skin/2: total displacement = 0.16 → sq = 0.0256 > 0.0225.
    atom.cl_x[base + CL_X_INDEX_3D(0)] += (MD_FLOAT)0.02;
    ASSERT_TRUE(needsReneigh(&atom, &param),
        "needsReneigh true when displacement > skin/2");

    // After storing the new reference the flag clears.
    storeReferencePositions(&atom);
    ASSERT_TRUE(!needsReneigh(&atom, &param),
        "needsReneigh false after re-storing reference positions");

    return 0;
}

static int test_displacement_threshold_exact(void)
{
    Parameter param;
    Atom atom;
    Neighbor neighbor;

    build_small_system(&param, &atom, &neighbor);

    if (atom.Nclusters_local == 0 || atom.iclusters[0].natoms == 0) {
        return 0;
    }

    // Use a known skin and verify threshold boundary precisely.
    param.skin = (MD_FLOAT)1.0;
    // threshold = (0.5)^2 = 0.25

    storeReferencePositions(&atom);

    // Displace by exactly skin/2 in each of two atoms (different clusters if available).
    int base0 = CI_VECTOR3_BASE_INDEX(0);
    atom.cl_x[base0 + CL_X_INDEX_3D(0)] += (MD_FLOAT)0.5;
    // Displacement sq = 0.25, which equals threshold — should NOT trigger (strict >).
    ASSERT_TRUE(!needsReneigh(&atom, &param),
        "needsReneigh false when displacement equals skin/2 exactly");

    // One epsilon over threshold.
    atom.cl_x[base0 + CL_X_INDEX_3D(0)] += (MD_FLOAT)1e-4;
    ASSERT_TRUE(needsReneigh(&atom, &param),
        "needsReneigh true when displacement infinitesimally exceeds skin/2");

    return 0;
}

int run_neighbor_tests(void)
{
    int rc = 0;

    tr_log("  neighbor: bbox vs list consistency");
    rc = test_neighbor_vs_bruteforce_bounding_boxes();
    if (rc) return rc;

    tr_log("  neighbor: displacement reneigh store and check");
    rc = test_store_and_check_displacement();
    if (rc) return rc;

    tr_log("  neighbor: displacement reneigh threshold boundary");
    rc = test_displacement_threshold_exact();
    if (rc) return rc;

    return 0;
}
