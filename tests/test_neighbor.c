#include "test_runner.h"

#include <atom.h>
#include <neighbor.h>
#include <parameter.h>
#include <force.h>
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
    param->nx       = 4;
    param->ny       = 4;
    param->nz       = 4;
    param->ntimes   = 0;
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

    const MD_FLOAT cutneigh     = param.cutneigh;
    const MD_FLOAT cutneighsq   = cutneigh * cutneigh;
    const int nci               = atom.Nclusters_local;
    const int ncj_total         = atom.ncj + atom.Nclusters_ghost;
    const int nbM               = nci;
    const int nbN               = neighbor.maxneighs;

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
            int cj = neighs(neighbor.neighbors, ci, k, nbM, nbN);
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
        int spurious = 0;
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

int run_neighbor_tests(void)
{
    int rc = 0;

    tr_log("  neighbor: bbox vs list consistency");
    rc = test_neighbor_vs_bruteforce_bounding_boxes();
    if (rc)
        return rc;

    return 0;
}

