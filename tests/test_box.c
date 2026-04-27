#include "test_runner.h"

#include <box.h>
#include <parameter.h>

/*
 * Unit tests for the box overlap functions used by MPI domain decomposition.
 * These test overlapBox() and overlapFullBox() from box.c without requiring
 * an actual MPI runtime.
 */

/* Test 1: non-overlapping boxes.
 * Two boxes that are far apart should return -100 (no overlap).
 */
static int test_overlapBox_no_overlap(void)
{
    Box mybox = { .id = 0, .lo = { 0.0, 0.0, 0.0 }, .hi = { 5.0, 5.0, 5.0 } };
    Box other = { .id = 1, .lo = { 10.0, 10.0, 10.0 }, .hi = { 15.0, 15.0, 15.0 } };
    Box cut;
    MD_FLOAT xprd     = 30.0; /* large enough that PBC images don't overlap either */
    MD_FLOAT cutneigh = 1.0;

    int pbc = overlapBox(0, 0, &mybox, &other, &cut, xprd, cutneigh);
    ASSERT_INT_EQ(pbc, -100, "non-overlapping boxes return -100");

    return 0;
}

/* Test 2: adjacent boxes with cutoff overlap.
 * Box [0,5] and [5,10] along x with cutneigh=1 should overlap.
 * The backward direction (dir=1) extends other->lo by cutneigh into mybox.
 */
static int test_overlapBox_adjacent_overlap(void)
{
    Box mybox = { .id = 0,
        .lo                   = { 0.0, 0.0, 0.0 },
        .hi                   = { 5.0, 5.0, 5.0 } };
    Box other = { .id = 1,
        .lo                   = { 5.0, 0.0, 0.0 },
        .hi                   = { 10.0, 5.0, 5.0 } };
    Box cut;
    MD_FLOAT xprd     = 20.0;
    MD_FLOAT cutneigh = 1.0;

    /* Backward direction (dir=1): extends other->lo by -cutneigh, creating
     * an overlap region [4,5] x [0,5] x [0,5] inside mybox. */
    int pbc = overlapBox(0, 1, &mybox, &other, &cut, xprd, cutneigh);
    ASSERT_TRUE(pbc != -100, "adjacent boxes should overlap with cutoff");
    ASSERT_INT_EQ(pbc, 0, "adjacent non-periodic overlap => pbc=0");

    /* Also verify the cut box covers the expected overlap region. */
    ASSERT_NEAR(cut.lo[0], 4.0, 1e-12, "cut lo[0] = mybox->hi - cutneigh");
    ASSERT_NEAR(cut.hi[0], 5.0, 1e-12, "cut hi[0] = mybox->hi");

    return 0;
}

/* Test 3: periodic overlap.
 * Box at the left edge of the domain should overlap with box at the right edge
 * via PBC wrapping.
 */
static int test_overlapBox_periodic(void)
{
    MD_FLOAT xprd     = 10.0;
    MD_FLOAT cutneigh = 1.0;

    /* mybox is at the left end, other is at the right end.
     * In the forward direction (dir=0), the periodic image of other
     * is shifted left by xprd: [8-10, 10-10] = [-2, 0], which
     * is within cutneigh of mybox. */
    Box mybox = { .id = 0,
        .lo                   = { 0.0, 0.0, 0.0 },
        .hi                   = { 2.0, 5.0, 5.0 } };
    Box other = { .id = 1,
        .lo                   = { 8.0, 0.0, 0.0 },
        .hi                   = { 10.0, 5.0, 5.0 } };
    Box cut;

    int pbc = overlapBox(0, 0, &mybox, &other, &cut, xprd, cutneigh);
    ASSERT_TRUE(pbc != -100, "periodic boxes should overlap");
    ASSERT_TRUE(pbc == 1 || pbc == -1, "periodic overlap should have non-zero pbc");

    return 0;
}

/* Test 4: overlapFullBox checks all 27 periodic images.
 * Nearby boxes should overlap, far-away ones should not.
 */
static int test_overlapFullBox_nearby(void)
{
    Parameter param;
    initParameter(&param);
    param.xprd = 10.0;
    param.yprd = 10.0;
    param.zprd = 10.0;

    MD_FLOAT cutneigh[3] = { 1.0, 1.0, 1.0 };

    /* Two adjacent boxes along x */
    Box mybox = { .id = 0,
        .lo                   = { 0.0, 0.0, 0.0 },
        .hi                   = { 5.0, 5.0, 5.0 } };
    Box nearby = { .id = 1,
        .lo                    = { 4.5, 0.0, 0.0 },
        .hi                    = { 9.5, 5.0, 5.0 } };

    int result = overlapFullBox(&param, cutneigh, &mybox, &nearby);
    ASSERT_INT_EQ(result, 1, "nearby boxes overlap");

    /* Box far away in all dimensions */
    Box faraway = { .id = 2,
        .lo                     = { 7.0, 7.0, 7.0 },
        .hi                     = { 8.0, 8.0, 8.0 } };

    result = overlapFullBox(&param, cutneigh, &mybox, &faraway);
    ASSERT_INT_EQ(result, 0, "distant boxes do not overlap");

    return 0;
}

/* Test 5: overlapFullBox detects periodic neighbor.
 * A box at [9,10] and a box at [0,1] should overlap periodically.
 */
static int test_overlapFullBox_periodic(void)
{
    Parameter param;
    initParameter(&param);
    param.xprd = 10.0;
    param.yprd = 10.0;
    param.zprd = 10.0;

    MD_FLOAT cutneigh[3] = { 1.5, 1.5, 1.5 };

    Box mybox = { .id = 0,
        .lo                   = { 9.0, 0.0, 0.0 },
        .hi                   = { 10.0, 5.0, 5.0 } };
    Box other = { .id = 1,
        .lo                   = { 0.0, 0.0, 0.0 },
        .hi                   = { 1.0, 5.0, 5.0 } };

    int result = overlapFullBox(&param, cutneigh, &mybox, &other);
    ASSERT_INT_EQ(result, 1, "periodic neighbor detected by overlapFullBox");

    return 0;
}

/* Test 6: self-overlap.
 * A box should not overlap with itself via the non-periodic path of overlapBox
 * (the same == 1 check skips the non-periodic path).
 */
static int test_overlapBox_self(void)
{
    Box mybox = { .id = 0,
        .lo                   = { 0.0, 0.0, 0.0 },
        .hi                   = { 5.0, 5.0, 5.0 } };
    Box cut;
    MD_FLOAT xprd     = 20.0;
    MD_FLOAT cutneigh = 1.0;

    int pbc = overlapBox(0, 0, &mybox, &mybox, &cut, xprd, cutneigh);
    /* With same box (id==id), the non-periodic path is skipped.
     * PBC path: lo_shifted = lo - xprd = -20, hi_shifted = hi - xprd + cutneigh = -14
     * This doesn't overlap [0,5], so pbc should be -100. */
    ASSERT_INT_EQ(pbc, -100, "self-overlap returns -100 (no PBC duplication needed)");

    return 0;
}

int run_box_tests(void)
{
    int rc = 0;

    tr_log("  box: non-overlapping boxes");
    rc = test_overlapBox_no_overlap();
    if (rc) return rc;

    tr_log("  box: adjacent boxes with cutoff overlap");
    rc = test_overlapBox_adjacent_overlap();
    if (rc) return rc;

    tr_log("  box: periodic overlap");
    rc = test_overlapBox_periodic();
    if (rc) return rc;

    tr_log("  box: overlapFullBox nearby");
    rc = test_overlapFullBox_nearby();
    if (rc) return rc;

    tr_log("  box: overlapFullBox periodic");
    rc = test_overlapFullBox_periodic();
    if (rc) return rc;

    tr_log("  box: self-overlap");
    rc = test_overlapBox_self();
    if (rc) return rc;

    return 0;
}
