#!/bin/bash
# Test script for LJ combination rules
#
# Test matrix:
#   verletlist: single/geometric/none x AVX512/AVX2/NONE (scalar)
#               - SIMD+none: expected compile error (SIMD kernel doesn't support none)
#   clusterpair: single/geometric/none x AVX512/AVX2
#               - Scalar not supported (clusterpair requires SIMD)
#               - SIMD+none: works (uses type-pair indexing)

set -e
cd "$(dirname "$0")/.."

PASS=0
FAIL=0

test_build_run() {
    local opt_scheme="$1"
    local simd="$2"
    local lj_comb_rule="$3"
    local expect_fail="$4"
    local name="${opt_scheme}-${simd}-${lj_comb_rule}"

    printf "%-35s " "$name"
    rm -rf ./build ./MDBench-*

    local build_output
    build_output=$(make -j4 TOOLCHAIN=GCC ISA=X86 SIMD=$simd DATA_TYPE=DP OPT_SCHEME=$opt_scheme LJ_COMB_RULE=$lj_comb_rule 2>&1)
    local build_result=$?

    if [ $build_result -eq 0 ]; then
        if [ "$expect_fail" = "yes" ]; then
            echo "FAIL (should not compile)"
            ((FAIL++))
            return 1
        fi

        local binary=$(ls MDBench-* 2>/dev/null | head -1)
        local output=$(./"$binary" -n 3 2>&1)
        local rule=$(echo "$output" | grep "LJ combination rule:" | awk '{print $4}')

        if [ "$rule" = "$lj_comb_rule" ]; then
            echo "PASS"
            ((PASS++))
        else
            echo "FAIL (expected '$lj_comb_rule', got '$rule')"
            ((FAIL++))
        fi
    else
        if [ "$expect_fail" = "yes" ]; then
            echo "PASS (expected compile error)"
            ((PASS++))
        else
            echo "FAIL (build failed)"
            echo "$build_output" | grep -E "error:" | head -3
            ((FAIL++))
        fi
    fi
}

echo "============================================"
echo "LJ Combination Rule Tests"
echo "============================================"
echo ""

echo "--- Verletlist SIMD (AVX512) ---"
test_build_run verletlist AVX512 single no
test_build_run verletlist AVX512 geometric no
test_build_run verletlist AVX512 none yes

echo ""
echo "--- Verletlist SIMD (AVX2) ---"
test_build_run verletlist AVX2 single no
test_build_run verletlist AVX2 geometric no
test_build_run verletlist AVX2 none yes

echo ""
echo "--- Verletlist Scalar ---"
test_build_run verletlist NONE single no
test_build_run verletlist NONE geometric no
test_build_run verletlist NONE none no

echo ""
echo "--- Clusterpair SIMD (AVX512) ---"
test_build_run clusterpair AVX512 single no
test_build_run clusterpair AVX512 geometric no
test_build_run clusterpair AVX512 none no

echo ""
echo "--- Clusterpair SIMD (AVX2) ---"
test_build_run clusterpair AVX2 single no
test_build_run clusterpair AVX2 geometric no
test_build_run clusterpair AVX2 none no

rm -rf ./build ./MDBench-*

echo ""
echo "============================================"
echo "Summary: $PASS passed, $FAIL failed"
echo "============================================"

[ $FAIL -eq 0 ]
