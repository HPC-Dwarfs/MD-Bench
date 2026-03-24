#!/bin/bash
# Test script for LJ combination rules with output verification
#
# Tests:
# 1. All configurations build and report correct LJ rule
# 2. Same configuration produces consistent output across runs
# 3. Single rule with 1 type == Geometric/None rules (same physics)

set -e
cd "$(dirname "$0")/.."

PASS=0
FAIL=0

# Reference values for output comparison
declare -A TEMP_REF
declare -A PRES_REF

build_and_run() {
    local opt_scheme="$1"
    local simd="$2"
    local lj_rule="$3"
    local timesteps="${4:-10}"

    rm -rf ./build ./MDBench-*
    if ! make -j4 TOOLCHAIN=GCC ISA=X86 SIMD=$simd DATA_TYPE=DP OPT_SCHEME=$opt_scheme LJ_COMB_RULE=$lj_rule >/dev/null 2>&1; then
        return 1
    fi

    local binary=$(ls MDBench-* 2>/dev/null | head -1)
    ./"$binary" -n $timesteps 2>&1
}

extract_final_temp() {
    grep "^  10 " | awk '{print $2}'
}

extract_final_pres() {
    grep "^  10 " | awk '{print $3}'
}

compare_floats() {
    local a="$1"
    local b="$2"
    local tol="${3:-1e-10}"

    # Use awk for float comparison with relative tolerance
    awk -v a="$a" -v b="$b" -v tol="$tol" 'BEGIN {
        diff = (a - b) / (a + 1e-30);
        if (diff < 0) diff = -diff;
        exit (diff > tol)
    }'
}

test_config() {
    local opt_scheme="$1"
    local simd="$2"
    local lj_rule="$3"
    local expect_fail="${4:-no}"
    local name="${opt_scheme}-${simd}-${lj_rule}"

    printf "%-40s " "$name"

    local output
    output=$(build_and_run "$opt_scheme" "$simd" "$lj_rule" 10 2>&1)
    local build_result=$?

    if [ $build_result -ne 0 ]; then
        if [ "$expect_fail" = "yes" ]; then
            echo "PASS (expected compile error)"
            ((PASS++))
            return 0
        else
            echo "FAIL (build failed)"
            ((FAIL++))
            return 1
        fi
    fi

    # Verify LJ rule is reported correctly
    local reported=$(echo "$output" | grep "LJ combination rule:" | awk '{print $4}')
    if [ "$reported" != "$lj_rule" ]; then
        echo "FAIL (rule: expected '$lj_rule', got '$reported')"
        ((FAIL++))
        return 1
    fi

    # Extract final temperature and pressure
    local temp=$(echo "$output" | extract_final_temp)
    local pres=$(echo "$output" | extract_final_pres)

    if [ -z "$temp" ] || [ -z "$pres" ]; then
        echo "FAIL (no output values)"
        ((FAIL++))
        return 1
    fi

    # Store for reference comparison
    TEMP_REF["$name"]="$temp"
    PRES_REF["$name"]="$pres"

    echo "PASS (T=$temp, P=$pres)"
    ((PASS++))
    return 0
}

test_output_consistency() {
    local name1="$1"
    local name2="$2"
    local desc="$3"

    printf "%-40s " "$desc"

    local t1="${TEMP_REF[$name1]}"
    local t2="${TEMP_REF[$name2]}"
    local p1="${PRES_REF[$name1]}"
    local p2="${PRES_REF[$name2]}"

    if [ -z "$t1" ] || [ -z "$t2" ]; then
        echo "SKIP (missing reference)"
        return 0
    fi

    if compare_floats "$t1" "$t2" 1e-6 && compare_floats "$p1" "$p2" 1e-6; then
        echo "PASS"
        ((PASS++))
    else
        echo "FAIL (T: $t1 vs $t2, P: $p1 vs $p2)"
        ((FAIL++))
    fi
}

# Like test_output_consistency but with a caller-specified tolerance.
# Use for cross-scheme comparisons (VL vs CP) where numerical differences
# are expected to be larger due to different neighbour-list structures.
test_output_consistency_tol() {
    local name1="$1"
    local name2="$2"
    local desc="$3"
    local tol="${4:-1e-4}"

    printf "%-40s " "$desc"

    local t1="${TEMP_REF[$name1]}"
    local t2="${TEMP_REF[$name2]}"
    local p1="${PRES_REF[$name1]}"
    local p2="${PRES_REF[$name2]}"

    if [ -z "$t1" ] || [ -z "$t2" ]; then
        echo "SKIP (missing reference)"
        return 0
    fi

    if compare_floats "$t1" "$t2" "$tol" && compare_floats "$p1" "$p2" "$tol"; then
        echo "PASS"
        ((PASS++))
    else
        echo "FAIL (T: $t1 vs $t2, P: $p1 vs $p2)"
        ((FAIL++))
    fi
}

echo "============================================"
echo "LJ Combination Rule Tests with Output Verification"
echo "============================================"
echo ""

echo "=== Build and Basic Output Tests ==="
echo ""

echo "--- Verletlist AVX512 ---"
test_config verletlist AVX512 single
test_config verletlist AVX512 geometric
test_config verletlist AVX512 none yes

echo ""
echo "--- Verletlist AVX2 ---"
test_config verletlist AVX2 single
test_config verletlist AVX2 geometric
test_config verletlist AVX2 none yes

echo ""
echo "--- Verletlist Scalar ---"
test_config verletlist NONE single
test_config verletlist NONE geometric
test_config verletlist NONE none

echo ""
echo "--- Clusterpair AVX512 ---"
test_config clusterpair AVX512 single
test_config clusterpair AVX512 geometric
test_config clusterpair AVX512 none

echo ""
echo "--- Clusterpair AVX2 ---"
test_config clusterpair AVX2 single
test_config clusterpair AVX2 geometric
test_config clusterpair AVX2 none

echo ""
echo "=== Output Consistency Tests ==="
echo "(With 1 atom type, single/geometric/none should produce identical results)"
echo ""

# Same scheme, different SIMD width should match
test_output_consistency "verletlist-AVX512-single" "verletlist-AVX2-single" \
    "VL single: AVX512 vs AVX2"
test_output_consistency "verletlist-AVX512-geometric" "verletlist-AVX2-geometric" \
    "VL geometric: AVX512 vs AVX2"
test_output_consistency "verletlist-AVX2-single" "verletlist-NONE-single" \
    "VL single: AVX2 vs scalar"
test_output_consistency "verletlist-AVX2-geometric" "verletlist-NONE-geometric" \
    "VL geometric: AVX2 vs scalar"

# With 1 atom type, all rules should match (same physics)
test_output_consistency "verletlist-NONE-single" "verletlist-NONE-geometric" \
    "VL scalar: single vs geometric"
test_output_consistency "verletlist-NONE-geometric" "verletlist-NONE-none" \
    "VL scalar: geometric vs none"

# Clusterpair consistency
test_output_consistency "clusterpair-AVX512-single" "clusterpair-AVX2-single" \
    "CP single: AVX512 vs AVX2"
test_output_consistency "clusterpair-AVX512-geometric" "clusterpair-AVX2-geometric" \
    "CP geometric: AVX512 vs AVX2"
test_output_consistency "clusterpair-AVX2-single" "clusterpair-AVX2-geometric" \
    "CP AVX2: single vs geometric"

# Cross-scheme: VL scalar vs CP AVX512 for each rule.
# Uses a looser tolerance (1e-4) because VL and CP build different neighbour lists.
test_output_consistency_tol "verletlist-NONE-geometric" "clusterpair-AVX512-geometric" \
    "VL scalar vs CP AVX512: geometric" 1e-4
test_output_consistency_tol "verletlist-NONE-single" "clusterpair-AVX512-single" \
    "VL scalar vs CP AVX512: single" 1e-4

rm -rf ./build ./MDBench-*

echo ""
echo "============================================"
echo "Summary: $PASS passed, $FAIL failed"
echo "============================================"

[ $FAIL -eq 0 ]
