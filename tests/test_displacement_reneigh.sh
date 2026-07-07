#!/usr/bin/env bash
set -euo pipefail

# Test displacement-based neighbor-list rebuild strategy.
#
# Three checks, each for both Verlet-list and Cluster-pair schemes:
#
#   1. Equivalence: with a large skin (displacement threshold never crossed),
#      enabling --displacement-reneigh produces identical results to fixed-
#      interval-only mode.
#
#   2. Stability: pure displacement mode (--reneigh-every 0) runs to
#      completion and produces physically plausible thermodynamics (T > 0).
#
#   3. Combined mode: both strategies enabled simultaneously, results stay
#      within tolerances compared to the fixed-interval reference.
#
# Usage:
#   ./tests/test_displacement_reneigh.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-NONE}"
# Clusterpair requires a SIMD width that maps to a valid cluster N (2, 4, or 8).
# Default to AVX2 (N=4) when the VL SIMD is NONE.
SIMD_CP="${SIMD_CP:-AVX2}"
DATA_TYPE="${DATA_TYPE:-DP}"
LJ_COMB_RULE="${LJ_COMB_RULE:-geometric}"

cd "${ROOT_DIR}"

# ── helpers ──────────────────────────────────────────────────────────────────

make_bin() {
    local scheme="$1"
    local simd="$2"
    make clean >/dev/null 2>&1 || true
    make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME="${scheme}" \
         SIMD="${simd}" DATA_TYPE="${DATA_TYPE}" \
         LJ_COMB_RULE="${LJ_COMB_RULE}" >/dev/null
}

bin_name() {
    local prefix="$1"
    local simd="$2"
    local suffix="${3:-}"
    if [ "${simd}" = "NONE" ]; then
        echo "./MDBench-${prefix}-${TOOLCHAIN}-${ISA}-${DATA_TYPE}${suffix}"
    else
        echo "./MDBench-${prefix}-${TOOLCHAIN}-${ISA}-${simd}-${DATA_TYPE}${suffix}"
    fi
}

get_last_tp() {
    grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "$1" | tail -n 1 || true
}

# Compare two T/P pairs; print diffs and fail if beyond tolerance.
check_equiv() {
    local label="$1" ref_T="$2" ref_P="$3" got_T="$4" got_P="$5"
    local tol_T="${6:-1e-6}" tol_P="${7:-1e-6}"
    python3 - "${ref_T}" "${ref_P}" "${got_T}" "${got_P}" \
             "${tol_T}" "${tol_P}" "${label}" << 'PY'
import sys, math
ref_T, ref_P, got_T, got_P = map(float, sys.argv[1:5])
tol_T, tol_P = map(float, sys.argv[5:7])
label = sys.argv[7]

def rel(a, b):
    return 0.0 if a == 0.0 and b == 0.0 else abs(b - a) / (abs(a) if a != 0.0 else 1.0)

dT = rel(ref_T, got_T)
dP = rel(ref_P, got_P)
print(f"  {label}: dT={dT:.2e} (tol {tol_T:.0e}), dP={dP:.2e} (tol {tol_P:.0e})")
if dT > tol_T or dP > tol_P:
    sys.stderr.write(f"FAILED {label}: dT={dT:.2e}, dP={dP:.2e}\n")
    sys.exit(1)
PY
}

check_plausible() {
    local label="$1" T="$2" P="$3"
    python3 - "${T}" "${P}" "${label}" << 'PY'
import sys, math
T, P = map(float, sys.argv[1:3])
label = sys.argv[3]
if not (math.isfinite(T) and math.isfinite(P) and T > 0):
    sys.stderr.write(f"FAILED {label}: T={T}, P={P} not physically plausible\n")
    sys.exit(1)
print(f"  {label}: T={T:.4g}, P={P:.4g} - plausible")
PY
}

# ── write a param file with a given reneigh_every ────────────────────────────

make_params() {
    local reneigh_every="$1"
    local skin="$2"
    local out
    out="$(mktemp "${TMPDIR:-/tmp}/mdbench_params.XXXXXX")"
    sed -e "s/^reneigh_every .*/reneigh_every ${reneigh_every}/" \
        -e "s/^skin .*/skin ${skin}/" \
        "${DATA_DIR}/mdbench_params.conf" > "${out}"
    echo "${out}"
}

NSTEPS=200
# Fewer steps than reneigh_every (100) and short enough that typical argon
# displacements (< 0.02 nm over 20 steps) stay below skin/2 = 0.05 nm.
FEW_STEPS=20
NORMAL_SKIN=0.1   # argon default

PASS=0
FAIL=0

run_scheme() {
    local scheme="$1"
    local prefix="$2"
    local simd="$3"
    local suffix="${4:-}"
    local bin

    echo ""
    echo "=== Scheme: ${scheme} (SIMD=${simd}) ==="

    make_bin "${scheme}" "${simd}"
    bin="$(bin_name "${prefix}" "${simd}" "${suffix}")"

    if [[ ! -x "${bin}" ]]; then
        echo "Binary '${bin}' not found or not executable" >&2
        FAIL=$((FAIL + 1))
        return
    fi

    # ── Check 1: equivalence when neither reneigh trigger fires ─────────────
    # Use FEW_STEPS < reneigh_every (100) so the fixed-interval trigger never
    # fires, and the run is short enough that argon displacements stay well
    # below skin/2 so the displacement trigger also never fires.
    # Both runs must therefore produce bit-for-bit identical output.
    echo "--- Check 1: equivalence (no rebuild triggered in either run) ---"

    params_norm_check1="$(make_params 100 "${NORMAL_SKIN}")"

    ref_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_ref.XXXXXX")"
    disp_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_disp.XXXXXX")"

    "${bin}" -i "${DATA_DIR}/input.gro" -p "${params_norm_check1}" -n "${FEW_STEPS}" \
        > "${ref_log}"
    "${bin}" -i "${DATA_DIR}/input.gro" -p "${params_norm_check1}" -n "${FEW_STEPS}" \
        --displacement-reneigh > "${disp_log}"

    ref_line="$(get_last_tp "${ref_log}")"
    disp_line="$(get_last_tp "${disp_log}")"

    if [[ -z "${ref_line}" || -z "${disp_line}" ]]; then
        echo "Could not extract thermo output" >&2
        FAIL=$((FAIL + 1))
    else
        ref_T=$(echo "${ref_line}" | awk '{print $2}')
        ref_P=$(echo "${ref_line}" | awk '{print $3}')
        disp_T=$(echo "${disp_line}" | awk '{print $2}')
        disp_P=$(echo "${disp_line}" | awk '{print $3}')

        if check_equiv "${scheme} no-rebuild equivalence" \
                       "${ref_T}" "${ref_P}" "${disp_T}" "${disp_P}" 0 0; then
            echo "  PASSED"
            PASS=$((PASS + 1))
        else
            FAIL=$((FAIL + 1))
        fi
    fi

    rm -f "${ref_log}" "${disp_log}" "${params_norm_check1}"

    # ── Check 2: pure displacement mode (reneigh_every=0) stability ──────────
    echo "--- Check 2: pure displacement mode stability ---"

    params_disp="$(make_params 0 "${NORMAL_SKIN}")"
    stab_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_stab.XXXXXX")"

    "${bin}" -i "${DATA_DIR}/input.gro" -p "${params_disp}" -n "${NSTEPS}" \
        --displacement-reneigh > "${stab_log}"

    stab_line="$(get_last_tp "${stab_log}")"
    if [[ -z "${stab_line}" ]]; then
        echo "Could not extract thermo output from pure-displacement run" >&2
        FAIL=$((FAIL + 1))
    else
        stab_T=$(echo "${stab_line}" | awk '{print $2}')
        stab_P=$(echo "${stab_line}" | awk '{print $3}')

        if check_plausible "${scheme} pure-displacement" "${stab_T}" "${stab_P}"; then
            echo "  PASSED"
            PASS=$((PASS + 1))
        else
            FAIL=$((FAIL + 1))
        fi
    fi

    rm -f "${stab_log}" "${params_disp}"

    # ── Check 3: combined mode stays within loose physical tolerance ──────────
    echo "--- Check 3: combined mode (fixed-interval + displacement) ---"

    params_norm="$(make_params 20 "${NORMAL_SKIN}")"
    ref2_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_ref2.XXXXXX")"
    comb_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_comb.XXXXXX")"

    "${bin}" -i "${DATA_DIR}/input.gro" -p "${params_norm}" -n "${NSTEPS}" \
        > "${ref2_log}"
    "${bin}" -i "${DATA_DIR}/input.gro" -p "${params_norm}" -n "${NSTEPS}" \
        --displacement-reneigh > "${comb_log}"

    ref2_line="$(get_last_tp "${ref2_log}")"
    comb_line="$(get_last_tp "${comb_log}")"

    if [[ -z "${ref2_line}" || -z "${comb_line}" ]]; then
        echo "Could not extract thermo output from combined-mode run" >&2
        FAIL=$((FAIL + 1))
    else
        ref2_T=$(echo "${ref2_line}" | awk '{print $2}')
        ref2_P=$(echo "${ref2_line}" | awk '{print $3}')
        comb_T=$(echo "${comb_line}" | awk '{print $2}')
        comb_P=$(echo "${comb_line}" | awk '{print $3}')

        # Combined mode may rebuild more often → trajectory can diverge; use loose tol.
        if check_equiv "${scheme} combined mode" \
                       "${ref2_T}" "${ref2_P}" "${comb_T}" "${comb_P}" 0.05 0.10; then
            echo "  PASSED"
            PASS=$((PASS + 1))
        else
            FAIL=$((FAIL + 1))
        fi
    fi

    rm -f "${ref2_log}" "${comb_log}" "${params_norm}"
}

run_scheme "verletlist"  "VL"      "${SIMD}"    "-${LJ_COMB_RULE}"
run_scheme "clusterpair" "CP-auto" "${SIMD_CP}"

echo ""
echo "Results: ${PASS} passed, ${FAIL} failed."

if [ "${FAIL}" -gt 0 ]; then
    exit 1
fi

echo "Displacement reneigh tests PASSED."
