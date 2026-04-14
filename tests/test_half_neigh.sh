#!/usr/bin/env bash
set -euo pipefail

# Compare full vs half neighbor-list force computation for Verlet Lists.
#
# Verifies that -half 0 (full) and -half 1 (half with Newton's third law)
# produce equivalent thermodynamic output within tolerance.
#
# Usage:
#   ./tests/test_half_neigh.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-NONE}"
DATA_TYPE="${DATA_TYPE:-DP}"
LJ_COMB_RULE="${LJ_COMB_RULE:-geometric}"

cd "${ROOT_DIR}"

echo "Building Verlet-list binary (SIMD=${SIMD}, LJ_COMB_RULE=${LJ_COMB_RULE})..."
make clean >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=verletlist SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" LJ_COMB_RULE="${LJ_COMB_RULE}" >/dev/null

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi
VL_BIN="./MDBench-VL-${TOOL_TAG}-${DATA_TYPE}"

if [[ ! -x "${VL_BIN}" ]]; then
  echo "Binary '${VL_BIN}' is not executable" >&2
  exit 1
fi

FULL_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_full.XXXXXX")"
HALF_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_half.XXXXXX")"

echo "Running with full neighbor lists (-half 0)..."
"${VL_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 200 -half 0 >"${FULL_LOG}"

echo "Running with half neighbor lists (-half 1)..."
"${VL_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 200 -half 1 >"${HALF_LOG}"

get_last_tp() {
  local file="$1"
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${file}" | tail -n 1 || true
}

full_line="$(get_last_tp "${FULL_LOG}")"
half_line="$(get_last_tp "${HALF_LOG}")"

if [[ -z "${full_line}" || -z "${half_line}" ]]; then
  echo "Could not extract thermo lines from outputs." >&2
  echo "Full neigh log: ${FULL_LOG}" >&2
  echo "Half neigh log: ${HALF_LOG}" >&2
  exit 1
fi

full_T=$(echo "${full_line}" | awk '{print $2}')
full_P=$(echo "${full_line}" | awk '{print $3}')
half_T=$(echo "${half_line}" | awk '{print $2}')
half_P=$(echo "${half_line}" | awk '{print $3}')

echo "Full: T=${full_T}, P=${full_P}"
echo "Half: T=${half_T}, P=${half_P}"

python3 - "$full_T" "$half_T" "$full_P" "$half_P" << 'PY'
import sys, math
full_T, half_T, full_P, half_P = map(float, sys.argv[1:])

def rel(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    if a == 0.0:
        return abs(b)
    return abs(b - a) / abs(a)

# Half-neigh accumulates forces in different order, so allow small
# floating-point divergence over 200 steps.
tol_T = 1e-4
tol_P = 1e-3

drift_T = rel(full_T, half_T)
drift_P = rel(full_P, half_P)

print(f"Relative T diff={drift_T:.2e} (tolerance: {tol_T:.2e})")
print(f"Relative P diff={drift_P:.2e} (tolerance: {tol_P:.2e})")

if drift_T > tol_T or drift_P > tol_P:
    sys.stderr.write(f"Half neigh equivalence failed: dT={drift_T:.2e}, dP={drift_P:.2e}\n")
    sys.exit(1)
PY

echo "Full vs half neighbor-list equivalence PASSED."

rm -f "${FULL_LOG}" "${HALF_LOG}"
