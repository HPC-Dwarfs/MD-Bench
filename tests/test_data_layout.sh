#!/usr/bin/env bash
set -euo pipefail

# Compare AOS vs SOA atom data layouts.
#
# Verifies that ATOM_DATA_LAYOUT=AOS and ATOM_DATA_LAYOUT=SOA produce
# identical thermodynamic output (positions are accessed via macros,
# so the layout change should be transparent to the physics).
#
# Usage:
#   ./tests/test_data_layout.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
DATA_TYPE="${DATA_TYPE:-DP}"
OPT_SCHEME="${OPT_SCHEME:-verletlist}"

# Clusterpair requires SIMD != NONE (CLUSTER_N must be 2, 4 or 8).
# Default to NONE for verletlist, AVX2 for clusterpair.
if [ "${OPT_SCHEME}" = "clusterpair" ]; then
  SIMD="${SIMD:-AVX2}"
else
  SIMD="${SIMD:-NONE}"
fi

cd "${ROOT_DIR}"

echo "Building ${OPT_SCHEME} with AOS layout..."
make clean >/dev/null 2>&1 || true
rm -rf build/ >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME="${OPT_SCHEME}" SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" ATOM_DATA_LAYOUT=AOS >/dev/null

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi
if [ "${OPT_SCHEME}" = "verletlist" ]; then
  OPT_TAG="VL"
else
  OPT_TAG="CP-${CLUSTER_PAIR_KERNEL:-auto}"
fi

BIN_TAG="${OPT_TAG}-${TOOL_TAG}-${DATA_TYPE}"
AOS_BIN="./MDBench-${BIN_TAG}"

if [[ ! -x "${AOS_BIN}" ]]; then
  echo "AOS binary '${AOS_BIN}' not found or not executable" >&2
  exit 1
fi
mv "${AOS_BIN}" "${AOS_BIN}-aos"
AOS_BIN="${AOS_BIN}-aos"

echo "Building ${OPT_SCHEME} with SOA layout..."
make clean >/dev/null 2>&1 || true
rm -rf build/ >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME="${OPT_SCHEME}" SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" ATOM_DATA_LAYOUT=SOA >/dev/null

SOA_BIN="./MDBench-${BIN_TAG}"
if [[ ! -x "${SOA_BIN}" ]]; then
  echo "SOA binary '${SOA_BIN}' not found or not executable" >&2
  exit 1
fi
mv "${SOA_BIN}" "${SOA_BIN}-soa"
SOA_BIN="${SOA_BIN}-soa"

AOS_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_aos.XXXXXX")"
SOA_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_soa.XXXXXX")"

echo "Running AOS variant..."
"${AOS_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 200 >"${AOS_LOG}"

echo "Running SOA variant..."
"${SOA_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 200 >"${SOA_LOG}"

get_last_tp() {
  local file="$1"
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${file}" | tail -n 1 || true
}

aos_line="$(get_last_tp "${AOS_LOG}")"
soa_line="$(get_last_tp "${SOA_LOG}")"

if [[ -z "${aos_line}" || -z "${soa_line}" ]]; then
  echo "Could not extract thermo lines from outputs." >&2
  echo "AOS log: ${AOS_LOG}" >&2
  echo "SOA log: ${SOA_LOG}" >&2
  exit 1
fi

aos_T=$(echo "${aos_line}" | awk '{print $2}')
aos_P=$(echo "${aos_line}" | awk '{print $3}')
soa_T=$(echo "${soa_line}" | awk '{print $2}')
soa_P=$(echo "${soa_line}" | awk '{print $3}')

echo "AOS: T=${aos_T}, P=${aos_P}"
echo "SOA: T=${soa_T}, P=${soa_P}"

# AOS and SOA should produce bit-identical results (same operations, same order).
python3 - "$aos_T" "$soa_T" "$aos_P" "$soa_P" << 'PY'
import sys, math
aos_T, soa_T, aos_P, soa_P = map(float, sys.argv[1:])

def rel(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    if a == 0.0:
        return abs(b)
    return abs(b - a) / abs(a)

# AOS and SOA should be identical (or nearly so).
tol_T = 1e-6
tol_P = 1e-6

diff_T = rel(aos_T, soa_T)
diff_P = rel(aos_P, soa_P)

print(f"Relative T diff: {diff_T:.2e} (tolerance: {tol_T:.2e})")
print(f"Relative P diff: {diff_P:.2e} (tolerance: {tol_P:.2e})")

if diff_T > tol_T or diff_P > tol_P:
    sys.stderr.write(f"AOS vs SOA mismatch: dT={diff_T:.2e}, dP={diff_P:.2e}\n")
    sys.exit(1)
PY

echo "AOS vs SOA data layout equivalence PASSED (${OPT_SCHEME})."

rm -f "${AOS_LOG}" "${SOA_LOG}" "${AOS_BIN}" "${SOA_BIN}"
