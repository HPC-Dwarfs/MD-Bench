#!/usr/bin/env bash
set -euo pipefail

# Compare the explicit (i,j)/(ci,cj) pair-list neighbor representation
# (NBLIST_DATA_LAYOUT=PAIRLIST) against the standard per-owner AOS list.
#
# Verifies:
#   1. PAIRLIST and AOS produce equivalent thermodynamic output.
#   2. PAIRLIST results are consistent across OpenMP thread counts (this is
#      the sharpest check on the per-thread force-buffer + reduction path
#      used to let multiple threads contribute partial force to the same
#      atom/cluster).
#
# Usage:
#   OPT_SCHEME=verletlist ./tests/test_pairlist.sh
#   OPT_SCHEME=clusterpair ./tests/test_pairlist.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
DATA_TYPE="${DATA_TYPE:-DP}"
OPT_SCHEME="${OPT_SCHEME:-verletlist}"

# Clusterpair requires SIMD != NONE (CLUSTER_N must be 2, 4 or 8).
SIMD="${SIMD:-AVX2}"
# Verletlist's PAIRLIST mode dispatches to the SIMD mixed-lane kernel when
# USE_SIMD_KERNEL=true, and to the scalar reference kernel otherwise;
# clusterpair's PAIRLIST mode only has the scalar reference kernel so this
# is a no-op there. Default to exercising the SIMD kernel.
USE_SIMD_KERNEL="${USE_SIMD_KERNEL:-true}"

cd "${ROOT_DIR}"

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

echo "Building ${OPT_SCHEME} with AOS layout..."
make clean >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME="${OPT_SCHEME}" SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" NBLIST_DATA_LAYOUT=AOS >/dev/null
AOS_BIN="./MDBench-${BIN_TAG}"
if [[ ! -x "${AOS_BIN}" ]]; then
  echo "AOS binary '${AOS_BIN}' not found or not executable" >&2
  exit 1
fi
mv "${AOS_BIN}" "${AOS_BIN}-aos"
AOS_BIN="${AOS_BIN}-aos"

echo "Building ${OPT_SCHEME} with PAIRLIST layout (USE_SIMD_KERNEL=${USE_SIMD_KERNEL})..."
make clean >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME="${OPT_SCHEME}" SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" NBLIST_DATA_LAYOUT=PAIRLIST USE_SIMD_KERNEL="${USE_SIMD_KERNEL}" >/dev/null
PL_BIN="./MDBench-${BIN_TAG}"
if [[ ! -x "${PL_BIN}" ]]; then
  echo "PAIRLIST binary '${PL_BIN}' not found or not executable" >&2
  exit 1
fi
mv "${PL_BIN}" "${PL_BIN}-pairlist"
PL_BIN="${PL_BIN}-pairlist"

AOS_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_aos.XXXXXX")"
PL1_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_pl1.XXXXXX")"
PL4_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_pl4.XXXXXX")"

echo "Running AOS variant..."
"${AOS_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 500 >"${AOS_LOG}"

echo "Running PAIRLIST variant (1 thread)..."
OMP_NUM_THREADS=1 "${PL_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 500 >"${PL1_LOG}"

echo "Running PAIRLIST variant (4 threads)..."
OMP_NUM_THREADS=4 "${PL_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 500 >"${PL4_LOG}"

get_last_tp() {
  local file="$1"
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${file}" | tail -n 1 || true
}

aos_line="$(get_last_tp "${AOS_LOG}")"
pl1_line="$(get_last_tp "${PL1_LOG}")"
pl4_line="$(get_last_tp "${PL4_LOG}")"

if [[ -z "${aos_line}" || -z "${pl1_line}" || -z "${pl4_line}" ]]; then
  echo "Could not extract thermo lines from outputs." >&2
  echo "AOS log: ${AOS_LOG}" >&2
  echo "PAIRLIST(1) log: ${PL1_LOG}" >&2
  echo "PAIRLIST(4) log: ${PL4_LOG}" >&2
  exit 1
fi

aos_T=$(echo "${aos_line}" | awk '{print $2}')
aos_P=$(echo "${aos_line}" | awk '{print $3}')
pl1_T=$(echo "${pl1_line}" | awk '{print $2}')
pl1_P=$(echo "${pl1_line}" | awk '{print $3}')
pl4_T=$(echo "${pl4_line}" | awk '{print $2}')
pl4_P=$(echo "${pl4_line}" | awk '{print $3}')

echo "AOS:            T=${aos_T}, P=${aos_P}"
echo "PAIRLIST (1t):  T=${pl1_T}, P=${pl1_P}"
echo "PAIRLIST (4t):  T=${pl4_T}, P=${pl4_P}"

python3 - "$aos_T" "$pl1_T" "$aos_P" "$pl1_P" "$pl4_T" "$pl4_P" << 'PY'
import sys
aos_T, pl1_T, aos_P, pl1_P, pl4_T, pl4_P = map(float, sys.argv[1:])

def rel(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    if a == 0.0:
        return abs(b)
    return abs(b - a) / abs(a)

# The per-thread force-buffer reduction sums in a different order than the
# direct in-place AOS kernel, so allow a little more slack than the bit-exact
# AOS-vs-SOA data-layout test, but still tight enough to catch a real bug.
tol = 1e-4

d_aos_pl1_T = rel(aos_T, pl1_T)
d_aos_pl1_P = rel(aos_P, pl1_P)
d_pl1_pl4_T = rel(pl1_T, pl4_T)
d_pl1_pl4_P = rel(pl1_P, pl4_P)

print(f"AOS vs PAIRLIST(1t):     dT={d_aos_pl1_T:.2e}, dP={d_aos_pl1_P:.2e}")
print(f"PAIRLIST(1t) vs (4t):    dT={d_pl1_pl4_T:.2e}, dP={d_pl1_pl4_P:.2e}")

if max(d_aos_pl1_T, d_aos_pl1_P, d_pl1_pl4_T, d_pl1_pl4_P) > tol:
    sys.stderr.write("PAIRLIST equivalence/thread-invariance check FAILED\n")
    sys.exit(1)
PY

echo "PAIRLIST equivalence and thread-count invariance PASSED (${OPT_SCHEME})."

rm -f "${AOS_LOG}" "${PL1_LOG}" "${PL4_LOG}" "${AOS_BIN}" "${PL_BIN}"
