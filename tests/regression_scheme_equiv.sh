#!/usr/bin/env bash
set -euo pipefail

# Compare Verlet-list and Cluster-pair schemes on the same argon testcase.
#
# Usage:
#   ./tests/regression_scheme_equiv.sh
#

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-AVX512}"
DATA_TYPE="${DATA_TYPE:-DP}"
CLUSTER_PAIR_KERNEL="${CLUSTER_PAIR_KERNEL:-auto}"
LJ_COMB_RULE="${LJ_COMB_RULE:-geometric}"

cd "${ROOT_DIR}"

echo "Building Verlet-list binary (LJ_COMB_RULE=${LJ_COMB_RULE})..."
make clean >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=verletlist SIMD="${SIMD}" DATA_TYPE="${DATA_TYPE}" LJ_COMB_RULE="${LJ_COMB_RULE}" >/dev/null
if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi
VL_TAG="VL-${TOOL_TAG}-${DATA_TYPE}"
VL_BIN="./MDBench-${VL_TAG}"

echo "Building Cluster-pair binary (LJ_COMB_RULE=${LJ_COMB_RULE})..."
make clean >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=clusterpair SIMD="${SIMD}" DATA_TYPE="${DATA_TYPE}" CLUSTER_PAIR_KERNEL="${CLUSTER_PAIR_KERNEL}" LJ_COMB_RULE="${LJ_COMB_RULE}" >/dev/null
CP_TAG="CP-${CLUSTER_PAIR_KERNEL}-${TOOL_TAG}-${DATA_TYPE}"
CP_BIN="./MDBench-${CP_TAG}"

for bin in "${VL_BIN}" "${CP_BIN}"; do
  if [[ ! -x "${bin}" ]]; then
    echo "Binary '${bin}' is not executable" >&2
    exit 1
  fi
done

VL_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_vl.XXXXXX")"
CP_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_cp.XXXXXX")"

echo "Running Verlet-list argon testcase with ${VL_BIN}"
"${VL_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" -n 200 >"${VL_LOG}"

echo "Running Cluster-pair argon testcase with ${CP_BIN}"
"${CP_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" -n 200 >"${CP_LOG}"

get_last_tp() {
  local file="$1"
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${file}" | tail -n 1 || true
}

vl_line="$(get_last_tp "${VL_LOG}")"
cp_line="$(get_last_tp "${CP_LOG}")"

if [[ -z "${vl_line}" || -z "${cp_line}" ]]; then
  echo "Could not extract thermo lines from outputs." >&2
  echo "VL log: ${VL_LOG}" >&2
  echo "CP log: ${CP_LOG}" >&2
  exit 1
fi

vl_T=$(echo "${vl_line}" | awk '{print $2}')
vl_P=$(echo "${vl_line}" | awk '{print $3}')
cp_T=$(echo "${cp_line}" | awk '{print $2}')
cp_P=$(echo "${cp_line}" | awk '{print $3}')

echo "VL: T=${vl_T}, P=${vl_P}"
echo "CP: T=${cp_T}, P=${cp_P}"

python - "$vl_T" "$cp_T" "$vl_P" "$cp_P" << 'PY'
import sys, math
vl_T, cp_T, vl_P, cp_P = map(float, sys.argv[1:])
def rel(a,b):
    return 0.0 if a == 0.0 else abs(b-a)/abs(a)
tol_T = 0.02
tol_P = 0.05
drift_T = rel(vl_T, cp_T)
drift_P = rel(vl_P, cp_P)
print(f"Relative T diff={drift_T}, relative P diff={drift_P}")
if drift_T > tol_T or drift_P > tol_P:
    sys.stderr.write(f"Scheme equivalence failed: dT={drift_T}, dP={drift_P}\n")
    sys.exit(1)
PY

echo "Verlet-list vs Cluster-pair scheme equivalence regression PASSED (LJ_COMB_RULE=${LJ_COMB_RULE})."

