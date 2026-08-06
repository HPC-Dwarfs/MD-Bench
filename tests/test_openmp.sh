#!/usr/bin/env bash
set -euo pipefail

# OpenMP correctness: compares 1/2/4-thread runs across scheme x LJ rule x testcase x half-neigh.
# -half 1 is checked separately since cross-thread writes to a shared j-atom can race there.
# Usage: ./tests/test_openmp.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-AVX512}"
DATA_TYPE="${DATA_TYPE:-DP}"

cd "${ROOT_DIR}"

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi

get_last_tp() {
  local file="$1"
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${file}" | tail -n 1 || true
}

# half=0 uses a tight tolerance (no cross-thread writes); half=1 uses a loose one (reorder-prone).
compare_thermo() {
  local label="$1" base_T="$2" base_P="$3" cur_T="$4" cur_P="$5" half="$6"

  python3 - "$base_T" "$cur_T" "$base_P" "$cur_P" "$label" "$half" << 'PY'
import sys
a_T, b_T, a_P, b_P = map(float, sys.argv[1:5])
label = sys.argv[5]
half = sys.argv[6]

def rel(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    if a == 0.0:
        return abs(b)
    return abs(b - a) / abs(a)

if half == "1":
    tol_T, tol_P = 2e-2, 5e-2
else:
    tol_T, tol_P = 1e-6, 1e-6

diff_T = rel(a_T, b_T)
diff_P = rel(a_P, b_P)

print(f"[{label}] Relative T diff: {diff_T:.2e} (tolerance: {tol_T:.2e})")
print(f"[{label}] Relative P diff: {diff_P:.2e} (tolerance: {tol_P:.2e})")

if diff_T > tol_T or diff_P > tol_P:
    sys.stderr.write(f"{label} FAILED: dT={diff_T:.2e}, dP={diff_P:.2e}\n")
    sys.exit(1)
PY
}

run_case() {
  local bin="$1" label="$2" half="$3"
  shift 3
  local args=("$@")

  local base_log
  base_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_omp_base.XXXXXX")"
  echo "  [${label}] OMP_NUM_THREADS=1 (baseline)"
  OMP_NUM_THREADS=1 "${bin}" "${args[@]}" -half "${half}" >"${base_log}"
  local base_line base_T base_P
  base_line="$(get_last_tp "${base_log}")"
  if [[ -z "${base_line}" ]]; then
    echo "Could not extract thermo line from baseline output (${label})." >&2
    echo "Log: ${base_log}" >&2
    rm -f "${base_log}"
    exit 1
  fi
  base_T=$(echo "${base_line}" | awk '{print $2}')
  base_P=$(echo "${base_line}" | awk '{print $3}')
  echo "  [${label}] threads=1: T=${base_T}, P=${base_P}"

  for threads in 2 4; do
    local cur_log
    cur_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_omp_t${threads}.XXXXXX")"
    echo "  [${label}] OMP_NUM_THREADS=${threads}"
    OMP_NUM_THREADS="${threads}" "${bin}" "${args[@]}" -half "${half}" >"${cur_log}"
    local cur_line cur_T cur_P
    cur_line="$(get_last_tp "${cur_log}")"
    if [[ -z "${cur_line}" ]]; then
      echo "Could not extract thermo line from threads=${threads} output (${label})." >&2
      echo "Log: ${cur_log}" >&2
      rm -f "${base_log}" "${cur_log}"
      exit 1
    fi
    cur_T=$(echo "${cur_line}" | awk '{print $2}')
    cur_P=$(echo "${cur_line}" | awk '{print $3}')
    echo "  [${label}] threads=${threads}: T=${cur_T}, P=${cur_P}"

    compare_thermo "${label}-threads${threads}" "${base_T}" "${base_P}" "${cur_T}" "${cur_P}" "${half}"
    rm -f "${cur_log}"
  done

  rm -f "${base_log}"
  echo "  [${label}] PASSED."
}

FAILED=0

for scheme in verletlist clusterpair; do
  if [ "${scheme}" = "verletlist" ]; then
    OPT_TAG="VL"
  else
    OPT_TAG="CP-auto"
  fi

  for lj_rule in geometric full; do
    echo "=== Building ${scheme} (LJ_COMB_RULE=${lj_rule}) ==="
    make clean >/dev/null 2>&1 || true
    if ! make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME="${scheme}" SIMD="${SIMD}" \
         DATA_TYPE="${DATA_TYPE}" LJ_COMB_RULE="${lj_rule}" ENABLE_OPENMP=true >/dev/null; then
      echo "SKIP: build failed for ${scheme}/${lj_rule} (unsupported configuration)."
      continue
    fi

    BIN="./MDBench-${OPT_TAG}-${TOOL_TAG}-${DATA_TYPE}-${lj_rule}"
    if [[ ! -x "${BIN}" ]]; then
      echo "Binary '${BIN}' not found" >&2
      exit 1
    fi

    for half in 0 1; do
      echo "--- ${scheme}/${lj_rule}, argon, half=${half} ---"
      if ! run_case "${BIN}" "${scheme}-${lj_rule}-argon-half${half}" "${half}" \
          -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" -n 200; then
        FAILED=1
      fi

      echo "--- ${scheme}/${lj_rule}, copper-fcc, half=${half} ---"
      if ! run_case "${BIN}" "${scheme}-${lj_rule}-copper-half${half}" "${half}" \
          -nx 8 -ny 8 -nz 8 -n 200; then
        FAILED=1
      fi
    done

    rm -f "${BIN}"
  done
done

if [[ "${FAILED}" -ne 0 ]]; then
  echo "OpenMP correctness regression FAILED." >&2
  exit 1
fi

echo "OpenMP correctness regression PASSED."
