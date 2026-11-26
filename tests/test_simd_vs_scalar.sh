#!/usr/bin/env bash
set -euo pipefail

# Compare the verlet-list intrinsics kernels and neighbor-list build against the plain C
# (autovec) ones of the same binary (--force-simd, --neighbor-simd). Verifies:
#   1. The selected kernel is reported
#   2. Temperature and pressure match within tolerance

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

# Build configuration
TOOLCHAIN="${TOOLCHAIN:-GCC}"
SIMD="${SIMD:-AVX2}"
DATA_TYPE="${DATA_TYPE:-DP}"
# LJ combination rule of the runs (--lj-comb)
LJ_COMB="${LJ_COMB:-geometric}"

cd "${ROOT_DIR}"

CFG=(TOOLCHAIN="${TOOLCHAIN}" SIMD="${SIMD}" DATA_TYPE="${DATA_TYPE}" OPT_SCHEME=verletlist)
echo "Building verlet-list binary (SIMD=${SIMD})..."
if ! make -j "${CFG[@]}" >/dev/null 2>&1; then
  echo "ERROR: Failed to build the binary" >&2
  make "${CFG[@]}"
  exit 1
fi
BIN="$(make -s "${CFG[@]}" print-target)"
if [[ ! -x "${BIN}" ]]; then
  echo "Binary '${BIN}' is not executable" >&2
  exit 1
fi

LOG_DIR="$(mktemp -d "${TMPDIR:-/tmp}/mdbench_simd.XXXXXX")"
trap 'rm -rf "${LOG_DIR}"' EXIT

# Runs one mode, checks the reported kernel and prints the final "T P"
run_mode() {
  local name="$1" kernel="$2"
  shift 2
  local log="${LOG_DIR}/${name}.log"
  "${BIN}" --lj-comb "${LJ_COMB}" -i "${DATA_DIR}/input.gro" \
    -p "${DATA_DIR}/mdbench_params.conf" -n 500 "$@" >"${log}"
  local reported
  reported=$(grep "Kernel:" "${log}" | grep -v "Computational" | awk '{print $2}')
  if [[ "${reported}" != "${kernel}" ]]; then
    echo "ERROR: ${name} did not use the ${kernel} kernel (got: ${reported})" >&2
    exit 1
  fi
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${log}" | tail -n 1 |
    awk '{print $2, $3}'
}

ref="$(run_mode autovec autovec --force-simd autovec)"
echo "autovec:                 T P = ${ref}"
status=0
for mode in "intrinsics intrinsics --force-simd intrinsics" \
  "intrinsics-compress intrinsics-compress --force-simd intrinsics-compress" \
  "neighbor-intrinsics autovec --neighbor-simd intrinsics"; do
  # shellcheck disable=SC2086
  set -- ${mode}
  name="$1"
  out="$(run_mode "$@")"
  printf "%-24s T P = %s\n" "${name}:" "${out}"
  # SIMD and plain C should be very close (instruction reordering, rounding)
  if ! python3 - ${ref} ${out} <<'PY'; then
import sys
ref_T, ref_P, T, P = map(float, sys.argv[1:])
rel = lambda a, b: 0.0 if a == b else abs(b - a) / max(abs(a), 1e-30)
dT, dP = rel(ref_T, T), rel(ref_P, P)
print(f"  relative T diff {dT:.2e}, P diff {dP:.2e} (tolerance 1e-6)")
sys.exit(dT > 1e-6 or dP > 1e-6)
PY
    echo "ERROR: ${name} does not match autovec" >&2
    status=1
  fi
done

if [[ ${status} -ne 0 ]]; then
  exit 1
fi
echo ""
echo "Intrinsics vs autovec validation PASSED"
