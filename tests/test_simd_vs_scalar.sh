#!/usr/bin/env bash
set -euo pipefail

# Compare SIMD vs SCALAR kernel implementations for Verlet Lists
# Verifies:
#   1. Correct kernel is selected (SIMD vs SCALAR)
#   2. Temperature and pressure match within tolerance

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

# Build configuration
TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-AVX2}"
DATA_TYPE="${DATA_TYPE:-DP}"

cd "${ROOT_DIR}"

# Build SIMD variant
echo "Building SIMD kernel binary..."
rm -rf build/ >/dev/null 2>&1 || true  # Force clean build directory
make clean >/dev/null 2>&1 || true
if ! make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" OPT_SCHEME=verletlist \
     USE_SIMD_KERNEL=true >/dev/null 2>&1; then
  echo "ERROR: Failed to build SIMD binary" >&2
  make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" \
       DATA_TYPE="${DATA_TYPE}" OPT_SCHEME=verletlist \
       USE_SIMD_KERNEL=true
  exit 1
fi

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi

SIMD_BIN="./MDBench-VL-${TOOL_TAG}-${DATA_TYPE}"

# Rename SIMD binary to preserve it before building SCALAR variant
if [[ ! -x "${SIMD_BIN}" ]]; then
  echo "Binary '${SIMD_BIN}' is not executable after SIMD build" >&2
  exit 1
fi
mv "${SIMD_BIN}" "${SIMD_BIN}-simd"
SIMD_BIN="${SIMD_BIN}-simd"

# Build SCALAR variant
echo "Building SCALAR kernel binary..."
rm -rf build/ >/dev/null 2>&1 || true  # Force clean build directory
make clean >/dev/null 2>&1 || true
if ! make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" OPT_SCHEME=verletlist \
     USE_SIMD_KERNEL=false >/dev/null 2>&1; then
  echo "ERROR: Failed to build SCALAR binary" >&2
  make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" \
       DATA_TYPE="${DATA_TYPE}" OPT_SCHEME=verletlist \
       USE_SIMD_KERNEL=false
  exit 1
fi

SCALAR_BIN="./MDBench-VL-${TOOL_TAG}-${DATA_TYPE}"
# Rename to avoid confusion with SIMD binary
if [[ ! -x "${SCALAR_BIN}" ]]; then
  echo "Binary '${SCALAR_BIN}' is not executable after SCALAR build" >&2
  exit 1
fi
mv "${SCALAR_BIN}" "${SCALAR_BIN}-scalar"
SCALAR_BIN="${SCALAR_BIN}-scalar"

# Verify binaries exist
for bin in "${SIMD_BIN}" "${SCALAR_BIN}"; do
  if [[ ! -x "${bin}" ]]; then
    echo "Binary '${bin}' is not executable" >&2
    exit 1
  fi
done

# Create temp files for output
SIMD_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_simd.XXXXXX")"
SCALAR_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_scalar.XXXXXX")"

# Run SIMD variant
echo "Running SIMD kernel..."
"${SIMD_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" -n 500 >"${SIMD_LOG}"

# Run SCALAR variant
echo "Running SCALAR kernel..."
"${SCALAR_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" -n 500 >"${SCALAR_LOG}"

# Verify kernel selection
simd_kernel=$(grep "Kernel:" "${SIMD_LOG}" | grep -v "Computational" | awk '{print $2}')
scalar_kernel=$(grep "Kernel:" "${SCALAR_LOG}" | grep -v "Computational" | awk '{print $2}')

echo "SIMD binary reports kernel: ${simd_kernel}"
echo "SCALAR binary reports kernel: ${scalar_kernel}"

if [[ "${simd_kernel}" != "SIMD" ]]; then
  echo "ERROR: SIMD binary did not use SIMD kernel (got: ${simd_kernel})" >&2
  exit 1
fi

if [[ "${scalar_kernel}" != "SCALAR" ]]; then
  echo "ERROR: SCALAR binary did not use SCALAR kernel (got: ${scalar_kernel})" >&2
  exit 1
fi

echo "✓ Kernel selection verified"

# Extract final T and P values
get_last_tp() {
  local file="$1"
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${file}" | tail -n 1 || true
}

simd_line="$(get_last_tp "${SIMD_LOG}")"
scalar_line="$(get_last_tp "${SCALAR_LOG}")"

if [[ -z "${simd_line}" || -z "${scalar_line}" ]]; then
  echo "Could not extract thermo lines from outputs." >&2
  echo "SIMD log: ${SIMD_LOG}" >&2
  echo "SCALAR log: ${SCALAR_LOG}" >&2
  exit 1
fi

simd_T=$(echo "${simd_line}" | awk '{print $2}')
simd_P=$(echo "${simd_line}" | awk '{print $3}')
scalar_T=$(echo "${scalar_line}" | awk '{print $2}')
scalar_P=$(echo "${scalar_line}" | awk '{print $3}')

echo "SIMD:   T=${simd_T}, P=${simd_P}"
echo "SCALAR: T=${scalar_T}, P=${scalar_P}"

# Compare with tolerance
python - "$simd_T" "$scalar_T" "$simd_P" "$scalar_P" << 'PY'
import sys, math
simd_T, scalar_T, simd_P, scalar_P = map(float, sys.argv[1:])

def rel_diff(a, b):
    """Calculate relative difference between a and b"""
    if a == 0.0 and b == 0.0:
        return 0.0
    if a == 0.0:
        return abs(b)
    return abs(b - a) / abs(a)

# Tolerances - SIMD and SCALAR should be very close
# Allow small differences due to instruction reordering, rounding
tol_T = 1e-6  # 0.0001% for temperature
tol_P = 1e-6  # 0.0001% for pressure

diff_T = rel_diff(scalar_T, simd_T)
diff_P = rel_diff(scalar_P, simd_P)

print(f"Relative T diff: {diff_T:.2e} (tolerance: {tol_T:.2e})")
print(f"Relative P diff: {diff_P:.2e} (tolerance: {tol_P:.2e})")

if diff_T > tol_T or diff_P > tol_P:
    sys.stderr.write(f"SIMD vs SCALAR mismatch: dT={diff_T:.2e}, dP={diff_P:.2e}\n")
    sys.exit(1)
PY

echo "✓ Temperature and pressure match within tolerance"
echo ""
echo "SIMD vs SCALAR kernel validation PASSED"

# Cleanup
rm -f "${SIMD_LOG}" "${SCALAR_LOG}" "${SCALAR_BIN}" "${SIMD_BIN}"
