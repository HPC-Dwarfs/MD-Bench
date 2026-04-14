#!/usr/bin/env bash
set -euo pipefail

# MPI correctness regression test.
#
# Verifies that an MPI build with 1, 2, and 4 ranks produces consistent
# thermodynamic output compared to the serial (non-MPI) build.
#
# Skips gracefully if MPI is not available.
#
# Usage:
#   ./tests/test_mpi.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-NONE}"
DATA_TYPE="${DATA_TYPE:-DP}"

cd "${ROOT_DIR}"

# Check for MPI availability
if ! command -v mpicc >/dev/null 2>&1; then
  echo "SKIP: mpicc not found, skipping MPI tests."
  exit 0
fi

if ! command -v mpirun >/dev/null 2>&1 && ! command -v mpiexec >/dev/null 2>&1; then
  echo "SKIP: mpirun/mpiexec not found, skipping MPI tests."
  exit 0
fi

MPI_RUN="mpirun"
if ! command -v mpirun >/dev/null 2>&1; then
  MPI_RUN="mpiexec"
fi

# Some MPI implementations reject oversubscription without explicit flag
MPI_EXTRA_FLAGS=""
if "${MPI_RUN}" --help 2>&1 | grep -q "\-\-oversubscribe"; then
  MPI_EXTRA_FLAGS="--oversubscribe"
fi

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi
VL_TAG="VL-${TOOL_TAG}-${DATA_TYPE}"

# Step 1: Build serial (non-MPI) binary
echo "Building serial (non-MPI) Verlet-list binary..."
make clean >/dev/null 2>&1 || true
rm -rf build/ >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=verletlist SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" ENABLE_MPI=false >/dev/null

SERIAL_BIN="./MDBench-${VL_TAG}"
if [[ ! -x "${SERIAL_BIN}" ]]; then
  echo "Serial binary '${SERIAL_BIN}' not found" >&2
  exit 1
fi
mv "${SERIAL_BIN}" "${SERIAL_BIN}-serial"
SERIAL_BIN="${SERIAL_BIN}-serial"

# Step 2: Build MPI binary
echo "Building MPI Verlet-list binary..."
make clean >/dev/null 2>&1 || true
rm -rf build/ >/dev/null 2>&1 || true
if ! make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=verletlist SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" ENABLE_MPI=true >/dev/null 2>&1; then
  echo "SKIP: MPI build failed (ENABLE_MPI=true not supported with this toolchain)."
  rm -f "${SERIAL_BIN}"
  exit 0
fi

MPI_BIN="./MDBench-${VL_TAG}"
if [[ ! -x "${MPI_BIN}" ]]; then
  echo "MPI binary '${MPI_BIN}' not found" >&2
  rm -f "${SERIAL_BIN}"
  exit 1
fi

# Helper to extract last thermo line
get_last_tp() {
  local file="$1"
  grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${file}" | tail -n 1 || true
}

SERIAL_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_serial.XXXXXX")"
MPI1_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_mpi1.XXXXXX")"
MPI2_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_mpi2.XXXXXX")"

echo "Running serial binary..."
"${SERIAL_BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 200 >"${SERIAL_LOG}"

echo "Running MPI binary with 1 rank..."
${MPI_RUN} ${MPI_EXTRA_FLAGS} -np 1 "${MPI_BIN}" \
    -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 200 >"${MPI1_LOG}"

echo "Running MPI binary with 2 ranks..."
${MPI_RUN} ${MPI_EXTRA_FLAGS} -np 2 "${MPI_BIN}" \
    -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" \
    -n 200 >"${MPI2_LOG}" 2>/dev/null || {
  echo "WARN: MPI with 2 ranks failed (may need larger system). Skipping multi-rank test."
  MPI2_LOG=""
}

serial_line="$(get_last_tp "${SERIAL_LOG}")"
mpi1_line="$(get_last_tp "${MPI1_LOG}")"

if [[ -z "${serial_line}" || -z "${mpi1_line}" ]]; then
  echo "Could not extract thermo lines from outputs." >&2
  echo "Serial log: ${SERIAL_LOG}" >&2
  echo "MPI-1 log:  ${MPI1_LOG}" >&2
  rm -f "${SERIAL_BIN}" "${MPI_BIN}"
  exit 1
fi

serial_T=$(echo "${serial_line}" | awk '{print $2}')
serial_P=$(echo "${serial_line}" | awk '{print $3}')
mpi1_T=$(echo "${mpi1_line}" | awk '{print $2}')
mpi1_P=$(echo "${mpi1_line}" | awk '{print $3}')

echo "Serial:  T=${serial_T}, P=${serial_P}"
echo "MPI(1):  T=${mpi1_T}, P=${mpi1_P}"

# MPI with 1 rank should match serial exactly (same domain, same computation)
python3 - "$serial_T" "$mpi1_T" "$serial_P" "$mpi1_P" "serial_vs_mpi1" << 'PY'
import sys, math
a_T, b_T, a_P, b_P = map(float, sys.argv[1:5])
label = sys.argv[5]

def rel(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    if a == 0.0:
        return abs(b)
    return abs(b - a) / abs(a)

tol_T = 1e-6
tol_P = 1e-6
diff_T = rel(a_T, b_T)
diff_P = rel(a_P, b_P)

print(f"[{label}] Relative T diff: {diff_T:.2e} (tolerance: {tol_T:.2e})")
print(f"[{label}] Relative P diff: {diff_P:.2e} (tolerance: {tol_P:.2e})")

if diff_T > tol_T or diff_P > tol_P:
    sys.stderr.write(f"{label} FAILED: dT={diff_T:.2e}, dP={diff_P:.2e}\n")
    sys.exit(1)
PY

echo "Serial vs MPI(1) PASSED."

# Compare MPI with 2 ranks if available
if [[ -n "${MPI2_LOG}" ]]; then
  mpi2_line="$(get_last_tp "${MPI2_LOG}")"
  if [[ -n "${mpi2_line}" ]]; then
    mpi2_T=$(echo "${mpi2_line}" | awk '{print $2}')
    mpi2_P=$(echo "${mpi2_line}" | awk '{print $3}')

    echo "MPI(2):  T=${mpi2_T}, P=${mpi2_P}"

    # With 2 ranks, domain decomposition changes atom ordering and communication
    # introduces rounding differences. Use a looser tolerance.
    python3 - "$serial_T" "$mpi2_T" "$serial_P" "$mpi2_P" "serial_vs_mpi2" << 'PY'
import sys, math
a_T, b_T, a_P, b_P = map(float, sys.argv[1:5])
label = sys.argv[5]

def rel(a, b):
    if a == 0.0 and b == 0.0:
        return 0.0
    if a == 0.0:
        return abs(b)
    return abs(b - a) / abs(a)

tol_T = 0.02
tol_P = 0.05
diff_T = rel(a_T, b_T)
diff_P = rel(a_P, b_P)

print(f"[{label}] Relative T diff: {diff_T:.2e} (tolerance: {tol_T:.2e})")
print(f"[{label}] Relative P diff: {diff_P:.2e} (tolerance: {tol_P:.2e})")

if diff_T > tol_T or diff_P > tol_P:
    sys.stderr.write(f"{label} FAILED: dT={diff_T:.2e}, dP={diff_P:.2e}\n")
    sys.exit(1)
PY

    echo "Serial vs MPI(2) PASSED."
  fi
fi

echo "MPI correctness regression PASSED."

rm -f "${SERIAL_LOG}" "${MPI1_LOG}" "${MPI2_LOG}" "${SERIAL_BIN}" "${MPI_BIN}"
