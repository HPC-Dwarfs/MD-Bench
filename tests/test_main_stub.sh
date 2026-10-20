#!/usr/bin/env bash
set -euo pipefail

# Smoke-test the synthetic micro-benchmark harnesses (main-stub.c) for both
# optimization schemes.
#
# These harnesses fabricate atoms/neighbor lists directly (no real
# simulation, no reference values), so there is nothing to numerically
# validate here. This test only checks that:
#   - both stub binaries build,
#   - they run to completion (exit 0) across access patterns and, for
#     clusterpair, across the -sc (super-clustering) override,
#   - their output has no NaN/Inf force-kernel garbage.
#
# CPU only: this test always builds with TOOLCHAIN=GCC. There is currently no
# GPU (NVCC/HIPCC) runner in CI, so the super-clustering GPU kernel path that
# main-stub.c's -sc flag targets is not exercised here.
#
# Usage:
#   ./tests/test_main_stub.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

TOOLCHAIN="GCC"
ISA="X86"
DATA_TYPE="DP"

check_output() {
  local log="$1"
  if grep -qiE 'nan|inf|error' "${log}"; then
    echo "Suspicious output in ${log}:" >&2
    cat "${log}" >&2
    exit 1
  fi
}

run_and_check() {
  local bin="$1"
  shift
  local log
  log="$(mktemp "${TMPDIR:-/tmp}/mdbench_stub.XXXXXX")"
  echo "  Running: ${bin} $*"
  "${bin}" "$@" >"${log}" 2>&1
  check_output "${log}"
  rm -f "${log}"
}

echo "=== verletlist stub ==="
SIMD="NONE"
make clean >/dev/null 2>&1 || true
rm -rf build/ >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=verletlist SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" >/dev/null
VL_BIN_TAG="VL-${TOOLCHAIN}-${ISA}-${DATA_TYPE}-${LJ_COMB_RULE:-geometric}"
VL_BIN="./MDBench-${VL_BIN_TAG}"
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=verletlist SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" "${VL_BIN}-stub" >/dev/null

for pattern in seq fix rand; do
  run_and_check "${VL_BIN}-stub" -n 5 -na 64 -nn 8 -p "${pattern}"
done
run_and_check "${VL_BIN}-stub" -n 5 -na 64 -f eam -e data/Cu_u6.eam

echo "=== clusterpair stub ==="
SIMD="AVX2"
make clean >/dev/null 2>&1 || true
rm -rf build/ >/dev/null 2>&1 || true
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=clusterpair SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" >/dev/null
CP_BIN_TAG="CP-${CLUSTER_PAIR_KERNEL:-auto}-${TOOLCHAIN}-${ISA}-${SIMD}-${DATA_TYPE}-${LJ_COMB_RULE:-geometric}"
CP_BIN="./MDBench-${CP_BIN_TAG}"
make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" OPT_SCHEME=clusterpair SIMD="${SIMD}" \
     DATA_TYPE="${DATA_TYPE}" "${CP_BIN}-stub" >/dev/null

for pattern in seq fix; do
  run_and_check "${CP_BIN}-stub" -n 5 -ni 64 -nn 8 -p "${pattern}"
done
run_and_check "${CP_BIN}-stub" -n 5 -ni 64 -nn 8 -p rand

# -sc must be settable on both sides and reflected in the printed banner. Note
# this only checks the CLI override path; it cannot exercise the compile-time
# default (derived from CLUSTERPAIR_KERNEL_GPU_SUPERCLUSTERS), which only
# differs on a GPU super-clustering build not available in this CPU-only test.
sc0_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_sc0.XXXXXX")"
sc1_log="$(mktemp "${TMPDIR:-/tmp}/mdbench_sc1.XXXXXX")"
"${CP_BIN}-stub" -n 5 -ni 64 -sc 0 >"${sc0_log}" 2>&1
"${CP_BIN}-stub" -n 5 -ni 64 -sc 1 >"${sc1_log}" 2>&1
check_output "${sc0_log}"
check_output "${sc1_log}"

if ! grep -q "Super-clustering: no" "${sc0_log}"; then
  echo "-sc 0 did not disable super-clustering:" >&2
  cat "${sc0_log}" >&2
  exit 1
fi
if ! grep -q "Super-clustering: yes" "${sc1_log}"; then
  echo "-sc 1 did not enable super-clustering:" >&2
  cat "${sc1_log}" >&2
  exit 1
fi
rm -f "${sc0_log}" "${sc1_log}"

echo "main-stub smoke test PASSED (CPU only; GPU super-clustering kernel path not covered)."

rm -f "${VL_BIN}" "${VL_BIN}-stub" "${CP_BIN}" "${CP_BIN}-stub"
