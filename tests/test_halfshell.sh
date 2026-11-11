#!/usr/bin/env bash
set -euo pipefail

# Forces of the cluster-pair kernels with half neighbor lists have to match the ones with
# full neighbor lists. Without MPI a half list keeps only half of the periodic ghost
# clusters (half shell) and the kernels write their reaction forces into the remaining
# ghost clusters, which reverse() folds back onto the real atoms.
#
# Compares the per-atom forces of one kernel pass (--dump-forces) of the argon gas and
# the GROMACS benchmark water system for the 4xN, 2xNN and 2xN kernels, with a single thread (half lists are not thread safe).
#
# Needs python3 with numpy, the test is skipped otherwise.
#
# Usage:
#   ./tests/test_halfshell.sh
#   SIMD=AVX2 ./tests/test_halfshell.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-AVX2}"
PYTHON="${PYTHON:-python3}"

if ! "${PYTHON}" -c "import numpy" >/dev/null 2>&1; then
  echo "SKIP: ${PYTHON} with numpy is required"
  exit 0
fi

cd "${ROOT_DIR}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mdbench_halfshell.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi

status=0
# kernel:data type:tolerance (rms error relative to the rms force)
CONFIGS=("4xN:DP:1e-9")
case "${SIMD}" in
  AVX2) CONFIGS+=("2xNN:SP:1e-5" "2xN:DP:1e-9") ;;
  AVX512) CONFIGS+=("2xNN:SP:1e-5") ;;
esac
for config in "${CONFIGS[@]}"; do
  IFS=: read -r kernel dtype tol <<<"${config}"
  make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" DATA_TYPE="${dtype}" \
       OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL="${kernel}" ACCURATE_RECIPROCAL=true \
       >/dev/null
  bin="./MDBench-CP-${kernel}-${TOOL_TAG}-${dtype}-geometric-accrecip"
  for case in argon gmx; do
    if [ "${case}" = "argon" ]; then
      args=(-i data/argon/input.gro -p data/argon/mdbench_params.conf)
    else
      args=(-gmxbench -size 1)
    fi
    for half in 0 1; do
      OMP_NUM_THREADS=1 "${bin}" "${args[@]}" -half "${half}" --kernel-bench 1 \
          --dump-forces "${WORK}/f_${half}.txt" >/dev/null
    done
    if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/f_0.txt" "${WORK}/f_1.txt" \
         --tol "${tol}" --tol-max "$(python3 -c "print(${tol}*10)")" \
         --label "${kernel} ${dtype} ${case} half vs full"; then
      status=1
    fi
  done
done

if [ "${status}" -eq 0 ]; then
  echo "Half vs full neighbor list force equivalence PASSED"
else
  echo "Half vs full neighbor list force equivalence FAILED" >&2
fi
exit "${status}"
