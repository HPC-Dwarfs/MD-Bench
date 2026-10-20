#!/usr/bin/env bash
set -euo pipefail

# Validates the electrostatics (reaction-field, Ewald) and molecular exclusions of the
# cluster-pair 4xN kernel against an exact float64 reference.
#
# The system is the water box of the GROMACS non-bonded benchmark (-gmxbench). The
# reference (tests/gromacs_nonbonded/nbnxm_reference.py) evaluates the same pair
# forces with erf/erfc, so it shares no approximation with the kernels. Both full and
# half neighbor lists are checked, in single and double precision.
#
# Needs python3 with numpy and scipy, the test is skipped otherwise.
#
# Usage:
#   ./tests/test_coulomb_exclusions.sh
#   SIMD=AVX2 ./tests/test_coulomb_exclusions.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-AVX2}"
PYTHON="${PYTHON:-python3}"

if ! "${PYTHON}" -c "import numpy, scipy" >/dev/null 2>&1; then
  echo "SKIP: ${PYTHON} with numpy and scipy is required"
  exit 0
fi

cd "${ROOT_DIR}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mdbench_coulomb.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi

status=0
for rule in rf ewald none; do
  "${PYTHON}" "${DIR}/nbnxm_reference.py" --size 1 --coulomb "${rule}" --out "${WORK}/ref_${rule}.txt"
  for dtype in DP SP; do
    if [ "${dtype}" = "DP" ]; then tol="1e-7"; else tol="5e-5"; fi
    excl=auto
    suffix="-${rule}-accrecip"
    if [ "${rule}" = "none" ]; then excl=true; suffix="-accrecip-excl"; fi
    make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" DATA_TYPE="${dtype}" \
         OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL=4xN COULOMB_RULE="${rule}" \
         EXCLUSIONS="${excl}" ACCURATE_RECIPROCAL=true >/dev/null
    bin="./MDBench-CP-4xN-${TOOL_TAG}-${dtype}-geometric${suffix}"
    for half in 0 1; do
      # Half lists scatter j-forces from several threads, which is not thread safe
      OMP_NUM_THREADS=1 "${bin}" -gmxbench -size 1 -half "${half}" --kernel-bench 1 \
          --dump-forces "${WORK}/f.txt" >/dev/null
      if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/ref_${rule}.txt" "${WORK}/f.txt" \
           --tol "${tol}" --tol-max "$(python3 -c "print(${tol}*10)")" \
           --label "${rule} ${dtype} half=${half}"; then
        status=1
      fi
    done
  done
done

if [ "${status}" -eq 0 ]; then
  echo "Coulomb and exclusion validation PASSED"
else
  echo "Coulomb and exclusion validation FAILED" >&2
fi
exit "${status}"
