#!/usr/bin/env bash
set -euo pipefail

# Validates the electrostatics (reaction-field, Ewald) and molecular exclusions of the
# cluster-pair SIMD kernels (4xN, 2xNN, 2xN as far as the SIMD back-end can run them), the
# cluster-pair reference kernel (REF) and the verlet-list kernels (VL) against an exact
# float64 reference.
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
# kernel:data type of every kernel the SIMD back-end can run (REF: cluster-pair
# reference kernel, VL: verlet lists)
CONFIGS=("4xN:DP" "4xN:SP" "REF:DP" "REF:SP" "VL:DP" "VL:SP")
case "${SIMD}" in
  AVX2) CONFIGS+=("2xNN:SP" "2xN:DP") ;;
  AVX512) CONFIGS+=("2xNN:SP") ;;
esac
for rule in rf ewald none; do
  "${PYTHON}" "${DIR}/nbnxm_reference.py" --size 1 --coulomb "${rule}" --out "${WORK}/ref_${rule}.txt"
  for config in "${CONFIGS[@]}"; do
    IFS=: read -r kernel dtype <<<"${config}"
    if [ "${dtype}" = "DP" ]; then tol="1e-7"; else tol="5e-5"; fi
    excl=auto
    suffix="-${rule}"
    if [ "${rule}" = "none" ]; then excl=true; suffix=""; fi
    case "${kernel}" in
      VL)
        make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" DATA_TYPE="${dtype}" \
             OPT_SCHEME=verletlist COULOMB_RULE="${rule}" EXCLUSIONS="${excl}" >/dev/null
        bin="./MDBench-VL-${TOOL_TAG}-${dtype}-geometric${suffix}"
        ;;
      REF)
        make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" DATA_TYPE="${dtype}" \
             OPT_SCHEME=clusterpair USE_REFERENCE_KERNEL=true COULOMB_RULE="${rule}" \
             EXCLUSIONS="${excl}" >/dev/null
        bin="./MDBench-CP-auto-${TOOL_TAG}-${dtype}-geometric${suffix}"
        ;;
      *)
        make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" DATA_TYPE="${dtype}" \
             OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL="${kernel}" COULOMB_RULE="${rule}" \
             EXCLUSIONS="${excl}" ACCURATE_RECIPROCAL=true >/dev/null
        bin="./MDBench-CP-${kernel}-${TOOL_TAG}-${dtype}-geometric${suffix}-accrecip"
        ;;
    esac
    if [ "${rule}" = "none" ]; then bin="${bin}-excl"; fi
    for half in 0 1; do
      OMP_NUM_THREADS=4 "${bin}" -gmxbench -size 1 -half "${half}" --kernel-bench 1 \
          --dump-forces "${WORK}/f.txt" >/dev/null
      if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/ref_${rule}.txt" "${WORK}/f.txt" \
           --tol "${tol}" --tol-max "$(python3 -c "print(${tol}*10)")" \
           --label "${kernel} ${rule} ${dtype} half=${half}"; then
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
