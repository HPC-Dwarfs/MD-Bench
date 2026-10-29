#!/usr/bin/env bash
set -euo pipefail

# Validates the electrostatics (reaction-field, Ewald) and molecular exclusions of the
# cluster-pair GPU kernels (gpusimple and supercluster) against an exact float64
# reference, like test_coulomb_exclusions.sh does for the CPU SIMD kernels.
#
# The system is the water box of the GROMACS non-bonded benchmark (-gmxbench), the
# reference is tests/gromacs_nonbonded/nbnxm_reference.py. Full and half neighbor lists
# are checked, in single and double precision.
#
# Needs nvcc, a CUDA GPU and python3 with numpy and scipy, the test is skipped otherwise.
#
# Usage:
#   ./tests/test_gpu_coulomb.sh
#   CUDA_ARCH=sm_90 ./tests/test_gpu_coulomb.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

ISA="${ISA:-X86}"
CUDA_ARCH="${CUDA_ARCH:-sm_80}"
PYTHON="${PYTHON:-python3}"

if ! command -v nvcc >/dev/null 2>&1; then
  echo "SKIP: nvcc not found"
  exit 0
fi
if ! nvidia-smi -L >/dev/null 2>&1; then
  echo "SKIP: no CUDA GPU found"
  exit 0
fi
if ! "${PYTHON}" -c "import numpy, scipy" >/dev/null 2>&1; then
  echo "SKIP: ${PYTHON} with numpy and scipy is required"
  exit 0
fi

cd "${ROOT_DIR}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mdbench_gpu_coulomb.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

status=0
for rule in rf ewald none; do
  "${PYTHON}" "${DIR}/nbnxm_reference.py" --size 1 --coulomb "${rule}" --out "${WORK}/ref_${rule}.txt"
  for kernel in gpusimple supercluster; do
    for dtype in SP DP; do
      if [ "${dtype}" = "DP" ]; then tol="1e-7"; else tol="5e-5"; fi
      excl=auto
      suffix="-${rule}"
      if [ "${rule}" = "none" ]; then excl=true; suffix="-excl"; fi
      make TOOLCHAIN=NVCC ISA="${ISA}" CUDA_ARCH="${CUDA_ARCH}" DATA_TYPE="${dtype}" \
           OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL="${kernel}" COULOMB_RULE="${rule}" \
           EXCLUSIONS="${excl}" >/dev/null
      bin="./MDBench-CP-${kernel}-NVCC-${ISA}-${dtype}-geometric${suffix}"
      for half in 0 1; do
        "${bin}" -gmxbench -size 1 -half "${half}" -n 0 \
            --dump-forces "${WORK}/f.txt" >/dev/null
        if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/ref_${rule}.txt" "${WORK}/f.txt" \
             --tol "${tol}" --tol-max "$(python3 -c "print(${tol}*10)")" \
             --label "${kernel} ${rule} ${dtype} half=${half}"; then
          status=1
        fi
      done
    done
  done
done

if [ "${status}" -eq 0 ]; then
  echo "GPU Coulomb and exclusion validation PASSED"
else
  echo "GPU Coulomb and exclusion validation FAILED" >&2
fi
exit "${status}"
