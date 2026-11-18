#!/usr/bin/env bash
set -euo pipefail

# Forces/energies of the GPU kernels (Coulomb, exclusions, -energy) on -gmxbench vs float64 reference
#
# Usage:
#   ./tests/test_gpu_coulomb.sh
#   CUDA_ARCH=sm_90 ./tests/test_gpu_coulomb.sh
#   TOOLCHAIN=HIPCC GPU_ARCH=gfx942 ./tests/test_gpu_coulomb.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

TOOLCHAIN="${TOOLCHAIN:-NVCC}"
CUDA_ARCH="${CUDA_ARCH:-sm_80}"
PYTHON="${PYTHON:-python3}"

if [ "${TOOLCHAIN}" = "HIPCC" ]; then
  if ! command -v hipcc >/dev/null 2>&1; then
    echo "SKIP: hipcc not found"
    exit 0
  fi
  if [ -z "${GPU_ARCH:-}" ]; then
    echo "SKIP: set GPU_ARCH for HIP (e.g. gfx942)"
    exit 0
  fi
  ARCH_ARG="GPU_ARCH=${GPU_ARCH}"
else
  if ! command -v nvcc >/dev/null 2>&1; then
    echo "SKIP: nvcc not found"
    exit 0
  fi
  if ! nvidia-smi -L >/dev/null 2>&1; then
    echo "SKIP: no CUDA GPU found"
    exit 0
  fi
  ARCH_ARG="CUDA_ARCH=${CUDA_ARCH}"
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
  for size in 1 8; do
    "${PYTHON}" "${DIR}/nbnxm_reference.py" --size "${size}" --coulomb "${rule}" \
        --out "${WORK}/ref_${rule}_${size}.txt"
  done
  for kernel in gpusimple supercluster VL; do
    for dtype in SP DP; do
      if [ "${dtype}" = "DP" ]; then tol="1e-7"; etol="1e-9"; else tol="5e-5"; etol="2e-5"; fi
      cfg=(TOOLCHAIN="${TOOLCHAIN}" "${ARCH_ARG}" DATA_TYPE="${dtype}"
           COULOMB_RULE="${rule}")
      if [ "${kernel}" = "VL" ]; then
        cfg+=(OPT_SCHEME=verletlist)
      else
        cfg+=(OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL="${kernel}")
      fi
      make "${cfg[@]}" >/dev/null
      bin="$(make -s "${cfg[@]}" print-target)"
      for size in 1 8; do
        for half in 0 1; do
          for opts in "-n 0" "--kernel-bench 2" "--kernel-bench 2 -energy"; do
            # shellcheck disable=SC2086
            "${bin}" -gmxbench -size "${size}" -half "${half}" ${opts} \
                --dump-forces "${WORK}/f.txt" >/dev/null
            if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/ref_${rule}_${size}.txt" \
                 "${WORK}/f.txt" --tol "${tol}" --tol-max "$(python3 -c "print(${tol}*10)")" \
                 --etol "${etol}" \
                 --label "${kernel} ${rule} ${dtype} size=${size} half=${half} ${opts}"; then
              status=1
            fi
          done
        done
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
