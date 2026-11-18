#!/usr/bin/env bash
set -euo pipefail

# Forces/energies of -gmxbench with MPI (1, 2, 4 ranks) vs the float64 reference
#
# Usage:
#   ./tests/test_mpi_gmxbench.sh
#   SIMD=AVX2 ./tests/test_mpi_gmxbench.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
SIMD="${SIMD:-AVX2}"
PYTHON="${PYTHON:-python3}"
SIZE="${SIZE:-8}"

if ! command -v mpicc >/dev/null 2>&1 || ! command -v mpirun >/dev/null 2>&1; then
  echo "SKIP: mpicc/mpirun not found"
  exit 0
fi
if ! "${PYTHON}" -c "import numpy, scipy" >/dev/null 2>&1; then
  echo "SKIP: ${PYTHON} with numpy and scipy is required"
  exit 0
fi

MPI_EXTRA_FLAGS=()
# Open MPI refuses more processes than slots (one task of a batch job) without it
if mpirun --help 2>&1 | grep -q -- "--oversubscribe" || mpirun --version 2>&1 | grep -q "Open MPI"; then
  MPI_EXTRA_FLAGS=(--oversubscribe)
fi

cd "${ROOT_DIR}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mdbench_mpi_gmx.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

# The SIMD kernel of the back-end in single precision
case "${SIMD}" in
  AVX512) SIMD_KERNEL=2xNN ;;
  *) SIMD_KERNEL=4xN ;;
esac

"${PYTHON}" "${DIR}/nbnxm_reference.py" --size "${SIZE}" --coulomb ewald --out "${WORK}/ref.txt"

status=0
for kernel in "${SIMD_KERNEL}" REF VL; do
  kernel_opts=()
  case "${kernel}" in
    VL) args=(OPT_SCHEME=verletlist) ;;
    REF) args=(OPT_SCHEME=clusterpair USE_REFERENCE_KERNEL=true) ;;
    *)
      args=(OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL="${kernel}")
      kernel_opts=(--accurate-recip)
      ;;
  esac
  cfg=(TOOLCHAIN="${TOOLCHAIN}" SIMD="${SIMD}" DATA_TYPE=SP ENABLE_MPI=true
       COULOMB_RULE=ewald PREFIX_DIR="${WORK}/${kernel}" "${args[@]}")
  make "${cfg[@]}" >/dev/null
  bin="$(make -s "${cfg[@]}" print-target)"
  for np in 1 2 4; do
    for half in 0 1; do
      OMP_NUM_THREADS=1 mpirun "${MPI_EXTRA_FLAGS[@]}" -np "${np}" "${bin}" -gmxbench \
          -size "${SIZE}" -half "${half}" --kernel-bench 1 -energy "${kernel_opts[@]}" \
          --dump-forces "${WORK}/f.txt" >"${WORK}/run.log" 2>&1 || {
        echo "FAILED: ${kernel} np=${np} half=${half}"
        tail -n 5 "${WORK}/run.log"
        status=1
        continue
      }
      if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/ref.txt" "${WORK}/f.txt" \
           --tol 5e-5 --tol-max 5e-4 --etol 2e-5 \
           --label "${kernel} np=${np} half=${half}"; then
        status=1
      fi
    done
  done
done

if [ "${status}" -eq 0 ]; then
  echo "MPI GROMACS benchmark validation PASSED"
else
  echo "MPI GROMACS benchmark validation FAILED" >&2
fi
exit "${status}"
