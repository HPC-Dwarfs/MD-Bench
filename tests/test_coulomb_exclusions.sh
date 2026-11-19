#!/usr/bin/env bash
set -euo pipefail

# Forces/energies of the CPU kernels (Coulomb, exclusions, -energy, -halflj) on -gmxbench vs float64 reference
#
# Usage:
#   ./tests/test_coulomb_exclusions.sh
#   SIMD=AVX2 ./tests/test_coulomb_exclusions.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
SIMD="${SIMD:-AVX2}"
PYTHON="${PYTHON:-python3}"

if ! "${PYTHON}" -c "import numpy, scipy" >/dev/null 2>&1; then
  echo "SKIP: ${PYTHON} with numpy and scipy is required"
  exit 0
fi

cd "${ROOT_DIR}"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/mdbench_coulomb.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

status=0
# kernel:data type of every kernel the SIMD back-end can run (REF: cluster-pair
# reference kernel, VL: verlet lists)
CONFIGS=("4xN:DP" "REF:DP" "REF:SP" "VL:DP" "VL:SP" "VLSOA:SP")
case "${SIMD}" in
  AVX512) ;;
  *) CONFIGS+=("4xN:SP") ;;
esac
case "${SIMD}" in
  AVX2) CONFIGS+=("2xNN:SP" "2xN:DP") ;;
  AVX512) CONFIGS+=("2xNN:SP") ;;
esac
for rule in rf ewald ewald_table none; do
  for size in 1 2; do
    "${PYTHON}" "${DIR}/nbnxm_reference.py" --size "${size}" --coulomb "${rule}" \
        --out "${WORK}/ref_${rule}_${size}.txt"
  done
  for config in "${CONFIGS[@]}"; do
    IFS=: read -r kernel dtype <<<"${config}"
    # Energies are sums over all pairs, their single precision error is larger
    if [ "${dtype}" = "DP" ]; then tol="1e-7"; etol="1e-9"; else tol="5e-5"; etol="2e-5"; fi
    # GROMACS sizes the Ewald table for single precision force accuracy, also in double
    if [ "${rule}" = "ewald_table" ] && [ "${dtype}" = "DP" ]; then tol="1e-6"; etol="1e-8"; fi
    cfg=(TOOLCHAIN="${TOOLCHAIN}" SIMD="${SIMD}" DATA_TYPE="${dtype}"
         COULOMB_RULE="${rule}")
    kernel_opts=()
    case "${kernel}" in
      VL) cfg+=(OPT_SCHEME=verletlist) ;;
      # The neighbor list layout is not part of the binary name: own directory, no rebuilds
      VLSOA) cfg+=(OPT_SCHEME=verletlist NBLIST_DATA_LAYOUT=SOA PREFIX_DIR="${WORK}/vlsoa") ;;
      REF) cfg+=(OPT_SCHEME=clusterpair USE_REFERENCE_KERNEL=true) ;;
      *)
        cfg+=(OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL="${kernel}")
        kernel_opts=(--accurate-recip)
        ;;
    esac
    make "${cfg[@]}" >/dev/null
    bin="$(make -s "${cfg[@]}" print-target)"
    # Each kernel option set with the geometric rule, the plain forces also with full
    OPTIONS=("--lj-comb geometric" "--lj-comb full" "--lj-comb geometric -energy")
    case "${kernel}" in
      VL*) ;;
      *) OPTIONS+=("--lj-comb geometric -halflj -energy") ;;
    esac
    for size in 1 2; do
      for half in 0 1; do
        for opts in "${OPTIONS[@]}"; do
          energy_args=()
          if [[ " ${opts} " == *" -energy "* ]]; then
            energy_args=(--require-energy)
          fi
          # shellcheck disable=SC2086
          OMP_NUM_THREADS=4 "${bin}" -gmxbench -size "${size}" -half "${half}" \
              --kernel-bench 1 ${opts} "${kernel_opts[@]}" --dump-forces "${WORK}/f.txt" >/dev/null
          if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/ref_${rule}_${size}.txt" \
               "${WORK}/f.txt" --tol "${tol}" --tol-max "$(python3 -c "print(${tol}*10)")" \
               --etol "${etol}" "${energy_args[@]}" \
               --label "${kernel} ${rule} ${dtype} size=${size} half=${half} ${opts}"; then
            status=1
          fi
        done
      done
    done
  done
done

if [ "${status}" -eq 0 ]; then
  echo "Coulomb and exclusion validation PASSED"
else
  echo "Coulomb and exclusion validation FAILED" >&2
fi
exit "${status}"
