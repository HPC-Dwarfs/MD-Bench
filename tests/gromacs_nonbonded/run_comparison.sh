#!/usr/bin/env bash
set -euo pipefail

# Forces, energies and kernel time of MD-Bench 4xN/2xNN vs GROMACS 4xM/2xMM (patched gmx)
#
# Usage:
#   GMX=/path/to/gmx ./tests/gromacs_nonbonded/run_comparison.sh [size ...]
#
# Environment:
#   GMX         gmx built with gromacs-force-dump.patch (AVX2_256 for 4xM, AVX_512 for 2xMM)
#   SIMD        MD-Bench SIMD (default AVX2), AVX512 compares 2xNN with 2xMM
#   ITER        kernel iterations for the timing (default 500)
#   PIN_CORE    pin both programs to this core with taskset (optional)
#   RULES       Coulomb treatments to run (default "ewald ewald_table rf")
#   VARIANTS    ';'-separated kernel options (default ";-energy;-halflj")

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
SIMD="${SIMD:-AVX2}"
PYTHON="${PYTHON:-python3}"
ITER="${ITER:-500}"
RULES="${RULES:-ewald ewald_table rf}"
VARIANTS="${VARIANTS:-;-energy;-halflj}"
SIZES=("$@")
[ ${#SIZES[@]} -eq 0 ] && SIZES=(1)

if [ "${SIMD}" = "AVX512" ]; then
  KERNEL=2xNN
  GMX_KERNEL=2xmm
else
  KERNEL=4xN
  GMX_KERNEL=4xm
fi

if [ -z "${GMX:-}" ] || [ ! -x "${GMX}" ]; then
  echo "Set GMX to a gmx binary built with gromacs-force-dump.patch" >&2
  exit 2
fi
GMX_HELP="$("${GMX}" nonbonded-benchmark -h 2>&1 || true)"
if ! grep -q -- "-forces" <<<"${GMX_HELP}"; then
  echo "${GMX} has no 'nonbonded-benchmark -forces' option, apply gromacs-force-dump.patch" >&2
  exit 2
fi

PIN=()
[ -n "${PIN_CORE:-}" ] && PIN=(taskset -c "${PIN_CORE}")

WORK="$(mktemp -d "${TMPDIR:-/tmp}/mdbench_gmxcmp.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT
cd "${ROOT_DIR}"

build() { # rule, prints the binary
  local cfg=(TOOLCHAIN="${TOOLCHAIN}" SIMD="${SIMD}" DATA_TYPE=SP
             OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL="${KERNEL}" COULOMB_RULE="$1")
  make "${cfg[@]}" >/dev/null
  make -s "${cfg[@]}" print-target
}
mdb_value() { # file label
  awk -F': *' -v k="$2" '$1 ~ k {gsub(/ us$/, "", $2); print $2}' "$1" | head -n 1
}
row() { # label time total useful
  printf '%-48s %12s %14s %14s\n' "$1" "$2" "$3" "$4"
}

status=0
IFS=';' read -r -a VARIANT_LIST <<<"${VARIANTS}"
for size in "${SIZES[@]}"; do
  for rule in ${RULES}; do
    gmx_coulomb=(-coulomb ewald)
    [ "${rule}" = "rf" ] && gmx_coulomb=(-coulomb reaction-field)
    [ "${rule}" = "ewald_table" ] && gmx_coulomb=(-coulomb ewald -table)
    echo
    echo "=== size ${size} ($((3000 * size)) atoms), Coulomb: ${rule}, ${KERNEL} vs ${GMX_KERNEL} ==="

    bin="$(build "${rule}")"

    rows=()
    for variant in "${VARIANT_LIST[@]}"; do
      read -r -a opts <<<"${variant}"
      energy_args=()
      if [[ " ${variant} " == *" -energy "* ]]; then
        energy_args=(--require-energy)
      fi
      name="${variant:-forces}"
      "${PIN[@]}" "${GMX}" nonbonded-benchmark -size "${size}" -nt 1 -simd "${GMX_KERNEL}" \
          "${gmx_coulomb[@]}" "${opts[@]}" -combrule geometric -cutoff 1.0 -iter "${ITER}" \
          -warmup 50 -time -forces "${WORK}/gmx_forces.txt" >"${WORK}/gmx.log" 2>&1
      gmx_line="$(grep -E '^(Ewald|RF) ' "${WORK}/gmx.log" | tail -n 1)"
      rows+=("$(row "GROMACS ${GMX_KERNEL} ${name}" \
          "$(echo "${gmx_line}" | awk '{print $(NF-2)}')" \
          "$(echo "${gmx_line}" | awk '{print $(NF-1)}')" \
          "$(echo "${gmx_line}" | awk '{print $NF}')")")

      for half in 1 0; do
        # The full list only for the plain force kernel
        [ "${half}" = "0" ] && [ -n "${variant}" ] && continue
        OMP_NUM_THREADS=1 "${PIN[@]}" "${bin}" -gmxbench -size "${size}" -half "${half}" \
            --kernel-bench "${ITER}" --accurate-recip "${opts[@]}" \
            --dump-forces "${WORK}/mdb_forces.txt" \
            >"${WORK}/mdb.log"
        # GROMACS sums energies in single precision (2e-4 off the exact RF energy at size 8)
        if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/gmx_forces.txt" \
             "${WORK}/mdb_forces.txt" --tol 1e-5 --tol-max 1e-4 --etol 1e-3 \
             "${energy_args[@]}" --label "MD-Bench vs GROMACS ${name}, half=${half}"; then
          status=1
        fi
        list="half"
        [ "${half}" = "0" ] && list="full"
        rows+=("$(row "MD-Bench ${KERNEL} ${list} ${name}, accurate 1/r^2" \
            "$(mdb_value "${WORK}/mdb.log" 'Time per iteration')" \
            "$(mdb_value "${WORK}/mdb.log" 'Computed pairs/us')" \
            "$(mdb_value "${WORK}/mdb.log" 'Useful pairs/us')")")
      done

      OMP_NUM_THREADS=1 "${PIN[@]}" "${bin}" -gmxbench -size "${size}" -half 1 \
          --kernel-bench "${ITER}" "${opts[@]}" >"${WORK}/mdb_fast.log"
      rows+=("$(row "MD-Bench ${KERNEL} half ${name}, 1/r^2 estimate" \
          "$(mdb_value "${WORK}/mdb_fast.log" 'Time per iteration')" \
          "$(mdb_value "${WORK}/mdb_fast.log" 'Computed pairs/us')" \
          "$(mdb_value "${WORK}/mdb_fast.log" 'Useful pairs/us')")")
    done

    echo
    row "kernel (1 thread, SP)" "us/iter" "pairs/us" "useful pairs/us"
    printf '%s\n' "${rows[@]}"
    printf '\nGROMACS pair list: %s\n' "$(grep '^# pairs' "${WORK}/gmx_forces.txt")"
    printf 'MD-Bench pair list: %s\n' "$(grep 'Cluster pairs' "${WORK}/mdb.log" | sed 's/^ *//')"
  done
done

if [ "${status}" -eq 0 ]; then
  echo
  echo "MD-Bench and GROMACS agree on the forces and energies"
else
  echo "Force or energy comparison FAILED" >&2
fi
exit "${status}"
