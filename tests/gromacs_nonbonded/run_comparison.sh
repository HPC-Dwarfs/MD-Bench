#!/usr/bin/env bash
set -euo pipefail

# Compares MD-Bench's 4xN cluster-pair kernel with GROMACS' NBNxM kernel on the
# GROMACS non-bonded benchmark (gmx nonbonded-benchmark), in forces and kernel time.
#
# Forces: both programs evaluate one kernel pass on the same system (3000 SPC/E atoms
# times -size, charges, intramolecular exclusions, cut-off 1.0, no pair-list buffer)
# and the per-atom forces are compared, so the physics, exclusions and the Coulomb
# treatment have to agree for the check to pass.
#
# Time: GROMACS' loop of kernel calls on one fixed pair list is mirrored by
# MD-Bench's --kernel-bench. Both use one thread and single precision.
#
# GROMACS' benchmark only prints timings, so it has to be built from a tree with
# tests/gromacs_nonbonded/gromacs-force-dump.patch applied, which adds a "-forces
# <file>" option (see the README in this directory).
#
# Usage:
#   GMX=/path/to/gmx ./tests/gromacs_nonbonded/run_comparison.sh [size ...]
#
# Environment:
#   GMX         gmx binary built with the force dump patch (required)
#   SIMD        MD-Bench SIMD (default AVX2); GROMACS' 4xM kernel is used
#   ITER        kernel iterations for the timing (default 500)
#   PIN_CORE    pin both programs to this core with taskset (optional)
#   RULES       Coulomb treatments to run (default "ewald rf")

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DIR="${ROOT_DIR}/tests/gromacs_nonbonded"

TOOLCHAIN="${TOOLCHAIN:-GCC}"
ISA="${ISA:-X86}"
SIMD="${SIMD:-AVX2}"
PYTHON="${PYTHON:-python3}"
ITER="${ITER:-500}"
RULES="${RULES:-ewald rf}"
SIZES=("$@")
[ ${#SIZES[@]} -eq 0 ] && SIZES=(1)

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

if [ "${SIMD}" = "NONE" ]; then
  TOOL_TAG="${TOOLCHAIN}-${ISA}"
else
  TOOL_TAG="${TOOLCHAIN}-${ISA}-${SIMD}"
fi

build() { # rule accurate
  make TOOLCHAIN="${TOOLCHAIN}" ISA="${ISA}" SIMD="${SIMD}" DATA_TYPE=SP \
       OPT_SCHEME=clusterpair CLUSTER_PAIR_KERNEL=4xN COULOMB_RULE="$1" \
       ACCURATE_RECIPROCAL="$2" >/dev/null
}
binary() { # rule accurate
  local suffix=""
  [ "$2" = "true" ] && suffix="-accrecip"
  echo "./MDBench-CP-4xN-${TOOL_TAG}-SP-geometric-$1${suffix}"
}
mdb_value() { # file label
  awk -F': *' -v k="$2" '$1 ~ k {gsub(/ us$/, "", $2); print $2}' "$1" | head -n 1
}

status=0
for size in "${SIZES[@]}"; do
  for rule in ${RULES}; do
    gmx_coulomb=ewald
    [ "${rule}" = "rf" ] && gmx_coulomb=reaction-field
    echo
    echo "=== size ${size} ($((3000 * size)) atoms), Coulomb: ${rule} ==="

    "${PIN[@]}" "${GMX}" nonbonded-benchmark -size "${size}" -nt 1 -simd 4xm \
        -coulomb "${gmx_coulomb}" -combrule geometric -cutoff 1.0 -iter "${ITER}" \
        -warmup 50 -time -forces "${WORK}/gmx_forces.txt" >"${WORK}/gmx.log" 2>&1
    gmx_line="$(grep -E '^(Ewald|RF) ' "${WORK}/gmx.log" | tail -n 1)"
    gmx_time="$(echo "${gmx_line}" | awk '{print $(NF-2)}')"
    gmx_total="$(echo "${gmx_line}" | awk '{print $(NF-1)}')"
    gmx_useful="$(echo "${gmx_line}" | awk '{print $NF}')"

    build "${rule}" true
    bin_acc="$(binary "${rule}" true)"
    for half in 1 0; do
      OMP_NUM_THREADS=1 "${PIN[@]}" "${bin_acc}" -gmxbench -size "${size}" -half "${half}" \
          --kernel-bench "${ITER}" --dump-forces "${WORK}/mdb_forces_${half}.txt" \
          >"${WORK}/mdb_acc_${half}.log"
      if ! "${PYTHON}" "${DIR}/compare_forces.py" "${WORK}/gmx_forces.txt" \
           "${WORK}/mdb_forces_${half}.txt" --tol 1e-5 --tol-max 1e-4 \
           --label "MD-Bench vs GROMACS, half=${half}"; then
        status=1
      fi
    done

    build "${rule}" false
    bin_fast="$(binary "${rule}" false)"
    OMP_NUM_THREADS=1 "${PIN[@]}" "${bin_fast}" -gmxbench -size "${size}" -half 1 \
        --kernel-bench "${ITER}" >"${WORK}/mdb_fast.log"

    printf '\n%-44s %12s %14s %14s\n' "kernel (1 thread, SP)" "us/iter" "pairs/us" "useful pairs/us"
    printf '%-44s %12s %14s %14s\n' "GROMACS 4xM, analytical corr." "${gmx_time}" "${gmx_total}" "${gmx_useful}"
    for variant in "acc_1:MD-Bench 4xN half, accurate 1/r^2" "acc_0:MD-Bench 4xN full, accurate 1/r^2" "fast:MD-Bench 4xN half, 1/r^2 estimate"; do
      f="${WORK}/mdb_${variant%%:*}.log"
      [ "${variant%%:*}" = "fast" ] && f="${WORK}/mdb_fast.log"
      printf '%-44s %12s %14s %14s\n' "${variant#*:}" \
        "$(mdb_value "${f}" 'Time per iteration')" \
        "$(mdb_value "${f}" 'Computed pairs/us')" \
        "$(mdb_value "${f}" 'Useful pairs/us')"
    done
    printf '\nGROMACS pair list: %s\n' "$(grep '^# pairs' "${WORK}/gmx_forces.txt")"
    printf 'MD-Bench pair list: %s\n' "$(grep 'Cluster pairs' "${WORK}/mdb_acc_1.log" | sed 's/^ *//')"
  done
done

if [ "${status}" -eq 0 ]; then
  echo
  echo "MD-Bench and GROMACS agree on the forces"
else
  echo "Force comparison FAILED" >&2
fi
exit "${status}"
