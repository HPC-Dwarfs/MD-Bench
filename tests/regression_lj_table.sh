#!/usr/bin/env bash
set -euo pipefail

# Equivalence regression: tabulated/spline LJ (-f lj_table) must reproduce the
# analytic LJ kernel (-f lj). The two share identical units, so over a modest
# number of timesteps the thermo trajectory must agree to within the table's
# interpolation accuracy (default 1000 knots).
#
# Usage:
#   ./tests/regression_lj_table.sh /path/to/MDBench-<TAG>
# or set MDBENCH_BIN in the environment.

BIN="${1:-${MDBENCH_BIN:-}}"

if [[ -z "${BIN}" ]]; then
    echo "Usage: $0 /path/to/MDBench-<TAG>  (or set MDBENCH_BIN)" >&2
    exit 1
fi

if [[ ! -x "${BIN}" ]]; then
    echo "Binary '${BIN}' is not executable" >&2
    exit 1
fi

NSTEPS="${NSTEPS:-100}"
RELTOL="${RELTOL:-1e-3}"

# Skip gracefully on builds that do not implement the tabulated kernel yet
# (e.g. the clusterpair scheme): the binary exits with "Unknown force field".
probe="$("${BIN}" -f lj_table -n 0 -nx 4 -ny 4 -nz 4 2>&1 || true)"
if echo "${probe}" | grep -qi "Unknown force field"; then
    echo "SKIP: '${BIN}' does not support -f lj_table (tabulated kernel not built)"
    exit 0
fi

LJ_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_lj.XXXXXX")"
TAB_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_lj_table.XXXXXX")"

echo "Running analytic LJ:  ${BIN} -f lj -n ${NSTEPS}"
"${BIN}" -f lj       -n "${NSTEPS}" -nx 12 -ny 12 -nz 12 >"${LJ_LOG}"
echo "Running tabulated LJ: ${BIN} -f lj_table -n ${NSTEPS}"
"${BIN}" -f lj_table -n "${NSTEPS}" -nx 12 -ny 12 -nz 12 >"${TAB_LOG}"

parse_last() {
    grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "$1" | tail -n 1 || true
}

lj_line="$(parse_last "${LJ_LOG}")"
tab_line="$(parse_last "${TAB_LOG}")"

if [[ -z "${lj_line}" || -z "${tab_line}" ]]; then
    echo "Could not find thermo lines in output." >&2
    echo "Analytic log:  ${LJ_LOG}" >&2
    echo "Tabulated log: ${TAB_LOG}" >&2
    exit 1
fi

T_lj=$(echo "${lj_line}"  | awk '{print $2}')
P_lj=$(echo "${lj_line}"  | awk '{print $3}')
T_tab=$(echo "${tab_line}" | awk '{print $2}')
P_tab=$(echo "${tab_line}" | awk '{print $3}')

echo "Analytic : T=${T_lj} P=${P_lj}"
echo "Tabulated: T=${T_tab} P=${P_tab}"

PYTHON="$(command -v python3 || command -v python)"
"${PYTHON}" - "$T_lj" "$T_tab" "$P_lj" "$P_tab" "$RELTOL" << 'PY'
import sys
T0, T1, P0, P1, tol = (float(x) for x in sys.argv[1:6])

def rel(a, b):
    return abs(a - b) / (abs(a) + 1e-300)

rt, rp = rel(T0, T1), rel(P0, P1)
print(f"relative diff: temperature={rt:.3e} pressure={rp:.3e} (tol={tol:.1e})")
if rt > tol or rp > tol:
    print("FAIL: tabulated LJ diverged from analytic beyond tolerance")
    sys.exit(1)
print("PASS: tabulated LJ matches analytic LJ")
PY

echo "Tabulated LJ regression completed successfully."
