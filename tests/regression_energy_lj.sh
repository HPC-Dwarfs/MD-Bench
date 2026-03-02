#!/usr/bin/env bash
set -euo pipefail

# Energy/temperature stability regression for default LJ solid copper testcase.
#
# Usage:
#   ./tests/regression_energy_lj.sh /path/to/MDBench-<TAG>
#

BIN="${1:-${MDBENCH_BIN:-}}"

if [[ -z "${BIN}" ]]; then
    echo "Usage: $0 /path/to/MDBench-<TAG>  (or set MDBENCH_BIN)" >&2
    exit 1
fi

if [[ ! -x "${BIN}" ]]; then
    echo "Binary '${BIN}' is not executable" >&2
    exit 1
fi

OUT_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_lj_energy.XXXXXX")"

echo "Running LJ energy regression with binary: ${BIN}"
# Default LJ solid testcase, extend timesteps for a more meaningful drift check.
"${BIN}" -n 1000 >"${OUT_LOG}"

echo "Simulation finished, parsing thermo output..."

thermo_lines="$(grep -E '^[0-9]+\s+[0-9.eE+-]+' "${OUT_LOG}" || true)"
if [[ -z "${thermo_lines}" ]]; then
    echo "Could not find thermo lines in output; check that MD-Bench printed stats." >&2
    echo "Full log at: ${OUT_LOG}" >&2
    exit 1
fi

first_line="$(echo "${thermo_lines}" | head -n 1)"
last_line="$(echo "${thermo_lines}" | tail -n 1)"

T0=$(echo "${first_line}" | awk '{print $2}')
T1=$(echo "${last_line}"  | awk '{print $2}')

drift=$(python - "$T0" "$T1" << 'PY'
import sys, math
T0=float(sys.argv[1]); T1=float(sys.argv[2])
if T0 == 0.0:
    print(0.0)
else:
    print(abs(T1-T0)/T0)
PY
)

echo "Initial T=${T0}, final T=${T1}, relative drift=${drift}"
echo "NOTE: This test currently reports drift but does not enforce a hard threshold."
echo "LJ energy/temperature regression completed (drift reported above)."

