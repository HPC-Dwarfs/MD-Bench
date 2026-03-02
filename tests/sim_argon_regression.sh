#!/usr/bin/env bash
set -euo pipefail

# Simple correctness-style regression test for the argon testcase.
#
# Usage:
#   ./tests/sim_argon_regression.sh /path/to/MDBench-<TAG>
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

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT_DIR}/data/argon"

OUT_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_argon.XXXXXX")"

echo "Running argon regression with binary: ${BIN}"
"${BIN}" -i "${DATA_DIR}/input.gro" -p "${DATA_DIR}/mdbench_params.conf" -n 500 >"${OUT_LOG}"

echo "Simulation finished, parsing thermo output..."

# Grab last thermo line (temperature and pressure are printed after 'step').
last_tp_line="$(grep -E '^[0-9]+\s+[0-9.eE+-]+' "${OUT_LOG}" | tail -n 1 || true)"

if [[ -z "${last_tp_line}" ]]; then
    echo "Could not find thermo line in output; check that MD-Bench printed stats." >&2
    echo "Full log at: ${OUT_LOG}" >&2
    exit 1
fi

step=$(echo "${last_tp_line}" | awk '{print $1}')
temp=$(echo "${last_tp_line}" | awk '{print $2}')
press=$(echo "${last_tp_line}" | awk '{print $3}')

echo "Last thermo: step=${step} temp=${temp} pressure=${press}"
echo "Argon regression completed successfully."

