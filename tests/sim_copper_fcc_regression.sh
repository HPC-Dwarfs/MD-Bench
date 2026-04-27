#!/usr/bin/env bash
set -euo pipefail

# Regression test for the default copper FCC lattice testcase (LJ, no input file).
#
# Usage:
#   ./tests/sim_copper_fcc_regression.sh /path/to/MDBench-<TAG>
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

# Reference values produced by a clean run at step 200 (ICX, AVX512, DP).
REF_TEMP="7.961676e-01"
REF_PRESS="6.721195e-01"
TEMP_TOL="1e-4"
PRESS_TOL="1e-3"

OUT_LOG="$(mktemp "${TMPDIR:-/tmp}/mdbench_copper_fcc.XXXXXX")"

echo "Running copper FCC regression with binary: ${BIN}"
"${BIN}" -n 200 >"${OUT_LOG}"

echo "Simulation finished, parsing thermo output..."

last_tp_line="$(grep -E '^[[:space:]]*[0-9]+[[:space:]]+[0-9.eE+-]+' "${OUT_LOG}" | tail -n 1 || true)"

if [[ -z "${last_tp_line}" ]]; then
    echo "Could not find thermo line in output; check that MD-Bench printed stats." >&2
    echo "Full log at: ${OUT_LOG}" >&2
    exit 1
fi

step=$(echo "${last_tp_line}" | awk '{print $1}')
temp=$(echo "${last_tp_line}" | awk '{print $2}')
press=$(echo "${last_tp_line}" | awk '{print $3}')

echo "Last thermo: step=${step} temp=${temp} pressure=${press}"

pass=1
temp_diff=$(awk -v a="${temp}" -v b="${REF_TEMP}" 'BEGIN { d=(a-b); if (d<0) d=-d; print d/b }')
press_diff=$(awk -v a="${press}" -v b="${REF_PRESS}" 'BEGIN { d=(a-b); if (d<0) d=-d; print d/b }')

temp_ok=$(awk -v d="${temp_diff}" -v t="${TEMP_TOL}" 'BEGIN { print (d <= t) ? "1" : "0" }')
press_ok=$(awk -v d="${press_diff}" -v t="${PRESS_TOL}" 'BEGIN { print (d <= t) ? "1" : "0" }')

if [[ "${temp_ok}" != "1" ]]; then
    echo "FAIL: temperature ${temp} deviates from reference ${REF_TEMP} by ${temp_diff} (tolerance ${TEMP_TOL})" >&2
    pass=0
fi

if [[ "${press_ok}" != "1" ]]; then
    echo "FAIL: pressure ${press} deviates from reference ${REF_PRESS} by ${press_diff} (tolerance ${PRESS_TOL})" >&2
    pass=0
fi

if [[ "${pass}" != "1" ]]; then
    echo "Full log at: ${OUT_LOG}" >&2
    exit 1
fi

rm -f "${OUT_LOG}"
echo "Copper FCC regression completed successfully."
