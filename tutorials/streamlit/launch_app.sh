#!/usr/bin/env bash
# Launch the OpenDiS Streamlit app.
# Headless: no browser window opens. Open the printed URL yourself, in a web browser.
set -euo pipefail
cd "$(dirname "$0")"

ADDR="${STREAMLIT_SERVER_ADDRESS:-localhost}"

# The app's env, opendis-vis (see README.md). The Python is, in this order: OPENDIS_PYTHON;
# the env at /opt/anaconda3; the env that conda knows by that name, wherever it is, so that a
# new terminal needs no activation; python3 on the PATH.
PY="${OPENDIS_PYTHON:-/opt/anaconda3/envs/opendis-vis/bin/python}"
if [[ ! -x "${PY}" ]] && command -v conda >/dev/null 2>&1; then
  PY="$(conda run -n opendis-vis python -c 'import sys; print(sys.executable)' 2>/dev/null \
    | awk 'NF { last = $0 } END { print last }' || true)"
fi
if [[ ! -x "${PY}" ]]; then
  PY="$(command -v python3 || command -v python)"
fi

# The port: STREAMLIT_SERVER_PORT if set; else 8505, or the next free port after it, so that the
# app runs beside other Streamlit apps, which start from Streamlit's default port, 8501.
free_port() {
  "${PY}" - "${ADDR}" <<'EOF'
import socket, sys
for port in range(8505, 8600):
    with socket.socket() as s:
        try:
            s.bind((sys.argv[1], port))
        except OSError:
            continue
    print(port)
    break
EOF
}
PORT="${STREAMLIT_SERVER_PORT:-$(free_port)}"
if [[ -z "${PORT}" ]]; then
  echo "No free port from 8505 to 8599: set STREAMLIT_SERVER_PORT." >&2
  exit 1
fi

echo "Starting OpenDiS (headless; no browser window opens)"
echo "  Python: ${PY}"
echo "  URL:    http://${ADDR}:${PORT}"
echo "  Open that URL in a web browser."
echo "  After editing the code, restart this script; reloading the tab is not enough."
echo

exec "${PY}" -m streamlit run app.py \
  --server.address "${ADDR}" \
  --server.port "${PORT}" \
  --server.headless true \
  --browser.gatherUsageStats false
