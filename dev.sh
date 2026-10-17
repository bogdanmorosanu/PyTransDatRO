#!/usr/bin/env bash
#
# PyTransDatRO development runner (macOS / Linux).
#
# Starts the FastAPI backend (uvicorn, port 8000, --reload) and the Vite
# frontend dev server (port 5173, proxies /api,/health,/docs to the backend),
# opens the web app in the default browser, and shuts both down cleanly on
# Ctrl+C. This is the POSIX counterpart to dev.ps1 (Windows/PowerShell).
#
# Usage:
#   ./dev.sh                 # start both dev servers + open browser
#   NO_BROWSER=1 ./dev.sh    # start without opening the browser
#   API_PORT=8001 ./dev.sh   # override the API port
#
# For a single-process unified run instead (serves the pre-built webapp/dist
# from the API, no hot reload), build the frontend then run the API alone:
#   npm --prefix webapp run build
#   .venv/bin/python -m uvicorn api.main:app --port 8000

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

API_PORT="${API_PORT:-8000}"

# Prefer the project venv interpreter; fall back to system python3.
if [[ -x "$ROOT_DIR/.venv/bin/python" ]]; then
  PYTHON="$ROOT_DIR/.venv/bin/python"
else
  PYTHON="python3"
  echo "warning: .venv not found, falling back to system python3 (needs Python >=3.10 with the [api] extras installed)" >&2
fi

echo "================================================="
echo "        PyTransDatRO Development Runner           "
echo "================================================="

API_PID=""
WEB_PID=""

cleanup() {
  echo ""
  echo "Stopping development servers..."
  # Kill child processes first (uvicorn's reload worker, npm's node child), then the parents.
  [[ -n "$API_PID" ]] && { pkill -P "$API_PID" 2>/dev/null || true; kill "$API_PID" 2>/dev/null || true; }
  [[ -n "$WEB_PID" ]] && { pkill -P "$WEB_PID" 2>/dev/null || true; kill "$WEB_PID" 2>/dev/null || true; }
  wait 2>/dev/null || true
  echo "All servers stopped. Goodbye!"
}
trap cleanup INT TERM EXIT

echo "[1/2] Starting REST API on http://127.0.0.1:${API_PORT} ..."
"$PYTHON" -m uvicorn api.main:app --host 127.0.0.1 --port "$API_PORT" --reload &
API_PID=$!

echo "[2/2] Starting Web App (Vite) on http://localhost:5173 ..."
npm --prefix webapp run dev &
WEB_PID=$!

# Give Vite a moment to bind before opening the browser.
sleep 2

if [[ -z "${NO_BROWSER:-}" ]]; then
  echo "Opening web app in default browser..."
  if command -v open >/dev/null 2>&1; then
    open "http://localhost:5173" || true          # macOS
  elif command -v xdg-open >/dev/null 2>&1; then
    xdg-open "http://localhost:5173" || true       # Linux
  fi
fi

echo ""
echo "Development servers are active:"
echo "  - REST API:     http://127.0.0.1:${API_PORT}"
echo "  - Swagger Docs: http://127.0.0.1:${API_PORT}/docs"
echo "  - Web App:      http://localhost:5173"
echo ""
echo "Press [Ctrl+C] to stop all servers and exit."

# Block until either server exits (then the EXIT trap stops the other).
# Portable poll loop — avoids `wait -n`, which needs bash >=4.3 (macOS ships 3.2).
while kill -0 "$API_PID" 2>/dev/null && kill -0 "$WEB_PID" 2>/dev/null; do
  sleep 1
done
