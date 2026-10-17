#!/usr/bin/env bash
#
# PyTransDatRO E2E Test Runner (macOS / Linux).
#
# Usage:
#   ./test-e2e.sh       # run headless Playwright suite
#   ./test-e2e.sh --ui  # open interactive visual runner
#

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

if [[ -x "$ROOT_DIR/.venv/bin/python" ]]; then
  export PYTHON="$ROOT_DIR/.venv/bin/python"
fi

UI_FLAG=""
if [[ "${1:-}" == "--ui" || "${1:-}" == "-ui" ]]; then
  UI_FLAG=":ui"
fi

echo "================================================="
echo "      PyTransDatRO E2E Test Runner (Playwright)  "
echo "================================================="

npm --prefix webapp run "test:e2e${UI_FLAG}"
