#!/usr/bin/env bash
#
# PyTransDatRO Full Verification Suite (macOS / Linux).
#
# Usage:
#   ./test-all.sh            # runs python, vitest, and playwright
#   ./test-all.sh --skip-e2e # runs python and vitest only
#

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

if [[ -x "$ROOT_DIR/.venv/bin/python" ]]; then
  PYTHON="$ROOT_DIR/.venv/bin/python"
else
  PYTHON="python3"
fi
export PYTHON

echo "================================================="
echo "       PyTransDatRO Full Verification Suite      "
echo "================================================="

echo -e "\n[1/3] Running Python Core & REST API Suite (pytest)..."
"$PYTHON" -m pytest -m "not legacy_live" -q

echo -e "\n[2/3] Running Frontend Unit & Parser Suite (vitest)..."
npm --prefix webapp test

if [[ "${1:-}" != "--skip-e2e" ]]; then
  echo -e "\n[3/3] Running Browser End-to-End Suite (playwright)..."
  npm --prefix webapp run test:e2e
fi

echo -e "\n================================================="
echo "   [✓] ALL SUITES PASSED SUCCESSFULLY           "
echo "================================================="
