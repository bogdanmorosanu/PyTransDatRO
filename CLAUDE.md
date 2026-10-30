# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

PyTransDatRO is a Romanian geodetic coordinate transformation engine implementing the official ANCPI algorithm for bidirectional conversion between **Stereo 70** (national projection, meters) and **ETRS89** (lat/lon). It ships as three decoupled layers: a pure-Python core library, a FastAPI service, and a TypeScript/Leaflet web map.

## Architecture: Three Strictly Decoupled Layers

The single most important invariant in this repo is **layer discipline**. Changes must land in the correct layer, and dependencies flow in one direction only: `webapp/` → `api/` (HTTP) → `pytransdatro/` (import).

- **`pytransdatro/`** — Pure-Python core geodetic library. **Standard library only** (`math`, `struct`, `functools`, `abc`, `pathlib`, `typing`). The *sole* exception is `numpy`, allowed **strictly** for deserializing the binary `.spg` grid files — never for processing. This layer must **NEVER** import from `api/`, `fastapi`, `pydantic`, or any web dependency.
  - `trans_ro.py` — `TransRO` coordinator class; public entry points `st70_to_etrs89()` / `etrs89_to_st70()` (2D & 3D). Caches grid cells and polynomial params across calls for batch throughput.
  - `proj_stereo.py` — Oblique Stereographic projection. `trans_grid.py` — bicubic spline 2D shift interpolation + 1D geoid height. `trans_helmert2d.py` — 4-param 2D Helmert. `spg_reader.py` — binary grid reader + auto-discovery. `constants.py` — geodetic/projection invariants. `telemetry.py` — async SQLite WAL logger. `exceptions.py` — domain errors (`OutOfGridErr`, `NoDataGridErr`).
- **`api/`** — FastAPI service consuming `pytransdatro` as a library. Must contain no client-side business logic. `config.py` (pydantic-settings), `dependencies.py` (DI providers, `get_settings`), `middleware/rate_limit.py`, `routes/` (point, batch, file, feedback, info, health, legacy), `schemas/` (Pydantic request/response). Entry point: `api.main:app`.
- **`webapp/`** — Static TypeScript + Vite + Leaflet SPA consuming the API over HTTP. `src/components/` (panels, charts), `src/ui/` (reusable vanilla-TS control factories), `src/style/` (modular CSS), `src/utils/` (parsers, exporters), `src/config.ts`.

### Critical invariants (do not violate)

- **Zero numerical regression**: coordinate results must stay byte-for-byte identical to ANCPI benchmarks (< 10⁻⁹ tolerance). Verify with the ANCPI comparison tool after any math change.
- **Do not modify** the loading structure or content of binary grid files in `pytransdatro/grids/`.
- **Legacy endpoint** (`/transdatonline/cooOpService`) is a 100% byte-for-byte emulation of TransDatOnline for GIS plugins — preserve exact output format.
- **Configuration hierarchy**: geodetic invariants → `pytransdatro/constants.py` (never in `.env`); grid params → `.spg` binary header (dynamic); server/throttling ops → `api/config.py` + `.env` (prefix `PYTRANSDAT_`); web client → `webapp/src/config.ts` (prefix `VITE_`). Never hardcode ports, caps, timeouts, or paths in business logic.
- **No silent failures**: never `except:` or `except Exception: pass`. Catch specific exceptions and log via `logging.getLogger(__name__)` (not `print`).

`GLOSSARY.md` at the repo root defines the project terminology.

## Common Commands

```bash
# Python environment (conda env named `pytransdat`), or editable install:
pip install -e ".[api,dev]"

# Python tests (unit, API, telemetry)
python -m pytest -q
python -m pytest tests/test_trans_ro.py -q          # single module
python -m pytest tests/test_trans_ro.py::test_name  # single test

# Lint & format (line length 100; rules E/W/F/I)
python -m ruff check pytransdatro api tests tools
python -m ruff format --check pytransdatro api tests tools

# Frontend (run from repo root with --prefix)
npm --prefix webapp test            # Vitest unit tests
npm --prefix webapp run test:e2e   # Playwright UI tests (needs PYTHON=<pytransdat python>; one-time: npx --prefix webapp playwright install chromium)
npm --prefix webapp run build       # tsc typecheck + Vite production bundle (run before finishing any webapp task)
npm --prefix webapp run dev         # Vite dev server (http://localhost:5173)

# ANCPI benchmark regression (9,046 national control points, < 10⁻⁹ tolerance)
python tools/generate_ancpi_comparison_reports.py

# Run full stack locally
docker compose up --build           # unified service on :8000 (web + API + /docs)
python -m uvicorn api.main:app --reload --port 8000   # API only

# Local dev loop (two servers with hot reload)
./dev.sh                            # macOS/Linux launcher: API (:8000, --reload) + Vite (:5173); Ctrl+C stops both
dev.ps1                             # Windows/PowerShell equivalent with portable Python auto-detection
```

> **Dev launchers**: `dev.sh` (macOS/Linux) and `dev.ps1` (Windows) start the two-process dev loop — uvicorn with `--reload` plus the Vite dev server, which proxies `/api`, `/health`, `/docs`, `/openapi.json` to the API (see `webapp/vite.config.ts`). Use this for live editing. For a one-process run, build the frontend (`npm --prefix webapp run build`) and start the API alone — it mounts `webapp/dist` as static files and serves web + API + `/docs` on `:8000`.

## Testing Conventions

- Every feature or bugfix must ship with tests: `tests/*.py` (pytest) for core/API, and `webapp/src/**` (Vitest) for any change to `dmsParser.ts`, `pointParser.ts`, `gridBounds.ts`, `exporters.ts`, or `src/ui/` factories.
- Any change to a covered UI flow (single-point transform, batch paste, file upload/download, map centering/no-auto-pan, language switch, API error display) or to a rule in `.agents/rules/webapp-design.md` must add or update a Playwright e2e test (`webapp/e2e/`); pure styling/copy changes are exempt; Vitest remains the default for parsers and `src/ui/` factories.
- Parsers must tolerate empty lines, comments, and whitespace without throwing.
- **Windows sandbox note**: Starlette's `TestClient` uses loopback sockets; under a restricted Windows sandbox, running full `pytest` or any `TestClient`-importing module (`tests/test_api_*.py`, `tests/test_telemetry.py`) can hang — run those unsandboxed. Core/config tests, Ruff, and Vitest run fine sandboxed. (Not relevant on macOS/Linux.)

## Webapp Design Rules (see `.agents/rules/webapp-design.md`)

The web map has strict, opinionated UX invariants — read the rules file before UI work. Highlights: no decorative emojis in labels; no ad-hoc/inline styling (all styles derive from `variables.css` and modular stylesheets via `src/ui/` factories); 480px workbench panel width; coordinate transforms must **never** auto-pan/zoom the map; the app's own centering/zoom actions are bounded to Romania's envelope and never pan outside it for out-of-grid points (free user navigation is unrestricted); quasigeoid height is a measurement and must never be colored as a status (green/yellow/red are reserved for health/warning/error); strict 1-to-1 line correspondence in batch/file processing with 1-based `Point N` indexing.

## Repo-Specific Agent Workflow

- `.agents/rules/` holds the authoritative coding standards (code-quality, configuration, python-style, webapp-design) — consult them for any nontrivial change.
- `.agents/skills/` and `.agents/tasks/` are a pre-existing (Gemini-oriented) workflow: skills describe audit/verify/doc/commit routines, and each completed chat session is recorded as a numbered task file in `.agents/tasks/`. `GEMINI.md` is the counterpart to this file for that agent.
- `docs/` contains living architecture docs (`architecture.md`, `pytransdatro_algorithm.md`, `rest_api.md`, `grid.md`, `testing.md`, `telemetry.md`, `log.md`, `web_app.md`) — keep them in sync when architectural contracts change.
- The multi-phase refactor in `refactor_roadmap.md` is **100% complete** and is the immutable baseline; it does not constrain new feature work but its patterns must be upheld.
