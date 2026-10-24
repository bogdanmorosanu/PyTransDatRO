# Task 46: Production Deployment Hardening

## Objectives
- [x] Conduct an architectural audit of changes implemented between Task 27 and Task 45 to verify container deployment readiness.
- [x] Verify backward compatibility for URL query parameter hydration (`/?op=ETRS89ToStereo70&c1=46.00000000&c2=25.00000000&c3=100.000`).
- [x] Fix Nginx reverse proxy routing in `nginx.conf` to proxy `/health`, `/cooOpService`, `/transdatonline/*`, and `/redoc` to the backend.
- [x] Enable host `.env` configuration injection into the `backend` container via `docker-compose.yml` (`env_file: - path: .env; required: false`).
- [x] Add build-time `ARG` declarations for `VITE_APP_NAME` and `VITE_API_BASE_URL` in `Dockerfile.web`.
- [x] Clean up `requirements.txt` by purging the unused `pandas` dependency.
- [x] Update `.env.example` and `docs/deployment.md` with complete configuration references.
- [x] Expand `.github/workflows/ci.yml` to execute all API and telemetry test suites (`test_api_legacy_stats.py`, `test_logging.py`).
- [x] Verify complete quality gates locally: Ruff linter & formatter, 83 Vitest tests, Vite production bundle build, and 200 Python pytest tests.

---

## Implementation Details

1. **Production Reverse Proxy Routing (`nginx.conf`)**:
   - Expanded proxy match regex from `^/(api|docs|openapi.json)` to `^/(api|docs|redoc|openapi.json|health|cooOpService|transdatonline)`.
   - Guaranteed that third-party desktop clients and automated scripts calling legacy servlets receive direct backend JSON responses rather than being caught by Nginx's SPA fallback (`index.html`).
   - Allowed external health monitoring probes on the exposed web port to inspect `/health`.

2. **Container Orchestration & Environment Injection (`docker-compose.yml`, `Dockerfile.web`)**:
   - Added `env_file: - path: .env; required: false` to the `backend` service definition, enabling host `.env` configuration overrides to reach FastAPI without requiring `.env` inside the Docker image.
   - Added build-time arguments (`ARG VITE_APP_NAME`, `ARG VITE_API_BASE_URL`) to `Dockerfile.web` for flexible frontend packaging.

3. **Dependency & Configuration Hygiene (`requirements.txt`, `.env.example`, `docs/deployment.md`)**:
   - Removed unused `pandas` dependency from `requirements.txt`, keeping the backend image pure and saving ~40 MB.
   - Synchronized `.env.example` and `docs/deployment.md` with `PYTRANSDAT_LEGACY_DB_PATH`, `feedback.jsonl`, `telemetry_fallback.jsonl`, and optional SMTP parameters.

4. **CI Workflow Verification Completeness (`.github/workflows/ci.yml`)**:
   - Updated the API & Telemetry test step in GitHub Actions to include `tests/test_api_legacy_stats.py` and `tests/test_logging.py`.

---

## Verification & Testing

- **Docker Compose Config**: `docker-compose config` passed syntax and schema validation cleanly.
- **Ruff Code Audit**: `ruff check .` -> 0 errors; `ruff format --check` -> 55 files cleanly formatted.
- **Frontend Vitest Suite**: `npm --prefix webapp test` -> 83/83 passed across 6 test suites in 331ms.
- **Frontend Production Build**: `npm --prefix webapp run build` -> `dist/` bundle created in 490ms with zero errors.
- **Python Integration & Geodetic Parity Suite**: `pytest tests/` -> 200/200 passed in 14.62s.
