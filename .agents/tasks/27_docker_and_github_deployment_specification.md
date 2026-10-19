# Task 27: Docker & GitHub Deployment Specification and Implementation

## Objectives
- [x] Design container architecture and deployment strategy for PyTransDatRO on community infrastructure.
- [x] Establish open-source governance and multi-repository workflow between upstream (`bogdanmorosanu/PyTransDatRO`) and downstream (`geospatialorg/PyTransDatRO`).
- [x] Conduct deep technical review of reverse-proxy routing, rate-limiting IP preservation, upload limits, environment configuration, SQLite WAL backup safety, and log rotation.
- [x] Produce comprehensive specification document in `specs/deployment.md`.
- [x] Create living technical documentation in `docs/deployment.md` and link from `docs/architecture.md`.
- [x] Implement Phase 1: 2-container Docker stack (`Dockerfile`, `Dockerfile.web`, `nginx.conf`, `docker-compose.yml`, `.env.example`).
- [x] Implement Phase 2: GitHub Actions CI workflow (`.github/workflows/ci.yml`).
- [x] Verify complete test suite locally (Ruff linting, 68 Vitest frontend tests, 62 Python core geodetic tests).

## Implementation Details

1. **Production Nginx Reverse Proxy (`nginx.conf`)**:
   - Configured port 80 HTTP server with `gzip` compression for static web assets.
   - Set `client_max_body_size 50M;` to match `PYTRANSDAT_MAX_FILE_UPLOAD_BYTES`.
   - Reverse proxies `/api/*`, `/docs`, and `/openapi.json` to `http://backend:8000`.
   - Injected client IP headers (`X-Real-IP`, `X-Forwarded-For`, `X-Forwarded-Proto`) ensuring FastAPI rate limiter distinguishes unique clients.
   - Configured 300-second proxy timeout buffers for streaming file conversions.

2. **Web Container Build (`Dockerfile.web`)**:
   - Multi-stage build:
     - Stage 1: `node:20-alpine` runs `npm ci` and `npm run build` targeting `webapp/`.
     - Stage 2: `nginx:alpine` packages pre-compiled assets into `/usr/share/nginx/html` with `nginx.conf`.
   - Includes container healthcheck using `wget`.

3. **Backend Container Build (`Dockerfile`)**:
   - Streamlined into a dedicated Python 3.11-slim backend runtime.
   - Installs runtime dependencies, copies `pytransdatro/` (including binary grids) and `api/`.
   - Includes healthcheck testing `/health` endpoint to verify geodetic grid pre-warming.

4. **Multi-Container Orchestration (`docker-compose.yml`)**:
   - Defines `backend` and `web` services communicating over the internal Docker network.
   - Maps host port `${WEBAPP_PORT:-8000}:80` (default 8000 locally, configurable to 8008 on geospatial server).
   - Enforces healthcheck dependency (`web` waits for `backend` to report `service_healthy`).
   - Configures persistent volume mount `./logs:/app/logs`.
   - Implements Docker JSON log rotation (`max-size: 10m`, `max-file: 3`).

5. **GitHub Actions CI Automation (`.github/workflows/ci.yml`)**:
   - `python-ci` job: Runs Ruff linter, formatting checks, core ANCPI geodetic parity tests, and API integration tests on Python 3.11.
   - `frontend-ci` job: Runs TypeScript compile, Vite production bundle build, and Vitest suite on Node.js 20.

## Verification & Testing
- **Linter & Formatting**: `ruff check .` passed with zero errors.
- **Frontend Unit Tests**: `npm --prefix webapp test` -> 68/68 passed (5 test files).
- **Core Geodetic Suite**: `pytest tests/test_trans_ro.py tests/test_config.py` -> 62/62 passed in 3.34s.
