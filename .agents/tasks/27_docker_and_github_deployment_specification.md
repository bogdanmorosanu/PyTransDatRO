# Task 27: Docker & GitHub Deployment Specification

## Objectives
- [x] Design container architecture and deployment strategy for PyTransDatRO on community infrastructure.
- [x] Establish open-source governance and multi-repository workflow between upstream (`bogdanmorosanu/PyTransDatRO`) and downstream (`geospatialorg/PyTransDatRO`).
- [x] Conduct deep technical review of reverse-proxy routing, rate-limiting IP preservation, upload limits, environment configuration, SQLite WAL backup safety, and log rotation.
- [x] Produce comprehensive specification document in `specs/deployment.md` ready for geospatial administrator review.
- [x] Provide user with step-by-step guidance for SSH key generation and sharing with geospatial.org administrator.

## Implementation Details

1. **Deployment Architecture (`specs/deployment.md`)**:
   - **2-Container Topology:**
     - `web`: Multi-stage build running `nginx:alpine` serving Vite static assets and reverse-proxying `/api/*`, `/docs`, and `/openapi.json` to backend.
     - `backend`: `python:3.11-slim` running FastAPI with embedded in-process `pytransdatro` geodetic calculation engine, mounting persistent SQLite volume.
   - **Persistent Volume Mount:**
     - Host directory `./logs` mounted to `/app/logs` inside backend container for `usage.db` (auto-created) and `legacy_stats.db` (seeded once via SCP).
   - **Safety & Production Hardening:**
     - Nginx client body limit: `client_max_body_size 50M;` matching `PYTRANSDAT_MAX_FILE_UPLOAD_BYTES`.
     - Reverse-proxy IP headers: `X-Forwarded-For` and `X-Real-IP` forwarding to prevent FastAPI sliding-window rate limiter from clustering all requests under Nginx container IP.
     - Container healthcheck synchronization: Nginx waits for backend healthcheck before serving traffic.
     - Docker log rotation: `max-size: 10m`, `max-file: 3` to prevent host disk exhaustion.
     - Non-blocking SQLite WAL live backup command documented.

2. **Multi-Repository Governance & Fork Model**:
   - Upstream: [`bogdanmorosanu/PyTransDatRO`](https://github.com/bogdanmorosanu/PyTransDatRO) under MIT license.
   - Downstream: [`geospatialorg/PyTransDatRO`](https://github.com/geospatialorg/PyTransDatRO) as an official GitHub Fork.
   - Workflow: Native 1-click Pull Requests from upstream to downstream. Upstream remains sovereign; geospatial maintains autonomy to merge or diverge.
   - Controlled deployment: SSH on-server command `git pull && docker compose up -d --build` avoids storing sensitive server secrets in GitHub Actions.

3. **SSH Key Protocol**:
   - Guided generation of modern `ed25519` key pair in Windows PowerShell.
   - Verified only public key (`id_ed25519.pub`) is shared with geospatial administrator.

## Verification
- Comprehensive spec audit performed covering FastAPI settings prefix (`PYTRANSDAT_`), Nginx proxy directives, rate-limiter dispatch logic, and Docker Compose restart/health dependencies.
- Confirmed zero modifications to core geodetic library or mathematical invariants.
