# PyTransDatRO

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.110+-green.svg)](https://fastapi.tiangolo.com)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**PyTransDatRO** is a high-performance Romanian geodetic transformation engine and modern web service implementing the official algorithm published by the Romanian **National Agency for Cadastre and Land Registration (ANCPI)**.

It provides bidirectional coordinate transformations between **Stereo 70** (National Projection) and **ETRS89** (European Terrestrial Reference System 1989), featuring:
- **Pure Python Core Library (`pytransdatro`)**: Zero external compiled C/GIS dependencies for numerical processing; fast in-memory binary grid lookup (`.spg`).
- **Modern REST Web API (`api/`)**: High-throughput FastAPI service supporting single points, batch transformations, CSV/TXT file conversion, and 100% backward compatibility with legacy TransDatOnline endpoints.
- **Interactive Web Map (`webapp/`)**: Responsive TypeScript + Leaflet web application for visual point transformation, map-based coordinate inspection, and file processing.
- **Telemetry & Monitoring**: Zero-overhead SQLite WAL background telemetry tracking transformation usage, bounding boxes, and performance metrics.

---

## Quickstart with Docker (Recommended for Reviewers & Developers)

The easiest way to evaluate and run the complete PyTransDatRO stack (Web Map + REST API + Geodetic Engine) without configuring Python, Node.js, or virtual environments is with **Docker Compose**:

### 1. Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows / macOS) or Docker Engine + Docker Compose (Linux).

### 2. Launch the Application
Clone the repository (or download the source archive) and run:

```bash
docker compose up --build
```

Docker will automatically:
1. Compile the TypeScript/Vite frontend into a production distribution.
2. Build the Python 3.11 environment with all API dependencies.
3. Launch the unified service on port `8000`.

### 3. Access the Interfaces
Once started, open your browser:

- **Interactive Web Map**: [http://localhost:8000](http://localhost:8000)
- **Interactive API Documentation (Swagger UI)**: [http://localhost:8000/docs](http://localhost:8000/docs)
- **Alternative API Docs (ReDoc)**: [http://localhost:8000/redoc](http://localhost:8000/redoc)
- **Health Check Endpoint**: [http://localhost:8000/health](http://localhost:8000/health)

To stop the container:
```bash
docker compose down
```

---

## Core Library Usage (Python)

If you are using `pytransdatro` as a standalone Python library in your own code:

```python
import pytransdatro

# Initialize the transformation engine (loads the .spg binary grid into memory once)
t = pytransdatro.TransRO()

# ----------------------------------------------------
# 1. Stereo 70 -> ETRS89 (2D & 3D)
# Input: North (m), East (m)[, Height (m)]
# Output: Latitude (rad), Longitude (rad)[, Ellipsoidal Height (m)]
# ----------------------------------------------------
lat, lon = t.st70_to_etrs89(693771.731, 310723.518)
lat, lon, h = t.st70_to_etrs89(693771.731, 310723.518, 122.714)

# ----------------------------------------------------
# 2. ETRS89 -> Stereo 70 (2D & 3D)
# Input: Latitude (rad), Longitude (rad)[, Ellipsoidal Height (m)]
# Output: North (m), East (m)[, Height (m)]
# ----------------------------------------------------
n, e = t.etrs89_to_st70(0.832795488097716, 0.3922724455407368)
n, e, h = t.etrs89_to_st70(0.832795488097716, 0.3922724455407368, 162.0874)
```

> **Note on Performance**: The engine is optimized for batch transformations. Grid cells (11 km 2D cells and 0.0333° 3D elevation cells) and polynomial parameters are cached across calls on the `TransRO` instance, maximizing throughput when transforming points within the same region.

---

## REST API Overview

The service exposes modern JSON, streaming, and backward-compatible endpoints:

| Endpoint | Method | Description |
| :--- | :--- | :--- |
| `/api/v1/transform/point` | `POST` | Transform single coordinates (supports deg/rad and 2D/3D). |
| `/api/v1/transform/batch` | `POST` | High-throughput batch coordinate transformation. |
| `/api/v1/transform/file` | `POST` | Upload and transform CSV/TXT coordinate files with streaming download. |
| `/api/v1/feedback` | `POST` | Submit discrepancy reports and user feedback. |
| `/api/v1/grid/info` | `GET` | Active SPG grid metadata, distortion models, and version information. |
| `/api/v1/telemetry/timeline` | `GET` | Queryable time series analytics filtered by client source and date window. |
| `/api/v1/telemetry/spatial-matrix` | `GET` | Macroscopic 0.1° density matrix and quadrant distribution. |
| `/health` | `GET` | Service liveness, readiness, and grid verification. |
| `/transdatonline/cooOpService`<br/>`/cooOpService` | `GET`/`POST` | **Legacy TransDatOnline compatibility** (100% byte-for-byte emulation for GIS plugins and scripts). |

---

## Local Developer Setup & Quality Workflow

### 1. Python Environment Setup
Using Conda (recommended, matching project environment):
```bash
conda activate pytransdat
# or install dependencies directly:
pip install -e .[api,dev]
```

### 2. Run Test Suites
```bash
# Run all 98 Python unit, API, and telemetry tests
python -m pytest -q

# Run all 65 Frontend Vitest unit tests
npm --prefix webapp test
```

### 3. Run Static Analysis & Linter Audits
```bash
# Check code style and formatting across the repository
python -m ruff check pytransdatro api tests tools
python -m ruff format --check pytransdatro api tests tools

# Verify TypeScript compilation and production bundle
npm --prefix webapp run build
```

### 4. Run ANCPI Benchmark Verification
```bash
# Validate zero-regression (< 10^-9 tolerance) across 9,046 national control points
python tools/generate_ancpi_comparison_reports.py
```

---

## Project Structure

```
PyTransDatRO/
├── api/                      # Decoupled FastAPI web service layer
│   ├── config.py             # Centralized Pydantic BaseSettings & logging setup
│   ├── dependencies.py       # FastAPI dependency injection providers
│   ├── middleware/           # Sliding-window rate limiting & CORS
│   ├── routes/               # REST routers (point, batch, file, feedback, info, health, legacy)
│   ├── schemas/              # Pydantic models (request/response schemas, ErrorResponse)
│   └── services/             # Background services (mailer)
├── pytransdatro/             # Pure-Python core geodetic library (standard library only)
│   ├── constants.py          # WGS84 and Stereo 70 geodetic & projection constants
│   ├── grids/                # ANCPI binary transformation grid (.spg) and metadata
│   ├── trans_ro.py           # Coordinator transformation engine class (TransRO)
│   ├── trans_grid.py         # Bicubic spline interpolation (2D shift) & 1D geoid
│   ├── proj_stereo.py        # Oblique Stereographic projection engine
│   ├── trans_helmert2d.py    # 4-parameter 2D Helmert transformation
│   ├── spg_reader.py         # Binary grid reader & auto-discovery
│   └── telemetry.py          # Asynchronous SQLite WAL telemetry logger
├── webapp/                   # Interactive TypeScript + Vite + Leaflet web map
│   ├── src/components/       # UI workbench panels & analytics charts
│   ├── src/ui/               # Accessible TypeScript UI factory helpers
│   ├── src/style/            # Pure Vanilla Modular CSS (layout, forms, buttons, badges, menus, telemetry)
│   └── src/utils/            # Client parsers (DMS, delimiters), exporters (CSV, GeoJSON, DXF), bounds
├── tests/                    # Pytest test suite (98 tests) & reference datasets
├── tools/                    # Benchmarking, metadata extraction, and ANCPI comparison reports
├── docs/                     # Living technical architecture, API, and testing documentation
├── .env.example              # Documented configuration template
├── pyproject.toml            # Project build metadata and Ruff linting rules
├── Dockerfile                # Multi-stage production container build
├── docker-compose.yml        # One-command developer review launcher
└── requirements-api.txt      # FastAPI & web service dependencies
```

---

## Upcoming / Roadmap (v1.0.0 Release)

The following sections will be finalized for the public v1.0.0 release:

- [ ] **PyPI Package Distribution** (`pip install pytransdatro`)
- [ ] **ANCPI Official Benchmark & Accuracy Validation Report** (millimeter-level validation dataset results)
- [ ] **Pre-built Container Registry Images** (`ghcr.io/bogdanmorosanu/pytransdatro`)
- [ ] **QGIS Plugin Integration Guide**

---

## License

This project is licensed under the [MIT License](LICENSE).

