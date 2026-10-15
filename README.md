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

The service exposes modern JSON and streaming endpoints:

| Endpoint | Method | Description |
| :--- | :--- | :--- |
| `/api/v1/transform/point` | `POST` | Transform single coordinates (supports deg/rad and 2D/3D). |
| `/api/v1/transform/batch` | `POST` | High-throughput batch coordinate transformation. |
| `/api/v1/transform/file` | `POST` | Upload and transform CSV/TXT coordinate files. |
| `/api/v1/health` | `GET` | Service liveness, readiness, and grid verification. |
| `/api/v1/info` | `GET` | Grid metadata, distortion models, and version information. |
| `/transdat/service.php` | `GET`/`POST` | **Legacy TransDatOnline compatibility** for existing GIS plugins and scripts. |

---

## Project Structure

```
PyTransDatRO/
├── api/                  # FastAPI web service application & routing layer
│   ├── routes/           # REST endpoints (point, batch, file, legacy, health)
│   └── middleware/       # Sliding-window rate limiting & CORS
├── pytransdatro/         # Pure-Python core geodetic library
│   ├── grids/            # ANCPI binary transformation grid (.spg) and metadata
│   ├── trans_ro.py       # Main TransRO transformation engine class
│   └── telemetry.py      # Background SQLite WAL telemetry logger
├── webapp/               # Interactive TypeScript + Vite + Leaflet web map
├── tests/                # Unit, integration, and ANCPI accuracy tests
├── Dockerfile            # Multi-stage production container build
├── docker-compose.yml    # One-command developer review launcher
└── requirements-api.txt  # FastAPI & web service dependencies
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
