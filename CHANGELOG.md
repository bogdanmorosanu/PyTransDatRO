# Changelog

All notable changes to **PyTransDatRO** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] - 2026-10-30

### Initial Stable Release

#### Core Geodetic Engine (`pytransdatro`)
- **Pure Python Standard Library Architecture**: Zero external scientific C-extensions required for coordinate math.
- **Official ANCPI / RomGEO Parity**: Verified against 9,046 national benchmark points with zero regression ($< 10^{-9}$ numerical tolerance).
- **Transformation Pipelines**:
  - Stereo 70 (Dealul Piscului 1970 · EPSG:31700) $\longleftrightarrow$ ETRS89 (GRS80 · EPSG:4258).
  - 2D Planimetric Helmert & Interpolation and 3D Altimetric Quasigeoid ($\zeta = h - H$).
- **Standalone Packaging**: Bundled binary grid assets (`rom_grid3d_25.09.spg`), metadata, and comparison summaries for standard PyPI distribution.

#### REST API (`api`)
- **Modern REST Endpoints (`/api/v1`)**:
  - `POST /api/v1/transform/point`: Single-point high-speed transformation.
  - `POST /api/v1/transform/batch`: Vectorized batch point conversion with per-point diagnostics.
  - `POST /api/v1/transform/file`: Streaming file parser supporting CSV, TXT, and DAT formats.
  - `GET /api/v1/info/grid`: Active shift grid metadata and coverage envelope.
  - `GET /health`: Comprehensive health check with dynamic runtime operational limits.
- **Legacy Compatibility (`/transdatonline/cooOpService`)**:
  - 100% byte-for-byte emulation of the 2012 TransDatOnline Java servlet contract.
- **Resource Protection & Security**:
  - Sliding-window IP rate limiter with automated deque eviction.
  - Honeypot and mathematical challenge verification for public feedback submissions.
  - Asynchronous non-blocking telemetry logging with SQLite fallback spillover.

#### Geodetic Studio Web App (`webapp`)
- **Interactive Multi-Panel Interface**:
  - Leaflet-powered GIS map with Romania grid boundary visualization, orthophoto basemaps, and coordinate inspection.
  - Tabbed workflow supporting Single Point, Batch Paste, and Multi-File Upload.
  - Live cursor coordinates badge, distance measurement, and layer color customization.
- **Telemetry & Analytics Dashboard**:
  - Dedicated usage dashboard featuring computation volume timelines, historical TransdatOnline archive (2011–2026), and spatial density heatmaps.
- **Client-Side Export**:
  - One-click export to CSV, GeoJSON, and AutoCAD DXF formats.
- **Internationalization (i18n)**:
  - Complete bilingual Romanian (`ro`) and English (`en`) support with 100% key and placeholder parity.

#### Deployment & Developer Tooling
- **Single-Entrypoint Configuration**: Root `.env` as the sole configuration file across backend and web containers.
- **Container Architecture**: Docker Compose stack featuring multi-stage Python backend and Nginx reverse proxy frontend with dynamic `client_max_body_size`.
- **Portable Scripting**: Machine-agnostic PowerShell dev and test runners (`dev.ps1`, `test-all.ps1`, `test-e2e.ps1`).
- **Automated CI/CD**: GitHub Actions workflow running Ruff lint/format, full pytest integration suite, Vitest frontend unit tests, and Playwright E2E browser tests.
