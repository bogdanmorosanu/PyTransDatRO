# PyTransDatRO Architecture

## Overview
`pytransdatro` provides bidirectional coordinate transformations between **Stereo70** (Romanian projected CRS) and **ETRS89** (geographic CRS). The library is designed to match the output of the official CNC / Romgeo application byte-for-byte, supporting unified 2D/3D transformations via modern `.spg` grid packages.

## Class Structure
- **`TransRO`** (`trans_ro.py`): The main coordinator class. It composes instances of grid, Helmert, and projection classes to perform the full mathematical pipeline.
- **`SpgReader`** (`spg_reader.py`): Singleton reader that auto-discovers and deserializes the active `.spg` grid file, caching flat Python coordinate arrays and metadata parameters.
- **`Grid` (Abstract)** (`trans_grid.py`): Base class for grid coordinate validation and interpolation caching.
  - **`BiInterp`**: Nested class for evaluating the 16-coefficient bicubic spline polynomial on a 4x4 subgrid.
  - **`Grid1D`**: Concrete class for height/elevation corrections (Z-shift / quasigeoid undulation) using nearest-neighbor colocate indexing.
  - **`Grid2D`**: Concrete class for northing/easting planimetric distortion corrections (N, E shifts) using bicubic spline interpolation.
- **`Helmert2D`** (`trans_helmert2d.py`): Encapsulates the 4-parameter 2D Helmert transformation constants (Translations, PPM/Scale, Rotation) between Stereo70 and StereoGRS80, loaded dynamically from the `.spg` metadata.
- **`StereoProj`** (`proj_stereo.py`): Implements Oblique Stereographic projection on the WGS84 ellipsoid.
  - **`_Ellipsoid`**: Handles basic ellipsoid element calculations.
  - **`_ConfSph`**: Defines the conformal sphere and computes conformal latitudes and longitudes.

## Pipeline Flow
### Stereo70 -> ETRS89 (`st70_to_etrs89`)
1. **2D Grid Transformation**: Applies `geodetic_shifts` planar distortion shifts from the SPG grid (subtraction).
2. **2D Helmert Transformation**: Converts Stereo70 to StereoGRS80 using dynamic SPG parameters.
3. **Stereographic Projection**: Converts grid coordinates to geographic coordinates `(Lat, Long)`.
4. **Height Correction (Optional)**: If a `Z` coordinate is provided, the 1D geoid layer (`geoid_heights`) provides the quasigeoid undulation $N$ via nearest-neighbor lookup: $h = H + N$.

### ETRS89 -> Stereo70 (`etrs89_to_st70`)
1. **Stereographic Projection**: Transforms geographic coordinates to StereoGRS80 grid coordinates.
2. **2D Helmert Transformation (Inverse)**: Converts StereoGRS80 to Stereo70.
3. **2D Grid Transformation**: Applies `geodetic_shifts` planar distortion shifts (addition).
4. **Height Correction (Optional)**: Subtracts height correction via the 1D geoid layer: $H = h - N$.

## Grid Management
The system automatically discovers and loads the single `.spg` grid located in `pytransdatro/grids/`. If multiple grids or no grids are found, explicit descriptive exceptions (`AmbiguousGridError`, `MissingGridError`) are raised to prevent transformation ambiguity. Legacy binary `.grd` files have been retired.

## Geodetic Constants (`pytransdatro/constants.py`)
All physical geodetic parameters and projection invariants are extracted into pure standard-library constants:
- **WGS84 Ellipsoid**: Semi-major axis $a = 6,378,137.0\text{ m}$, flattening $f = 1 / 298.257223563$, eccentricity squared $e^2$.
- **Stereo 70 Projection Constants**: Central meridian $\lambda_0 = 25^\circ\text{ E}$, origin latitude $\phi_0 = 46^\circ\text{ N}$, scale reduction factor $k_0 = 0.99975$, false easting $E_0 = 500,000.0\text{ m}$, false northing $N_0 = 500,000.0\text{ m}$.
- Conformal sphere constants and radius of Gaussian mean curvature ($R = \sqrt{MN}$).
- See [docs/core_engine.md](core_engine.md) for the complete mathematical formulation, dual-CRS 6-point extended bounding boxes, 4-stage fault-isolated batch pipeline, and domain exception hierarchy.

## Web API Layer (`api/`)
A decoupled REST service built on FastAPI provides public access to the `pytransdatro` transformation pipeline:
- **Centralized Configuration (`api/config.py`)**: Environment-driven settings using Pydantic `BaseSettings` (`.env` file support), configuring server ports, CORS origins, SQLite telemetry directories, sliding-window rate limit thresholds, upload file sizes, and logging levels.
- **Structured Error Handling (`api/schemas/error.py`)**: Uniform `ErrorResponse` schema (`{"detail": "...", "code": "..."}`) with custom exception handlers for `StarletteHTTPException` and `RequestValidationError` (`code: "VALIDATION_ERROR"`).
- **Modern Endpoints (`/api/v1`)**: Fast JSON point, batch, streaming file (`/transform/file`), feedback (`/feedback`), and telemetry endpoints supporting radians and decimal degrees.
- **Legacy Compatibility (`/transdatonline/cooOpService` and `/cooOpService`)**: 100% byte-for-byte backward compatible emulation of the 2012 TransDatOnline Java service.
- **Dependency Injection (`api/dependencies.py`)**: Clean inversion of control for `Settings`, `TransRO` singleton engine, and `SqliteUsageLogger`.
- See [docs/rest_api.md](rest_api.md) for endpoint references and [docs/transdatonline/](transdatonline/README.md) for the historical reverse-engineering analysis.

## Web Application Layer (`webapp/`)
A responsive Web GIS client interface built with Vite, TypeScript, Vanilla CSS, and Leaflet:
- **Interactive Workbench**: Single point, multi-line batch, and drag-and-drop file upload workbenches.
- **Embedded Web GIS**: Leaflet map with basemap selector, Romania SPG grid boundary, custom segmented GIS scale bar, and marker clustering.
- **Modular Stylesheet Architecture (`webapp/src/style/`)**: Pure CSS custom properties decoupled across `layout.css`, `forms.css`, `buttons.css`, `badges.css`, `menus.css`, `telemetry.css`, and `variables.css`.
- **Reusable UI Factory Functions (`webapp/src/ui/`)**: Accessible component factories (`renderDropdownMenu`, `renderSegmentedControl`, `renderSelect`, `renderCopyBadge`, `renderCopyButton`) enforcing ARIA attributes and keyboard navigation.
- **Client-Side Utilities**: DMS parsing, auto-delimited point parsing, Point ID preservation, and instant multi-format exports (CSV, GeoJSON, AutoCAD DXF).
- **Production Serving**: Hosted via high-performance Nginx reverse proxy serving the Vite static bundle while proxying API endpoints to the backend container.
- See [docs/webapp/](webapp/README.md) for full architectural specifications, design system tokens, component archetypes, and actionable inventory, and [webapp/README.md](../webapp/README.md) for frontend development commands.

## Deployment & Container Architecture (`docs/deployment.md`)
Production deployments follow an environment-agnostic 2-container topology orchestrated by Docker Compose:
- **`web` (Nginx Alpine)**: Serves pre-compiled Vite static assets, enforces 50 MB file upload caps, and preserves client IP headers (`X-Forwarded-For`, `X-Real-IP`) for sliding-window rate limiting.
- **`backend` (Python 3.11-slim)**: Runs FastAPI with embedded in-process `pytransdatro` geodetic math and persistent SQLite WAL volume mount (`./logs`).
- **Universal Agnosticism**: Configuration is strictly decoupled using `.env.example` templates and host-local untracked `.env` files.
- See [docs/deployment.md](deployment.md) for full container topology, reverse proxy rules, and operational runbooks.

## Version Management & Release Automation
The project maintains a strict Single-Source-of-Truth (SSoT) pattern for application and library versions:
- **Core Version Definition**: Tracked synchronously across `pytransdatro/__init__.py` (`__version__`), `pyproject.toml` (`version`), and `webapp/package.json` (`version`) in accordance with Semantic Versioning 2.0.0.
- **Dynamic Runtime Exposure**:
  - **FastAPI Backend**: `FastAPI(..., version=pytransdatro.__version__)` and `/health` probe dynamically export `pytransdatro.__version__`.
  - **Vite Web Application**: `webapp/vite.config.ts` injects `__APP_VERSION__` from `package.json` at build time into `APP_CONFIG.appVersion`, dynamically rendered across the UI brand badge (`v${APP_CONFIG.appVersion}`) and feedback reports.
- **Synchronization Tooling**: `tools/bump_version.py` validates SemVer compliance, verifies cross-file parity (`--check`), and atomically updates all metadata files.
- **Release Git Tagging**: Managed via the `task-commit` agent skill, creating conventional release commits (`chore(release): v<version>`) and annotated Git tags (`v<version>`).

## Quality Assurance & Testing Matrix
The project enforces a cross-layer automated testing and verification strategy:
- **Core Geodetic Unit Suite (`tests/test_trans_ro.py`)**: 54 tests covering forward, inverse, 2D/3D roundtrips, half-grid verification, 6-point extended dual-CRS bounding boxes, `OutOfGrid2DErr` vs `OutOfGrid1DErr`, 4-stage `transform_batch` fault isolation, strict DMS grammar, and no-data boundary exceptions.
- **Web API & Telemetry Suite (`tests/test_api_*.py`, `tests/test_telemetry.py`, `tests/test_logging.py`, `tests/test_config.py`, `tests/test_feedback.py`)**: 146 tests covering modern REST endpoints, interval overlap period snapping, sentinel time windows, Option B batch fault isolation, legacy servlet emulation, sliding-window rate limiting, configuration defaults, feedback mailer, and SQLite WAL telemetry persistence (200 Python tests total).
- **Frontend Utility & UI Suite (`webapp/src/utils/*.test.ts`, `webapp/src/ui/*.test.ts`)**: 80 unit tests executed via Vitest verifying unspaced DMS parsing, coordinate auto-delimiters, Romania grid envelope bounds, pure string exporter generators, and accessible UI controls.
- **Official ANCPI Benchmark Verification (`tools/generate_ancpi_comparison_reports.py`)**: Automated verification against 9,046 national control points ensuring $< 10^{-9}$ tolerance and zero functional regression.
- **Static Analysis & Linters**: Zero warnings across Ruff (`ruff check .`, `ruff format --check .`) and TypeScript compiler (`tsc`).


