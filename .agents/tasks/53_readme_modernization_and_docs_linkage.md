# README Modernization & Documentation Alignment

## Objectives
- [x] Audit [`README.md`](../../README.md) for synchronization with actual repository implementation and architecture.
- [x] Eliminate hardcoded/dynamic test counts and sample numbers to prevent documentation rot across releases.
- [x] Retire obsolete "Upcoming / Roadmap (v1.0.0 Release)" items (PyPI distribution, QGIS plugin, pre-built container registry image, and standalone benchmark reports).
- [x] Establish an in-repo "Documentation & Technical Reference" section linking to authoritative living documentation in [`docs/`](../../docs/), centered on Option A (in-repo living docs with native GitHub Markdown / LaTeX math formula rendering).
- [x] Verify code quality gates and formatting across the project.

## Implementation Details

1. **Test & Sample Count Modernization**:
   - Removed hardcoded test quantities (`189 Python unit...`, `123 Frontend Vitest...`, `207 tests`) from the developer setup and project tree sections in [`README.md`](../../README.md).
   - Generalized ANCPI benchmark validation instructions to refer to official national control points rather than fixed point counts.

2. **Roadmap Deprecation & Alignment**:
   - Retired the "Upcoming / Roadmap" section. Evaluated roadmap candidates:
     - **PyPI Package**: Dropped in favor of direct Git repository dependency (`pip install git+https://...`) or cloning, eliminating external packaging overhead.
     - **QGIS Plugin**: Dropped as out of scope for the core engine; legacy `/transdatonline` HTTP compatibility already serves third-party clients.
     - **Pre-built Container Registry**: Clarified difference between local `docker compose up --build` and public registry hosting (`ghcr.io`); dropped until formal CI release publishing is desired.
     - **Benchmark Reports**: Retired from roadmap since comprehensive benchmarking against official ANCPI RomGEO data is already integrated in [`tools/generate_ancpi_comparison_reports.py`](../../tools/generate_ancpi_comparison_reports.py).

3. **Authoritative Living Documentation Linkage**:
   - Replaced the roadmap with a structured **Documentation & Technical Reference** section in [`README.md`](../../README.md):
     - **[Core Geodetic Mathematical Model](docs/core_engine.md)**: Bidirectional transformations (**Stereo 70 $\longleftrightarrow$ ETRS89**), Krasovsky 1940 / GRS80 ellipsoids, Oblique Stereographic projection, 2D Helmert transformation, 16-node bicubic spline interpolation, 1D quasigeoid undulation elevation handling ($H_{\text{MN75}} \longleftrightarrow h_{\text{ETRS89}}$), and 6-point dual-CRS boundary envelope theory.
     - **[System Architecture](docs/architecture.md)**: Decoupled three-tier layer architecture.
     - **[REST API Specification](docs/rest_api.md)**: JSON schemas, batch streaming, and TransDatOnline emulation.
     - **[Grid Lifecycle & Management](docs/grid_lifecycle.md)**: ANCPI `.spg` grid structure, auto-discovery, and update workflows.
     - **[Testing & Numerical Verification](docs/testing.md)**: Quality gates, automated test framework, and RomGEO ground-truth verification.
     - **[Web Application Guide](docs/webapp/README.md)**: TypeScript client workbench architecture and design system.

## Verification & Quality Gates
- **Linter & Formatting**: Executed `ruff check` across `pytransdatro`, `api`, `tests`, and `tools` (clean, 0 warnings/errors).
- **Test Integrity**: Verified pytest suite (207 passing tests) and frontend Vitest suite (123 passing tests).
- **Link Integrity**: Verified workspace-relative markdown paths resolve accurately.
