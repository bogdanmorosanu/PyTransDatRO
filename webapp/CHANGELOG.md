# Changelog - PyTransDatRO Web Application (`webapp/`)

All notable changes to the PyTransDatRO web application and its interactive Web GIS interface will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

---

## [1.0.0] - 2026-10-10

### Added
- **Interactive Leaflet Web GIS (`MapView.ts`)**:
  - Basemap selector supporting OpenStreetMap, CartoDB Positron, CartoDB Dark Matter, and ESRI World Imagery.
  - Automatic boundary overlay displaying the official Romanian TransDatRO grid bounding box (`[43.5, 20.2]` to `[48.3, 29.8]`).
  - Marker clustering via `leaflet.markercluster` for performant point rendering.
  - Automatic map panning and zooming to bounds when points are transformed.
  - Reverse geodetic click inspection: clicking on the map displays geographic (WGS84/ETRS89) and projected (Stereo 70) coordinates.
- **Single Point Converter (`PointConverter.ts`)**:
  - Live transformation between Stereo 70 and ETRS89 (both directions).
  - Flexible coordinate input format support (Decimal Degrees, Radians, DMS string parsing).
  - High-precision results display with quasigeoid undulation $\zeta$ and grid distortion shifts $\Delta E, \Delta N$.
  - "Show on Map" button and direct marker placement.
- **Batch Converter & Table Preview (`BatchConverter.ts`)**:
  - Multi-point text input with auto-delimiter detection (space, tab, comma, semicolon).
  - Interactive results table with instant map synchronization.
  - Client-side export to CSV, GeoJSON, and AutoCAD DXF.
- **File Upload (`FileUpload.ts`)**:
  - Drag-and-drop file processor with streaming progress bar.
- **Guardrails & Modals**:
  - `VolumeWarningModal.ts`: Performance warning prompt for transformations with > 2,000 points (full vs 1,000-point sample vs skip map).
  - `ReportModal.ts`: Diagnostic calculation report generator.
  - `FeedbackModal.ts`: Discrepancy reporting.
  - `StatsView.ts`: Telemetry visualization and grid model metadata viewer.
- **Design System & Internationalization**:
  - Vanilla CSS design system with Dark / Light / Auto mode support (`variables.css`, `base.css`, `components.css`, `map.css`).
  - Bilingual i18n support for Romanian (`ro`, default) and English (`en`).
- **REST Client & Proxy (`apiClient.ts` / `vite.config.ts`)**:
  - Type-safe REST client communicating with FastAPI endpoints under `/api/v1/`.
  - Vite reverse proxy forwarding `/api`, `/health`, and `/docs` to `http://127.0.0.1:8000`.
