---
name: doc-record
description: Creates or updates architectural and technical documentation in the docs/ folder. Use whenever modifying or refactoring major components (core engine, REST API, web app, logging, telemetry) to keep documentation in sync with the codebase.
---

# Documentation Record Skill (`doc-record`)

## Purpose
This skill ensures that living technical documentation in `docs/` remains the authoritative source of truth for the PyTransDatRO system. Unlike `.agents/tasks/` (which logs chronological session history), `docs/` defines the **current architectural reality, algorithms, schemas, and public contracts** of the application.

Use this skill whenever code refactoring, schema changes, or UI enhancements alter how a subsystem works.

---

## Documentation Mapping

When working on a subsystem, update the corresponding document in `docs/`:

| Subsystem / Layer | Authoritative Document | Contents to Maintain |
| :--- | :--- | :--- |
| **System Architecture** | [`docs/architecture.md`](../docs/architecture.md) | High-level layer topology, class hierarchy, pipeline flow, grid discovery. |
| **Core Geodetic Engine** | [`docs/core_engine.md`](../docs/core_engine.md) | Oblique Stereographic formulas, 2D Helmert, bicubic spline interpolation, quasigeoid undulation, dual-CRS 6-point extended bounding boxes, 4-stage batch pipeline, domain exceptions. |
| **Grid Management** | [`docs/grid.md`](../docs/grid.md) | `.spg` binary format specification, header offsets, planar shift subgrids. |
| **REST API Layer** | [`docs/rest_api.md`](../docs/rest_api.md) | Modern endpoints (`/api/v1/*`), Pydantic models, query parameters, error schemas, legacy `/transdatonline` emulation. |
| **Web Application & UI** | [`docs/web_app.md`](../docs/web_app.md)<br>[`webapp/design_guide.md`](../webapp/design_guide.md) | Web GIS Leaflet architecture, workbench tabs, UI factory helpers, modular CSS tokens, i18n support. |
| **Telemetry & Logging** | [`docs/logging.md`](../docs/logging.md) | SQLite asynchronous WAL logger schema, retention policies, Python `logging` dictConfig setup. |
| **Testing & Quality** | [`docs/testing.md`](../docs/testing.md) | Test matrix, numerical precision benchmarks, Vitest parser tests. |

---

## Operating Principles

1. **Keep in Sync (No Documentation Rot)**:
   - When adding or changing configuration options, updating endpoint responses, or refactoring classes, update the relevant `docs/*.md` file in the same PR or chat session.
2. **Clear Mathematical & Technical Precision**:
   - Preserve mathematical formulas in LaTeX markdown (e.g., $h = H + N$).
   - Provide concrete JSON request/response snippets for API endpoints.
   - Use Mermaid diagrams to visualize class relationships or pipeline sequences.
3. **Link Verification**:
   - Ensure markdown links point to valid relative paths (e.g. `[architecture.md](architecture.md)`).
4. **No Meta Chatter**:
   - Do not mention AI agents, chat numbers, or task files inside `docs/`. Write in an authoritative, technical, engineering tone.
