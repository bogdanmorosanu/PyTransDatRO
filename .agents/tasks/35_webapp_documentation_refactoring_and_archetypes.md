# Task 35: Web Application Documentation Refactoring & Component Archetypes

## Objectives
- [x] Create a dedicated, modular documentation directory `docs/webapp/` with a centralized `README.md` master index and task-based reading guide.
- [x] Author `docs/webapp/design_system.md` consolidating the GitHub Primer neutral dual-theme palette, Option #8 Subtle Surface+Border primary button, C2 Graticule brand logo, and map glass tokens.
- [x] Author `docs/webapp/archetypes.md` cataloging 14 component archetypes with member instances, DOM selectors, micro-typography, and styling invariants.
- [x] Author `docs/webapp/architecture.md` detailing the client-side SPA topology, unidirectional state flow, geodetic parsing pipelines, Leaflet multi-layer GIS, and i18n isolation.
- [x] Author `docs/webapp/actionable_inventory.md` migrating the complete living functional specification with canonical naming (`Name (#id)`) and cross-links.
- [x] Author `scratch/doc_reconciliation_report.md` with explicit data disposition tags (`[DATA CONSOLIDATED]`, `[DATA CONFLICT RESOLVED]`, `[DATA RETIRED - HISTORICAL]`) and rationale remarks for collaborative review.
- [x] Retire legacy fragmented documents: `webapp/changelog.md`, `docs/ui_actionable_inventory.md`, `webapp/design_guide.md`, and `docs/web_app.md`.
- [x] Align agentic framework and rules (`GEMINI.md`, `.agents/rules/webapp-design.md`, `.agents/skills/task-init/SKILL.md`, `.agents/skills/task-record/SKILL.md`, `webapp/README.md`, `docs/architecture.md`).
- [x] Validate all relative markdown links, machine-path portability, Vitest suite (80/80 passed), and Vite production bundle (0 errors).

---

## Implementation Details

1. **Modular Living Documentation Suite (`docs/webapp/`)**:
   - `docs/webapp/README.md`: Central navigation hub, directory sitemap, role-based reading guide for agents/engineers, canonical naming conventions, and golden architectural invariants.
   - `docs/webapp/design_system.md`: Visual tokens, Primer Neutral dual-theme specifications, Option #8 CTA button standard, C2 Map Projection Graticule vector logo, and Glassmorphism HUD tokens.
   - `docs/webapp/archetypes.md`: 14 component archetypes (Action CTA buttons, Segmented controls, Dropdown selects, Text inputs, Output read-only fields, Selectable copy badges, KPI summary cards, Visual status badges, Continuous presets, Multi-choice checkboxes, Date pickers, Export download buttons, Floating modals, and Chart series total switches).
   - `docs/webapp/architecture.md`: SPA client topology, reactive state lifecycle, geodetic math pipelines (Point, Batch, File), Leaflet 3-layer GIS architecture, and i18n locale number formatting.
   - `docs/webapp/actionable_inventory.md`: Complete functional inventory across Sections I–V and Section 7 Coverage Gap Analysis (14 items), bridging product requirements and automated E2E tests.

2. **Data Loss & Conflict Reconciliation Report (`scratch/doc_reconciliation_report.md`)**:
   - Itemized audit of 52 entries across 5 source files.
   - Explicit tags and remarks explaining data consolidation, conflict resolution (e.g. inline error banners eliminated for Two-Section Quality Report card; marker red color override eliminated for layer color preservation), and retirement of historical commit-level release notes.
   - Verified **zero accidental data loss** of operational invariants or domain knowledge.
   - Collaborative review conducted; user reviewed and formally approved the report dispositions and resolutions.

3. **Agentic System & Skill Alignment**:
   - Updated `GEMINI.md` to reference `docs/webapp/design_system.md` and removed obsolete `webapp/CHANGELOG.md` references from development lifecycle.
   - Updated `.agents/rules/webapp-design.md` to cross-link `docs/webapp/design_system.md`, `archetypes.md`, and `actionable_inventory.md`.
   - Updated `.agents/skills/task-init/SKILL.md` to index `docs/webapp/` in Step 1.
   - Updated `.agents/skills/task-record/SKILL.md` to maintain `docs/webapp/` and `webapp/README.md`.
   - Updated `webapp/README.md` and `docs/architecture.md` to reference `docs/webapp/`.

---

## Verification & Quality Gates

- **Relative Link Integrity**: Python scanner verified 100% of relative markdown links in `docs/`, `.agents/`, and `webapp/` resolve to existing targets.
- **Portability Invariant**: Verified zero host-specific absolute paths (`file:///`, `D:\ProiecteRealizate`) across repository files.
- **Vitest Unit Suite**: 80/80 passed in 391ms (`npm --prefix webapp test`).
- **Production Bundle**: `npm --prefix webapp run build` compiled cleanly in 522ms with zero TypeScript errors.
- **Ruff Linter**: `ruff check pytransdatro api tests tools` passed with all checks clean.
- **Core Geodetic Suite**: 54/54 tests passed in 3.33s (`tests/test_trans_ro.py`).
