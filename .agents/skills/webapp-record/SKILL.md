---
name: webapp-record
description: Validates and documents changes made to the pytransdat web application and web map. Updates webapp/CHANGELOG.md with concise categorized entries (Added, Changed, Fixed), refreshes webapp/README.md to reflect the latest state/architecture, and verifies the build via npm run build. Use whenever completing an edit, feature, or bugfix in webapp/.
---

# Webapp Record Skill

## Purpose
This skill captures and documents the changes made to the **PyTransDatRO Web Application** (`webapp/`) during the current conversation or work session. It ensures that:
1. The chronological development history in [`webapp/CHANGELOG.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/CHANGELOG.md) is updated with clear, concise bullet points.
2. The single-source-of-truth status in [`webapp/README.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/README.md) is updated if components, file structure, or architecture evolved.
3. Code integrity is verified with a clean TypeScript compile and production bundle build (`npm run build`).

---

## Operating Principles

1. **Concise, High-Signal Changelog**:
   - Focus on *what* was added, changed, or fixed and *why*.
   - Avoid conversational commentary or debugging intermediate steps.
   - Use standard Keep a Changelog categories: `### Added`, `### Changed`, `### Fixed`, `### Removed`.
2. **Synchronize Single Source of Truth**:
   - `webapp/README.md` must always reflect the *current, latest state* of `webapp/`.
   - Update the file tree in `webapp/README.md` Section 2 if files were created, moved, or deleted.
   - Update architecture or guardrail notes in Section 3 if behaviors changed (e.g., new map layer, changed thresholds, new API client methods).
3. **Build Verification Before Recording**:
   - Always run `npm run build` in `webapp/` before finalizing the documentation to guarantee no type errors or bundle issues were introduced.

---

## Workflow Steps

### Step 1: Inspect Changes via Git

Run git commands to identify all touched files and changes:
```powershell
git status -s webapp
git diff --stat webapp
```
Also inspect if any backend API routes (`api/routes/`, `api/schemas/`) were modified in tandem with the webapp.

### Step 2: Verify Build Integrity

In `webapp/`, verify that TypeScript compiles cleanly and Vite builds the bundle without errors:
```powershell
cd webapp
npm run build
```
If errors are reported:
- Resolve any TypeScript type errors or broken imports before recording the changes.

### Step 3: Update `webapp/CHANGELOG.md`

1. Open [`webapp/CHANGELOG.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/CHANGELOG.md).
2. Under the `## [Unreleased]` section (or a new dated version heading `## [X.Y.Z] - YYYY-MM-DD` if tagging a release), add concise bullet points:
   - **`### Added`**: New map layers, UI widgets, tools, export options, or API integration endpoints.
   - **`### Changed`**: UI layout adjustments, performance optimizations, updated clustering thresholds, or CSS token changes.
   - **`### Fixed`**: Bug fixes (e.g., coordinate reversal bugs, popup offset issues, leaflet container resize glitches).
   - **`### Removed`**: Deprecated code or unused assets.
3. Use `replace_file_content` to keep edits clean and atomic.

### Step 4: Update `webapp/README.md` (If Applicable)

Review [`webapp/README.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/README.md):
1. **Directory Structure**: If new components, styles, or utils were introduced or deleted, update the ASCII directory tree in Section 2.
2. **Key Architecture & Design Invariants**: If new map features (e.g. raster shift layers), guardrails, or API interactions were added, update Section 3.
3. Keep the documentation accurate, concise, and focused on the latest state.

### Step 5: Report Summary & Next Step

Provide a clear confirmation to the user:
- Summary of changes recorded in [`webapp/CHANGELOG.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/CHANGELOG.md).
- Any updates made to [`webapp/README.md`](file:///D:/ProiecteRealizate/pyTransDatRO/repo/PyTransDatRO/webapp/README.md).
- Status of `npm run build` verification.
- Proactively suggest running `task-commit` if the user is ready to commit the changes to Git.
