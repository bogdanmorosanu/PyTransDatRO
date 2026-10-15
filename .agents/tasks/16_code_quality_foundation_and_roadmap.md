# Task 16: Code Quality Foundation, Rules & Multi-Chat Refactoring Roadmap

## Objectives
- [x] Establish a comprehensive agentic AI foundation for the upcoming code refactoring iteration.
- [x] Create and refine multi-layer agent rules (`code-quality-standards.md`, `configuration-standards.md`, `python-style.md`, `webapp-design.md`).
- [x] Implement reusable agent skills (`code-audit`, `regression-verify`) for automated inspection and benchmark validation.
- [x] Create the master 6-phase refactoring roadmap (`refactor_roadmap.md`) to structure future chat sessions.
- [x] Verify existing test suites and frontend builds to establish a baseline before refactoring begins.

---

## Architectural & Process Decisions

1. **Lightweight Tooling Strategy**:
   - Tailored specifically for a single-developer workflow: Ruff for Python linting and formatting; standard `tsc` and lightweight Vitest for frontend utilities.
   - Heavy pre-commit hooks and enterprise CI suites were intentionally excluded to eliminate development friction.

2. **Web UI Architecture & Styling Unification**:
   - Maintain pure Vanilla TypeScript without heavy reactive framework migrations (Vue/React).
   - Deconstruct the 3,500+ line monolithic `webapp/src/style/components.css` into category-based stylesheets (`forms.css`, `buttons.css`, `badges.css`, `layout.css`, `telemetry.css`).
   - Create reusable TS UI factory functions in `webapp/src/ui/` (`renderSegmentedControl()`, `renderSelect()`, `renderCopyBadge()`) to eliminate ad-hoc element styling and guarantee visual consistency across tabs.

3. **Configuration & Parameter Hierarchy**:
   - Immutable mathematical constants (WGS84 ellipsoid parameters, Stereo 70 projection origins, false coordinates) reside in `pytransdatro/constants.py`.
   - Grid-specific Helmert 2D parameters remain dynamic, extracted from the `.spg` binary header via `SpgReader`.
   - Operational parameters (ports, database paths, rate limits, log levels) are managed via `.env` and Pydantic `BaseSettings` (`api/config.py`).
   - Client parameters (map defaults, thresholds) are accessed via `import.meta.env` (`webapp/src/config.ts`).

4. **Targeted Frontend Testing**:
   - Introduce Vitest specifically for critical parsing and serialization logic (`dmsParser.ts`, `pointParser.ts`, `gridBounds.ts`, `exporters.ts` for DXF/CSV/GeoJSON), avoiding the high maintenance cost of mocking Leaflet map interactions.

5. **Phased Multi-Chat Execution**:
   - Divide the entire refactoring iteration into 6 self-contained chats to maintain optimal context window efficiency and prevent code degradation.

---

## Artifacts Created & Modified

### 1. Master Refactoring Roadmap
- [`refactor_roadmap.md`](../../refactor_roadmap.md): Complete multi-chat master plan detailing the goals, checklists, and acceptance criteria for Phases 1 through 6.

### 2. Agent Rule Files (`.agents/rules/`)
- [`code-quality-standards.md`](../rules/code-quality-standards.md): Cross-layer rules enforcing zero regression, standard library purity in core, layer decoupling, and fail-safe exception handling.
- [`configuration-standards.md`](../rules/configuration-standards.md): Categorization rules for invariant constants vs. dynamic `.env` configuration.
- [`python-style.md`](../rules/python-style.md): Updated with universal type annotations, mathematical docstring requirements, and prohibitions against silent `except: pass` blocks.
- [`webapp-design.md`](../rules/webapp-design.md): Updated with rules prohibiting ad-hoc element styling and mandating modular CSS architecture and UI helper factories.

### 3. Reusable Agent Skills (`.agents/skills/`)
- [`code-audit/SKILL.md`](../skills/code-audit/SKILL.md): Automated static analysis workflow (Ruff checks, prohibited anti-pattern scans, boundary checks, TypeScript build).
- [`regression-verify/SKILL.md`](../skills/regression-verify/SKILL.md): Automated test runner executing core transformations, API integration tests, and official ANCPI benchmark comparisons.
- [`doc-record/SKILL.md`](../skills/doc-record/SKILL.md): Living documentation manager mapping code subsystems to authoritative guides in `docs/`.

---

## Verification Results

1. **Core Transformation Test Suite**:
   - Executed: `D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/test_trans_ro.py -q`
   - Result: `70 passed in 1.10s` (100% pass rate).
2. **API & Telemetry Test Suite**:
   - Executed: `D:\anaconda3\envs\pytransdat\python.exe -m pytest tests/test_api_modern.py tests/test_api_legacy.py tests/test_api_rate_limiter.py tests/test_telemetry.py -q`
   - Result: `32 passed in 10.02s` (100% pass rate).
3. **Frontend TypeScript & Bundler Verification**:
   - Executed: `npm --prefix webapp run build`
   - Result: Clean build with zero TypeScript compiler errors (`tsc && vite build` completed successfully).
