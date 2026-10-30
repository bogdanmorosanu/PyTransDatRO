# Contributing to PyTransDatRO

Thank you for your interest in contributing to **PyTransDatRO**! As an open-source geodetic engineering project licensed under the [MIT License](LICENSE), we welcome bug reports, suggestions, improvements, and code contributions.

---

## 1. Architectural Invariants

Before writing code, please review the immutable core invariants that govern this repository:

1. **Pure Python Core Geodetic Library (`pytransdatro/`)**:
   - The core mathematical and grid transformation engine is implemented in pure Python (standard library only).
   - No external scientific C/compiled dependencies (such as NumPy or SciPy) are permitted for coordinate computation.
   - **Exception**: `numpy` is permitted *strictly* for deserializing binary grid files (`.spg`).
   - The core library must **never** import from the web API layer (`api/`) or web frameworks.
2. **Zero Functional Regression Rule**:
   - All coordinate transformations must remain byte-for-byte identical to official ANCPI / RomGEO benchmarks with $< 10^{-9}$ numerical tolerance across all 9,046 national benchmark control points.
   - Always run the benchmark suite before submitting changes:
     ```bash
     python tools/generate_ancpi_comparison_reports.py
     ```
3. **Decoupled Three-Tier Architecture**:
   - Geodetic Math & Grids: `pytransdatro/` (pure Python).
   - REST API & Endpoints: `api/` (FastAPI, Pydantic, Uvicorn).
   - Interactive Studio & Maps: `webapp/` (TypeScript, Leaflet, Vite).
4. **Machine-Agnostic Portability**:
   - Never commit machine-specific absolute file paths (e.g. `D:\...` or `C:\Users\...`). All links and references must be repository-relative.

---

## 2. Development Setup

### Python Environment
1. Clone the repository:
   ```bash
   git clone https://github.com/bogdanmorosanu/PyTransDatRO.git
   cd PyTransDatRO
   ```
2. Create and activate a virtual environment (Python 3.10+):
   ```bash
   python -m venv .venv
   # Windows:
   .venv\Scripts\activate
   # Linux/macOS:
   source .venv/bin/activate
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt -r requirements-api.txt -r requirements-dev.txt
   ```

### Web Application Environment
1. Navigate to `webapp/` and install dependencies:
   ```bash
   npm --prefix webapp ci
   ```
2. Build the frontend or run in development mode:
   ```bash
   npm --prefix webapp run dev
   ```

---

## 3. Code Standards & Verification

Every pull request must pass the full automated quality suite:

### 1. Code Linting & Formatting
```bash
ruff check pytransdatro api tests tools
ruff format --check pytransdatro api tests tools
```

### 2. Python Tests & ANCPI Benchmarks
```bash
python tools/generate_ancpi_comparison_reports.py
pytest -m "not legacy_live" -q
```

### 3. Frontend TypeScript & Unit Tests
```bash
npm --prefix webapp test -- --run
npm --prefix webapp run build
```

### 4. End-to-End Browser Tests (Playwright)
```bash
# Windows:
.\test-e2e.ps1
# Or via npm:
npm --prefix webapp run test:e2e
```

### 5. Version Consistency
```bash
python tools/bump_version.py --check
```

---

## 4. Submitting Changes

1. **Branching**: Create a focused topic branch from `main`:
   ```bash
   git checkout -b feat/your-feature-name
   # or
   git checkout -b fix/your-bugfix-name
   ```
2. **Conventional Commits**: Write clear, descriptive commit messages:
   - `feat(...)`: New feature or capability
   - `fix(...)`: Bug fix or numerical correction
   - `refactor(...)`: Code cleanup without functional change
   - `test(...)`: Adding or updating test cases
   - `docs(...)`: Documentation updates
3. **Open a Pull Request**: Push your branch and open a PR against `main`. Provide a clear description of your changes and reference any related issues.
