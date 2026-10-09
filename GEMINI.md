# PyTransDatRO Project Invariants

- **Pure Python Core Geodetic Library (`pytransdatro/`)**:
  - The core geodetic library is implemented strictly in pure Python (standard library only), with no external scientific C/compiled dependencies (e.g., NumPy, SciPy) for processing.
  - **Exception: `numpy` is allowed STRICTLY for reading and unpickling the binary grid files.** All data must be processed using pure Python.
  - Core data processing and mathematical pipelines in `pytransdatro/` must never import from or depend on the web API layer (`api/`) or web frameworks.

- **Web API Layer (`api/`)**:
  - The web service layer is built on top of `pytransdatro` using standard modern web libraries (`fastapi`, `uvicorn`, `pydantic`, `python-multipart`, `httpx`).
  - These web dependencies are optional extras (`pytransdatro[api]`) and must never be imported inside `pytransdatro/`.

- **Binary Grid Files**: Do not modify the loading structure or content of the binary grid files in `grids/`.

- **Python Environment**: Always use the `pytransdat` conda environment located at `D:\anaconda3\envs\pytransdat` (Python executable: `D:\anaconda3\envs\pytransdat\python.exe`) for running all commands, tests, and scripts in this project.

For a detailed explanation of the transformation pipeline and class structure, see [architecture.md](docs/architecture.md).

