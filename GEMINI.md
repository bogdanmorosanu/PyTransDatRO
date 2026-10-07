# PyTransDatRO Project Invariants

- **Pure Python Implementation**: The project is implemented using pure Python only (standard library), with no external scientific C/compiled dependencies (e.g., NumPy, SciPy) for processing. **Exception: `numpy` is allowed STRICTLY for reading and unpickling the binary grid files.** All data must be processed using pure Python.
- **Binary Grid Files**: Do not modify the loading structure or content of the binary grid files in `grids/`.
- **Python Environment**: Always use the `pytransdat` conda environment located at `D:\anaconda3\envs\pytransdat` (Python executable: `D:\anaconda3\envs\pytransdat\python.exe`) for running all commands, tests, and scripts in this project.

For a detailed explanation of the transformation pipeline and class structure, see [architecture.md](docs/architecture.md).
