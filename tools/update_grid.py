#!/usr/bin/env python3
"""
Unified SPG Grid Update & Lifecycle Orchestration Tool.

This tool provides a complete, unified CLI workflow for updating, inspecting,
benchmarking, and verifying Romanian geodetic transformation grid files (.spg)
released by ANCPI / Centrul Național de Cartografie (CNC).

Commands:
    status      Audits current grid deployment, metadata, and benchmark consistency
    inspect     Inspects and validates an incoming .spg file without modifying repo
    apply       Archives current grid, installs candidate .spg, and generates metadata
    benchmark   Runs full 9,046-point ANCPI replication benchmark and updates reports
    verify      Runs pytest and Vitest verification suites against active grid
    rollback    Restores previously archived grid from pytransdatro/grids/archive/
"""

import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# Add repository root to Python path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

try:
    import numpy as np
except ImportError:
    print("Error: numpy is required to inspect and process .spg files.")
    sys.exit(1)

from tools.export_spg_metadata_json import export_spg_metadata
from tools.generate_ancpi_comparison_reports import generate_reports

GRIDS_DIR = PROJECT_ROOT / "pytransdatro" / "grids"
ARCHIVE_DIR = GRIDS_DIR / "archive"
METADATA_JSON_PATH = GRIDS_DIR / "spg_metadata.json"
COMPARISON_JSON_PATH = GRIDS_DIR / "ancpi_comparison_summary.json"
TESTS_DATA_DIR = PROJECT_ROOT / "tests" / "data"
GRID_BOUNDS_TS = PROJECT_ROOT / "webapp" / "src" / "utils" / "gridBounds.ts"


def compute_sha256(file_path: Path) -> str:
    """Computes SHA-256 hash of a file."""
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def format_bytes(num_bytes: int) -> str:
    """Formats raw bytes into a human-readable size string."""
    for unit in ["B", "KB", "MB", "GB"]:
        if abs(num_bytes) < 1024.0:
            return f"{num_bytes:3.2f} {unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:.2f} TB"


def get_active_spg_files() -> List[Path]:
    """Returns list of all .spg files directly in pytransdatro/grids/."""
    if not GRIDS_DIR.exists():
        return []
    return [p for p in GRIDS_DIR.glob("*.spg") if p.is_file()]


def get_active_spg_file() -> Optional[Path]:
    """Returns single active .spg file or None if zero or multiple exist."""
    files = get_active_spg_files()
    return files[0] if len(files) == 1 else None


def load_spg_data(spg_path: Path) -> Dict[str, Any]:
    """Loads and unpickles an SPG grid file."""
    raw = np.load(str(spg_path), allow_pickle=True)
    if not isinstance(raw, dict) and hasattr(raw, "item"):
        return raw.item()
    if isinstance(raw, dict):
        return raw
    raise ValueError(f"Unrecognized SPG data format in {spg_path}: {type(raw)}")


# -----------------------------------------------------------------------------
# COMMAND: status
# -----------------------------------------------------------------------------
def cmd_status(_args: argparse.Namespace) -> int:
    """Audits current active grid deployment, metadata freshness, and benchmark consistency."""
    print("=" * 70)
    print("PyTransDatRO Grid Environment & Synchronization Audit")
    print("=" * 70)

    spg_files = get_active_spg_files()
    active_spg = get_active_spg_file()

    # 1. Grid File Status
    print("\n[1] Active Binary Grid (.spg):")
    if not spg_files:
        print("  [ERROR] No .spg grid file found in pytransdatro/grids/!")
    elif len(spg_files) > 1:
        print(f"  [ERROR] Multiple .spg files detected in {GRIDS_DIR}:")
        for f in spg_files:
            print(f"    - {f.name} ({format_bytes(f.stat().st_size)})")
        print("  AmbiguousGridError will occur. Use 'apply' or archive extra files.")
    else:
        active_f = spg_files[0]
        size = active_f.stat().st_size
        sha = compute_sha256(active_f)
        mtime = datetime.datetime.fromtimestamp(active_f.stat().st_mtime, tz=datetime.timezone.utc)
        print(f"  [OK] File Name     : {active_f.name}")
        print(f"       File Size     : {format_bytes(size)} ({size:,} bytes)")
        print(f"       SHA-256       : {sha}")
        print(f"       Last Modified : {mtime.isoformat()}")

    # 2. Metadata JSON Status
    print("\n[2] Pre-compiled Grid Metadata (spg_metadata.json):")
    if not METADATA_JSON_PATH.exists():
        print("  [WARN] spg_metadata.json does not exist. (API will attempt dynamic fallback)")
    else:
        try:
            with open(METADATA_JSON_PATH, "r", encoding="utf-8") as f:
                meta_json = json.load(f)
            file_info = meta_json.get("file_info", {})
            meta_name = file_info.get("name")
            meta_sha = file_info.get("sha256")
            gen_at = meta_json.get("_generated_at", "Unknown")

            if active_spg and meta_name == active_spg.name and meta_sha == sha:
                print(f"  [OK] Synchronized with active grid: {meta_name}")
                print(f"       Generated At : {gen_at}")
                release = meta_json.get("metadata", {}).get("release", {})
                print(
                    f"       Grid Release : {release.get('major')}.{release.get('minor')} (Rev {release.get('revision')})"
                )
            else:
                print("  [WARN] spg_metadata.json is OUT OF SYNC with active grid:")
                print(f"         JSON points to : {meta_name} (SHA: {meta_sha[:16]}...)")
                print(
                    f"         Active grid is : {active_spg.name if active_spg else 'None'} (SHA: {sha[:16]}...)"
                )
                print("         Run 'update_grid.py apply' to regenerate.")
        except Exception as e:
            print(f"  [ERROR] Failed to read spg_metadata.json: {e}")

    # 3. ANCPI Comparison Summary Status
    print("\n[3] ANCPI Benchmark Replication Summary (ancpi_comparison_summary.json):")
    if not COMPARISON_JSON_PATH.exists():
        print(
            "  [WARN] ancpi_comparison_summary.json does not exist. Run 'update_grid.py benchmark'."
        )
    else:
        try:
            with open(COMPARISON_JSON_PATH, "r", encoding="utf-8") as f:
                comp_json = json.load(f)
            bench_grid = comp_json.get("grid_file")
            bench_status = comp_json.get("overall_status")
            bench_pts = comp_json.get("total_points")
            bench_gen = comp_json.get("_generated_at")

            status_tag = "[OK]" if bench_status == "PASS" else "[FAIL]"
            print(f"  {status_tag} Benchmark Status : {bench_status} ({bench_pts:,} points)")
            print(f"       Benchmarked Grid : {bench_grid}")
            print(f"       Generated At     : {bench_gen}")

            if active_spg and bench_grid and bench_grid != active_spg.name:
                print(
                    f"  [WARN] Benchmark was computed with '{bench_grid}', but active grid is '{active_spg.name}'!"
                )
                print("         Run 'update_grid.py benchmark' after verifying expected test data.")
        except Exception as e:
            print(f"  [ERROR] Failed to read ancpi_comparison_summary.json: {e}")

    # 4. Test Reference Datasets
    print("\n[4] Official Test Datasets (tests/data/):")
    in_e2s = TESTS_DATA_DIR / "etrs89_to_st70_input.csv"
    exp_e2s = TESTS_DATA_DIR / "etrs89_to_st70_expected.csv"
    in_s2e = TESTS_DATA_DIR / "st70_to_etrs89_input.csv"
    exp_s2e = TESTS_DATA_DIR / "st70_to_etrs89_expected.csv"

    test_files = [in_e2s, exp_e2s, in_s2e, exp_s2e]
    all_exist = all(tf.exists() for tf in test_files)
    if all_exist:
        exp_mtime = max(exp_e2s.stat().st_mtime, exp_s2e.stat().st_mtime)
        exp_dt = datetime.datetime.fromtimestamp(exp_mtime, tz=datetime.timezone.utc)
        print(f"  [OK] All 4 reference test CSVs present ({exp_e2s.stat().st_size:,} bytes).")
        print(f"       Expected Data Timestamp: {exp_dt.isoformat()}")
        if active_spg and active_spg.stat().st_mtime > exp_mtime:
            print("  [NOTE] Active grid was modified AFTER expected test data.")
            print("         Verify whether tests/data/*_expected.csv reflect the active grid.")
    else:
        print("  [ERROR] Missing test CSV files in tests/data/!")

    # 5. Archived Grids
    if ARCHIVE_DIR.exists():
        archives = list(ARCHIVE_DIR.glob("*.spg"))
        if archives:
            print(f"\n[5] Archived Grid Backups ({len(archives)} in {ARCHIVE_DIR.name}/):")
            for af in sorted(archives, key=lambda p: p.stat().st_mtime, reverse=True)[:5]:
                print(f"    - {af.name} ({format_bytes(af.stat().st_size)})")

    print("\n" + "=" * 70)
    return 0


# -----------------------------------------------------------------------------
# COMMAND: inspect
# -----------------------------------------------------------------------------
def cmd_inspect(args: argparse.Namespace) -> int:
    """Inspects and validates an incoming .spg candidate file without touching the repo."""
    target_path = Path(args.spg_file).resolve()
    if not target_path.exists():
        print(f"[ERROR] Specified file not found: {target_path}")
        return 1
    if not target_path.name.endswith(".spg"):
        print(f"[WARN] File '{target_path.name}' does not have standard '.spg' extension.")

    print("=" * 70)
    print(f"Inspecting Candidate SPG Grid: {target_path.name}")
    print("=" * 70)

    try:
        data = load_spg_data(target_path)
    except Exception as e:
        print(f"[ERROR] Failed to unpickle / parse grid file: {e}")
        return 1

    file_size = target_path.stat().st_size
    sha256 = compute_sha256(target_path)

    metadata = data.get("metadata", {})
    params = data.get("params", {})
    grids = data.get("grids", {})

    print("\n[Basic File Properties]")
    print(f"  File Name : {target_path.name}")
    print(f"  File Size : {format_bytes(file_size)} ({file_size:,} bytes)")
    print(f"  SHA-256   : {sha256}")

    print("\n[Metadata]")
    print(f"  Creator     : {metadata.get('created_by', 'N/A')}")
    print(f"  License     : {metadata.get('license', 'N/A')}")
    release = metadata.get("release", {})
    print(
        f"  Release     : v{release.get('major', '?')}.{release.get('minor', '?')} (Rev {release.get('revision', '?')})"
    )
    print(f"  Valid From  : {metadata.get('valid_from', 'N/A')}")
    print(f"  Notes       : {metadata.get('notes', 'N/A')}")

    print("\n[Parameters & Helmert]")
    interp = params.get("interpolation", {})
    print(
        f"  Interpolation Strategy : Horizontal={interp.get('horizontal')} (Bicubic), Vertical={interp.get('vertical')} (Colocate)"
    )
    helmert = params.get("helmert", {})
    os_st70 = helmert.get("os_st70", {})
    st70_os = helmert.get("st70_os", {})
    print(
        f"  Helmert os_st70 : tE={os_st70.get('tE')}, tN={os_st70.get('tN')}, dm={os_st70.get('dm')}, Rz={os_st70.get('Rz')}"
    )
    print(
        f"  Helmert st70_os : tE={st70_os.get('tE')}, tN={st70_os.get('tN')}, dm={st70_os.get('dm')}, Rz={st70_os.get('Rz')}"
    )

    print("\n[Grid Layers]")
    if "geodetic_shifts" in grids:
        g2d = grids["geodetic_shifts"]
        grid_arr = g2d.get("grid")
        meta_2d = g2d.get("metadata", {})
        shape_str = str(grid_arr.shape) if hasattr(grid_arr, "shape") else "N/A"
        crs = meta_2d.get("crs_type") or g2d.get("crs", "N/A")
        mine = meta_2d.get("mine", 0.0)
        maxe = meta_2d.get("maxe", 0.0)
        minn = meta_2d.get("minn", 0.0)
        maxn = meta_2d.get("maxn", 0.0)
        stepe = meta_2d.get("stepe", 0.0)
        stepn = meta_2d.get("stepn", 0.0)
        print(f"  Layer 1 (2D Shifts) : shape={shape_str}, CRS={crs}")
        print(f"    Easting Extents   : {mine:,.3f} m -> {maxe:,.3f} m (step {stepe} m)")
        print(f"    Northing Extents  : {minn:,.3f} m -> {maxn:,.3f} m (step {stepn} m)")
    else:
        print("  [ERROR] Missing 'geodetic_shifts' layer in grids dictionary!")

    if "geoid_heights" in grids:
        g1d = grids["geoid_heights"]
        grid_arr = g1d.get("grid")
        meta_1d = g1d.get("metadata", {})
        shape_str = str(grid_arr.shape) if hasattr(grid_arr, "shape") else "N/A"
        crs = meta_1d.get("crs_type") or g1d.get("crs", "N/A")
        minla = meta_1d.get("minla", 0.0)
        maxla = meta_1d.get("maxla", 0.0)
        minphi = meta_1d.get("minphi", 0.0)
        maxphi = meta_1d.get("maxphi", 0.0)
        stepla = meta_1d.get("stepla", 0.0)
        stepphi = meta_1d.get("stepphi", 0.0)
        print(f"  Layer 2 (1D Geoid)  : shape={shape_str}, CRS={crs}")
        print(f"    Longitude Extents : {minla:.6f}° -> {maxla:.6f}° (step {stepla:.6f}°)")
        print(f"    Latitude Extents  : {minphi:.6f}° -> {maxphi:.6f}° (step {stepphi:.6f}°)")
    else:
        print("  [ERROR] Missing 'geoid_heights' layer in grids dictionary!")

    # Compare with current active grid if present
    active_spg = get_active_spg_file()
    if active_spg and active_spg.resolve() != target_path:
        print(f"\n[Comparison with Current Active Grid ({active_spg.name})]")
        try:
            curr_data = load_spg_data(active_spg)
            curr_rel = curr_data.get("metadata", {}).get("release", {})
            cand_rel = release
            print(f"  Active Version    : v{curr_rel.get('major')}.{curr_rel.get('minor')}")
            print(f"  Candidate Version : v{cand_rel.get('major')}.{cand_rel.get('minor')}")
            curr_size = active_spg.stat().st_size
            diff_size = file_size - curr_size
            print(f"  Size Difference   : {diff_size:+d} bytes")
        except Exception as e:
            print(f"  Could not load active grid for comparison: {e}")

    print("\n[Integrity Evaluation]")
    is_valid = ("geodetic_shifts" in grids) and ("geoid_heights" in grids) and ("helmert" in params)
    if is_valid:
        print("  [PASS] Candidate file meets all PyTransDatRO structural invariants.")
    else:
        print("  [FAIL] Candidate file is missing mandatory internal layers or parameters.")

    print("=" * 70)
    return 0 if is_valid else 1


# -----------------------------------------------------------------------------
# COMMAND: apply
# -----------------------------------------------------------------------------
def cmd_apply(args: argparse.Namespace) -> int:
    """Applies a new SPG grid: archives current grid, installs candidate, and exports metadata."""
    candidate_path = Path(args.spg_file).resolve()
    if not candidate_path.exists():
        print(f"[ERROR] Candidate grid file not found: {candidate_path}")
        return 1

    print("=" * 70)
    print(f"Applying Grid Update: {candidate_path.name}")
    print("=" * 70)

    # 1. Inspect candidate
    try:
        data = load_spg_data(candidate_path)
        if "geodetic_shifts" not in data.get("grids", {}) or "geoid_heights" not in data.get(
            "grids", {}
        ):
            print("[ERROR] Candidate grid failed validation check (missing layers). Aborting.")
            return 1
    except Exception as e:
        print(f"[ERROR] Candidate grid could not be parsed: {e}. Aborting.")
        return 1

    cand_sha = compute_sha256(candidate_path)
    active_spg = get_active_spg_file()

    # 2. Check if identical to active grid
    if active_spg and compute_sha256(active_spg) == cand_sha and not args.force:
        print(f"[INFO] Candidate grid is byte-for-byte IDENTICAL to active grid: {active_spg.name}")
        print("Regenerating metadata JSON to ensure synchronization...")
        export_spg_metadata(active_spg, METADATA_JSON_PATH)
        print("[OK] Metadata refreshed.")
        return 0

    # 3. Archive existing grid(s) in pytransdatro/grids/
    ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    existing_spgs = get_active_spg_files()

    for old_spg in existing_spgs:
        if old_spg.resolve() == candidate_path.resolve():
            continue
        dest_archive = ARCHIVE_DIR / old_spg.name
        if dest_archive.exists():
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            dest_archive = ARCHIVE_DIR / f"{old_spg.stem}_{ts}{old_spg.suffix}"
        print(f"Archiving previous grid '{old_spg.name}' -> '{dest_archive.name}'...")
        shutil.move(str(old_spg), str(dest_archive))

    # 4. Copy candidate into pytransdatro/grids/
    dest_path = GRIDS_DIR / candidate_path.name
    if candidate_path.resolve() != dest_path.resolve():
        print(f"Installing new grid to '{dest_path}'...")
        shutil.copy2(str(candidate_path), str(dest_path))
    else:
        print(f"Candidate grid is already positioned at '{dest_path}'.")

    # 5. Extract fresh metadata JSON
    print(f"Extracting metadata and layer statistics -> '{METADATA_JSON_PATH.name}'...")
    export_spg_metadata(dest_path, METADATA_JSON_PATH)

    # 6. Check client bounds against webapp/src/utils/gridBounds.ts
    _check_client_bounds(data)

    # 7. Print high-visibility operator reminder
    print("\n" + "*" * 70)
    print("SUCCESS: New grid has been successfully installed into pytransdatro/grids/!")
    print("*" * 70)
    print("MANDATORY NEXT STEPS FOR FULL VERIFICATION:")
    print("1. Update Expected Benchmark Test Data:")
    print("   Transformation values will change with the new grid. To keep benchmark")
    print("   reports and regression tests passing, run the official ANCPI desktop application")
    print("   (TransDatRO or RomGEO) using this new grid on:")
    print("     - tests/data/etrs89_to_st70_input.csv  ->  tests/data/etrs89_to_st70_expected.csv")
    print("     - tests/data/st70_to_etrs89_input.csv  ->  tests/data/st70_to_etrs89_expected.csv")
    print("2. Regenerate Comparison Reports & ancpi_comparison_summary.json:")
    print("   Run:")
    print("     python tools/update_grid.py benchmark")
    print("3. Restart FastAPI & Vite Development Services:")
    print("   The Python backend caches the grid singleton in RAM. Restart the service:")
    print("     python -m uvicorn api.main:app --port 8000")
    print("*" * 70 + "\n")
    return 0


def _check_client_bounds(data: Dict[str, Any]) -> None:
    """Advises if new grid coverage deviates significantly from webapp client envelope."""
    try:
        g2d_meta = data.get("grids", {}).get("geodetic_shifts", {}).get("metadata", {})
        mine = float(g2d_meta.get("mine", 0.0))
        maxe = float(g2d_meta.get("maxe", 0.0))
        minn = float(g2d_meta.get("minn", 0.0))
        maxn = float(g2d_meta.get("maxn", 0.0))

        # Check Stereo70 envelope
        # Standard: minN: 213600, maxN: 785640, minE: 101260, maxE: 890790
        if minn < 210000 or maxn > 790000 or mine < 100000 or maxe > 895000:
            print(
                "[ADVISORY] Grid coverage extends beyond standard Stereo 70 bounds in gridBounds.ts."
            )
            print(f"           New Extents: N=[{minn}, {maxn}], E=[{mine}, {maxe}]")
            print(
                f"           Review {GRID_BOUNDS_TS.name} if client envelope preflight checks fail."
            )
    except Exception:
        pass


# -----------------------------------------------------------------------------
# COMMAND: benchmark
# -----------------------------------------------------------------------------
def cmd_benchmark(args: argparse.Namespace) -> int:
    """Executes the full 9,046-point ANCPI comparison benchmark and exports reports."""
    print("=" * 70)
    print("Executing ANCPI Replication Benchmark (9,046 Points)")
    print("=" * 70)

    active_spg = get_active_spg_file()
    if not active_spg:
        print("[ERROR] Exactly one active .spg grid file required in pytransdatro/grids/.")
        return 1

    print(f"Active Grid File : {active_spg.name}")
    print(f"Testing datasets in {TESTS_DATA_DIR}...")

    # Run the report generator
    generate_reports()

    # Read the generated summary
    if not COMPARISON_JSON_PATH.exists():
        print("[ERROR] Comparison JSON was not generated!")
        return 1

    with open(COMPARISON_JSON_PATH, "r", encoding="utf-8") as f:
        comp = json.load(f)

    status = comp.get("overall_status")
    print(f"\nBenchmark Overall Result: {status}")

    s2e = comp.get("st70_to_etrs89", {})
    e2s = comp.get("etrs89_to_st70", {})

    print(f"  Stereo70 -> ETRS89 : Status={s2e.get('status')}")
    print(f"    Max Abs Error    : {s2e.get('max_abs_error')}")
    print(f"    Mean Abs Error   : {s2e.get('mean_abs_error')}")

    print(f"  ETRS89 -> Stereo70 : Status={e2s.get('status')}")
    print(f"    Max Abs Error    : {e2s.get('max_abs_error')}")
    print(f"    Mean Abs Error   : {e2s.get('mean_abs_error')}")

    if status != "PASS":
        print("\n" + "!" * 70)
        print("[FAIL] Benchmark residuals EXCEEDED allowable ANCPI tolerances!")
        print("This typically indicates that tests/data/*_expected.csv were computed")
        print("using an older grid release. Update the expected CSVs with official ANCPI")
        print("RomGEO / TransDatRO output for this new grid.")
        print("!" * 70)
        if not getattr(args, "allow_fail", False):
            return 1

    print("\n[OK] Benchmark reports and ancpi_comparison_summary.json updated successfully.")
    return 0


# -----------------------------------------------------------------------------
# COMMAND: verify
# -----------------------------------------------------------------------------
def cmd_verify(_args: argparse.Namespace) -> int:
    """Runs core pytest and Vitest verification suites."""
    print("=" * 70)
    print("Running PyTransDatRO Verification Test Suites")
    print("=" * 70)

    python_exe = sys.executable

    # 1. Pytest
    print("\n[1/2] Running Python core geodetic & API tests (pytest)...")
    pytest_cmd = [
        python_exe,
        "-m",
        "pytest",
        "tests/test_trans_ro.py",
        "tests/test_api_modern.py",
        "-v",
    ]
    res_py = subprocess.run(pytest_cmd, cwd=str(PROJECT_ROOT))
    if res_py.returncode != 0:
        print("[ERROR] Pytest execution failed!")
        return res_py.returncode

    # 2. Vitest
    print("\n[2/2] Running frontend unit tests (Vitest)...")
    npm_cmd = "npm --prefix webapp test"
    res_npm = subprocess.run(npm_cmd, shell=True, cwd=str(PROJECT_ROOT))
    if res_npm.returncode != 0:
        print("[ERROR] Vitest execution failed!")
        return res_npm.returncode

    print("\n" + "=" * 70)
    print("[PASS] All verification suites passed with zero errors!")
    print("=" * 70)
    return 0


# -----------------------------------------------------------------------------
# COMMAND: rollback
# -----------------------------------------------------------------------------
def cmd_rollback(args: argparse.Namespace) -> int:
    """Rolls back to the most recently archived grid in pytransdatro/grids/archive/."""
    print("=" * 70)
    print("Grid Rollback Procedure")
    print("=" * 70)

    if not ARCHIVE_DIR.exists():
        print(f"[ERROR] Archive directory '{ARCHIVE_DIR}' does not exist. Nothing to restore.")
        return 1

    archives = sorted(ARCHIVE_DIR.glob("*.spg"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not archives:
        print(f"[ERROR] No archived .spg files found in '{ARCHIVE_DIR}'.")
        return 1

    target_archive = archives[0]
    print(f"Target archived grid to restore: '{target_archive.name}'")

    current_spgs = get_active_spg_files()
    for cur in current_spgs:
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_name = f"{cur.stem}_rollback_backup_{ts}{cur.suffix}"
        backup_dest = ARCHIVE_DIR / backup_name
        print(f"Moving current grid '{cur.name}' -> '{backup_name}'...")
        shutil.move(str(cur), str(backup_dest))

    restored_dest = GRIDS_DIR / target_archive.name
    print(f"Restoring '{target_archive.name}' to active grids directory...")
    shutil.move(str(target_archive), str(restored_dest))

    print("Regenerating metadata JSON for restored grid...")
    export_spg_metadata(restored_dest, METADATA_JSON_PATH)

    print("\n[OK] Rollback complete! Active grid restored to:", restored_dest.name)
    print("Please run 'python tools/update_grid.py benchmark' to re-verify benchmark data.")
    return 0


# -----------------------------------------------------------------------------
# CLI Entrypoint
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Unified SPG Grid Management & Lifecycle Tool for PyTransDatRO",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # status
    p_status = subparsers.add_parser("status", help="Audit current grid deployment & metadata")
    p_status.set_defaults(func=cmd_status)

    # inspect
    p_inspect = subparsers.add_parser(
        "inspect", help="Inspect candidate .spg file without modifying repo"
    )
    p_inspect.add_argument("spg_file", help="Path to candidate .spg file to inspect")
    p_inspect.set_defaults(func=cmd_inspect)

    # apply
    p_apply = subparsers.add_parser(
        "apply", help="Archive current grid, install candidate .spg, export metadata"
    )
    p_apply.add_argument("spg_file", help="Path to new .spg file to install")
    p_apply.add_argument(
        "--force", action="store_true", help="Force application even if file is identical"
    )
    p_apply.set_defaults(func=cmd_apply)

    # benchmark
    p_bench = subparsers.add_parser("benchmark", help="Run 9,046-point ANCPI replication benchmark")
    p_bench.add_argument(
        "--allow-fail", action="store_true", help="Do not exit with error code if benchmark fails"
    )
    p_bench.set_defaults(func=cmd_benchmark)

    # verify
    p_verify = subparsers.add_parser("verify", help="Run pytest and Vitest verification suites")
    p_verify.set_defaults(func=cmd_verify)

    # rollback
    p_rollback = subparsers.add_parser("rollback", help="Restore previously archived grid")
    p_rollback.set_defaults(func=cmd_rollback)

    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
