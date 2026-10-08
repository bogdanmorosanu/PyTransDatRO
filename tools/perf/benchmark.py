"""PyTransDatRO Benchmarking and Performance Tracking Harness

Features:
  - Measures Cold Startup / Initialization latency.
  - Tests 4 scenarios (S1 Micro, S2 Small, S3 Medium, S4 Large) in:
      * Both directions: Stereo70 -> ETRS89 and ETRS89 -> Stereo70.
      * Both dimensions: 2D (pure planimetric) and 3D (with 1D geoid height).
      * Both cache conditions: Cold Cache (cleared) and Warm Cache (steady-state).
  - Captures high-resolution timing (time.perf_counter), points/sec throughput,
    per-point latency (us/pt), and cache hits/misses.
  - Appends results into a structured history log file (JSON and Markdown table)
    so every iteration/improvement is tracked over time against the baseline.
"""

import os
import sys
import csv
import time
import json
import platform
import argparse
from datetime import datetime
from pathlib import Path

# Add project root to sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pytransdatro

PERF_DIR = REPO_ROOT / "tools" / "perf"
DATA_DIR = PERF_DIR / "data"
LOG_JSON_FILE = PERF_DIR / "perf_history.json"
LOG_MD_FILE = PERF_DIR / "perf_history.md"

SCENARIO_FILES = {
    "S1_micro": ("s1_micro_clustered.csv", 10),
    "S2_small": ("s2_small_clustered.csv", 250),
    "S3_medium": ("s3_medium_corridor.csv", 2500),
    "S4_large": ("s4_large_dispersed.csv", 25000),
}


def load_dataset(filename):
    """Loads CSV points into lists of tuples: (pid, n, e, z)."""
    filepath = DATA_DIR / filename
    pts = []
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        for row in reader:
            pts.append((row[0], float(row[1]), float(row[2]), float(row[3])))
    return pts


def precompute_etrs89_coords(t, st70_pts):
    """Precomputes ETRS89 (lat_rad, lon_rad, z) for running reverse transformation tests."""
    etrs89_pts = []
    for pid, n, e, z in st70_pts:
        lat, lon, h = t.st70_to_etrs89(n, e, z)
        etrs89_pts.append((pid, lat, lon, h))
    return etrs89_pts


def clear_trans_cache(t):
    """Clears all LRU caches inside TransRO instance."""
    if hasattr(t._t_gr2d, "_init_interp"):
        t._t_gr2d._init_interp.cache_clear()
    if hasattr(t._t_gr2d, "get_vs_at_idxs_cached"):
        t._t_gr2d.get_vs_at_idxs_cached.cache_clear()
    if hasattr(t._t_gr1d, "get_vs_at_idxs_cached"):
        t._t_gr1d.get_vs_at_idxs_cached.cache_clear()


def get_cache_stats(t):
    """Returns dictionary of cache statistics."""
    stats = {}
    if hasattr(t._t_gr2d, "_init_interp"):
        info = t._t_gr2d._init_interp.cache_info()
        stats["gr2d_interp_hits"] = info.hits
        stats["gr2d_interp_misses"] = info.misses
    if hasattr(t._t_gr2d, "get_vs_at_idxs_cached"):
        info = t._t_gr2d.get_vs_at_idxs_cached.cache_info()
        stats["gr2d_vs_hits"] = info.hits
        stats["gr2d_vs_misses"] = info.misses
    return stats


def measure_initialization():
    """Measures cold instantiation time of TransRO (including SPG file deserialization)."""
    # Reset singleton reader to measure authentic cold start if desired
    pytransdatro.spg_reader.SpgReader._instance = None
    start = time.perf_counter()
    t = pytransdatro.TransRO()
    elapsed_ms = (time.perf_counter() - start) * 1000.0
    return t, elapsed_ms


def run_batch_st70_to_etrs89(t, pts, with_z=True):
    """Executes st70_to_etrs89 for a batch of points."""
    n_arr = [pt[1] for pt in pts]
    e_arr = [pt[2] for pt in pts]
    if with_z:
        z_arr = [pt[3] for pt in pts]
        t.st70_to_etrs89(n_arr, e_arr, z_arr)
    else:
        t.st70_to_etrs89(n_arr, e_arr, None)


def run_batch_etrs89_to_st70(t, pts, with_z=True):
    """Executes etrs89_to_st70 for a batch of points."""
    lat_arr = [pt[1] for pt in pts]
    lon_arr = [pt[2] for pt in pts]
    if with_z:
        h_arr = [pt[3] for pt in pts]
        t.etrs89_to_st70(lat_arr, lon_arr, h_arr)
    else:
        t.etrs89_to_st70(lat_arr, lon_arr, None)


def benchmark_single_case(t, runner_fn, pts, with_z, runs=3, cold_cache=False):
    """Runs a test case multiple times and computes min, mean, throughput."""
    timings = []
    for _ in range(runs):
        if cold_cache:
            clear_trans_cache(t)
        start = time.perf_counter()
        runner_fn(t, pts, with_z=with_z)
        elapsed = time.perf_counter() - start
        timings.append(elapsed)

    best_time = min(timings)
    avg_time = sum(timings) / len(timings)
    pt_count = len(pts)
    throughput = pt_count / best_time if best_time > 0 else 0
    latency_us = (best_time / pt_count) * 1e6 if pt_count > 0 else 0

    return {
        "pt_count": pt_count,
        "best_sec": round(best_time, 5),
        "avg_sec": round(avg_time, 5),
        "throughput_pts_sec": round(throughput, 1),
        "latency_us_per_pt": round(latency_us, 2),
    }


def update_history_files(run_record):
    """Appends benchmark run record to perf_history.json and updates perf_history.md."""
    # 1. Update JSON
    history = []
    if LOG_JSON_FILE.exists():
        try:
            with open(LOG_JSON_FILE, "r", encoding="utf-8") as f:
                history = json.load(f)
        except Exception:
            history = []

    history.append(run_record)
    with open(LOG_JSON_FILE, "w", encoding="utf-8") as f:
        json.dump(history, f, indent=2)

    # 2. Update Markdown
    md_content = generate_markdown_report(history)
    with open(LOG_MD_FILE, "w", encoding="utf-8") as f:
        f.write(md_content)


def generate_markdown_report(history):
    """Generates a comprehensive markdown log and comparative delta table."""
    md = []
    md.append("# PyTransDatRO Performance & Optimization History Log\n")
    md.append("This document tracks performance metrics over time across versions and optimization experiments.\n")
    md.append(f"*Last Updated:* `{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}`\n\n")

    md.append("## 1. Run Registry\n\n")
    md.append("| Run ID | Version / Tag | Timestamp | Description | Init (ms) |\n")
    md.append("|---|---|---|---|---|\n")
    for r in history:
        md.append(
            f"| `{r['run_id']}` | **{r['version_tag']}** | {r['timestamp'][:19]} | {r['description']} | {r['init_time_ms']:.2f} |\n"
        )
    md.append("\n---\n")

    md.append("## 2. Latest Benchmark Summary (`" + history[-1]["version_tag"] + "`)\n\n")
    latest = history[-1]
    md.append(f"**Description:** {latest['description']}  \n")
    md.append(f"**Python:** {latest['env']['python_version']} ({latest['env']['platform']})  \n")
    md.append(f"**Cold Initialization:** `{latest['init_time_ms']:.2f} ms`\n\n")

    md.append("### Detailed Scenario Metrics\n\n")
    md.append("| Scenario | Points | Direction | Mode | Cache | Best Time (s) | Throughput (pts/s) | Latency (µs/pt) |\n")
    md.append("|---|---|---|---|---|---|---|---|\n")

    for k, res in latest["results"].items():
        md.append(
            f"| {res['scenario']} | {res['points']} | {res['direction']} | {res['mode']} | {res['cache']} | "
            f"`{res['best_sec']:.4f}` | **{res['throughput_pts_sec']:,.0f}** | {res['latency_us_per_pt']:.1f} |\n"
        )

    # 3. Comparative Delta Table against Baseline (Run 0)
    if len(history) > 1:
        baseline = history[0]
        md.append("\n---\n")
        md.append(f"## 3. Comparative Evolution vs Baseline (`{baseline['version_tag']}`)\n\n")
        md.append("| Scenario & Mode | Direction | Cache | Baseline Throughput | Latest Throughput | Speedup (Delta %) |\n")
        md.append("|---|---|---|---|---|---|\n")

        for k, curr_res in latest["results"].items():
            if k in baseline["results"]:
                base_res = baseline["results"][k]
                base_t = base_res["throughput_pts_sec"]
                curr_t = curr_res["throughput_pts_sec"]
                pct = ((curr_t - base_t) / base_t) * 100.0 if base_t > 0 else 0
                sign = "+" if pct >= 0 else ""
                md.append(
                    f"| {curr_res['scenario']} ({curr_res['mode']}) | {curr_res['direction']} | {curr_res['cache']} | "
                    f"{base_t:,.0f} pts/s | {curr_t:,.0f} pts/s | **{sign}{pct:.1f}%** |\n"
                )

    return "".join(md)


def run_benchmark(tag="baseline", description="Initial unoptimized baseline", runs=3):
    print("=" * 70)
    print(f" PyTransDatRO Benchmark Harness - Tag: [{tag}]")
    print("=" * 70)

    # 1. Cold start
    print("\n[Step 1/3] Measuring Cold Initialization (SPG load + instantiation)...")
    t, init_ms = measure_initialization()
    print(f"  Cold Init Time: {init_ms:.2f} ms")

    # 2. Loading datasets
    print("\n[Step 2/3] Loading scenario datasets and preparing reverse coordinates...")
    datasets = {}
    for sc_key, (csv_name, expected_pts) in SCENARIO_FILES.items():
        st70_pts = load_dataset(csv_name)
        etrs89_pts = precompute_etrs89_coords(t, st70_pts)
        datasets[sc_key] = {"st70": st70_pts, "etrs89": etrs89_pts}
        print(f"  {sc_key}: Loaded {len(st70_pts)} points.")

    # 3. Running benchmark matrix
    print(f"\n[Step 3/3] Running Benchmark Matrix (Runs per test: {runs})...\n")
    results = {}

    test_matrix = [
        # (Scenario, Direction, Mode, CacheCondition)
        # S1 Micro (10 pts)
        ("S1_micro", "st70_to_etrs89", "3D", "Warm"),
        ("S1_micro", "etrs89_to_st70", "3D", "Warm"),
        # S2 Small (250 pts)
        ("S2_small", "st70_to_etrs89", "3D", "Cold"),
        ("S2_small", "st70_to_etrs89", "3D", "Warm"),
        ("S2_small", "etrs89_to_st70", "3D", "Cold"),
        ("S2_small", "etrs89_to_st70", "3D", "Warm"),
        # S3 Medium (2,500 pts)
        ("S3_medium", "st70_to_etrs89", "3D", "Cold"),
        ("S3_medium", "st70_to_etrs89", "3D", "Warm"),
        ("S3_medium", "etrs89_to_st70", "3D", "Cold"),
        ("S3_medium", "etrs89_to_st70", "3D", "Warm"),
        # S4 Large (25,000 pts)
        ("S4_large", "st70_to_etrs89", "2D", "Warm"),
        ("S4_large", "st70_to_etrs89", "3D", "Cold"),
        ("S4_large", "st70_to_etrs89", "3D", "Warm"),
        ("S4_large", "etrs89_to_st70", "2D", "Warm"),
        ("S4_large", "etrs89_to_st70", "3D", "Cold"),
        ("S4_large", "etrs89_to_st70", "3D", "Warm"),
    ]

    for sc_key, direction, mode, cache_state in test_matrix:
        with_z = (mode == "3D")
        is_cold = (cache_state == "Cold")
        pts = datasets[sc_key]["st70"] if direction == "st70_to_etrs89" else datasets[sc_key]["etrs89"]
        runner = run_batch_st70_to_etrs89 if direction == "st70_to_etrs89" else run_batch_etrs89_to_st70

        metrics = benchmark_single_case(t, runner, pts, with_z=with_z, runs=runs, cold_cache=is_cold)

        test_key = f"{sc_key}__{direction}__{mode}__{cache_state}"
        results[test_key] = {
            "scenario": sc_key,
            "points": metrics["pt_count"],
            "direction": direction,
            "mode": mode,
            "cache": cache_state,
            **metrics,
        }

        dir_label = "ST70->ETRS89" if direction == "st70_to_etrs89" else "ETRS89->ST70"
        print(
            f"  [{sc_key:<9}] {dir_label:<13} | {mode} | {cache_state:<4} Cache: "
            f"{metrics['best_sec']:7.4f}s | {metrics['throughput_pts_sec']:8,.0f} pts/s | "
            f"{metrics['latency_us_per_pt']:6.1f} us/pt"
        )

    # 4. Record Run
    run_record = {
        "run_id": f"run_{int(time.time())}",
        "version_tag": tag,
        "description": description,
        "timestamp": datetime.now().isoformat(),
        "init_time_ms": round(init_ms, 2),
        "env": {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "processor": platform.processor(),
        },
        "results": results,
    }

    update_history_files(run_record)
    print("\n" + "=" * 70)
    print(f" Benchmark complete! Log updated in:")
    print(f"   - {LOG_JSON_FILE}")
    print(f"   - {LOG_MD_FILE}")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PyTransDatRO Benchmark Runner")
    parser.add_argument("--tag", default="baseline", help="Version or optimization tag (e.g., baseline, opt_minor)")
    parser.add_argument("--desc", default="Baseline performance measurement", help="Description of the run")
    parser.add_argument("--runs", type=int, default=3, help="Number of repetitions per test case (default: 3)")
    args = parser.parse_args()

    run_benchmark(tag=args.tag, description=args.desc, runs=args.runs)
