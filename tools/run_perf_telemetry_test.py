"""Run transformations on perf dataset scenarios (S1 to S4) in both directions with telemetry enabled,
recording both spatial usage and daily call/point volume metrics to ./logs/usage.db.
"""

import sys
import csv
import time
from pathlib import Path

# Ensure repo root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from pytransdatro import TransRO
from pytransdatro.telemetry import SqliteUsageLogger

DATA_DIR = REPO_ROOT / "tools" / "perf" / "data"
LOGS_DIR = REPO_ROOT / "logs"

# Mapping scenarios to source IDs (1=web_map, 2=rest_api, 3=desktop_cli)
SCENARIO_FILES = [
    ("S1_micro", "s1_micro_clustered.csv", 1, "web_map"),
    ("S2_small", "s2_small_clustered.csv", 2, "rest_api"),
    ("S3_medium", "s3_medium_corridor.csv", 2, "rest_api"),
    ("S4_large", "s4_large_dispersed.csv", 3, "desktop_cli"),
]


def load_csv_data(filepath):
    """Loads points from CSV into (northings, eastings, z_values)."""
    pids, northings, eastings, z_values = [], [], [], []
    with open(filepath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            pids.append(row[0])
            northings.append(float(row[1]))
            eastings.append(float(row[2]))
            z_values.append(float(row[3]))
    return northings, eastings, z_values


def main():
    print(f"Initializing SqliteUsageLogger targeting: {LOGS_DIR} ...")
    logger = SqliteUsageLogger(log_dir=str(LOGS_DIR))
    tr = TransRO(telemetry_logger=logger)

    total_pts_transformed = 0
    start_time = time.perf_counter()

    for name, filename, source_id, source_name in SCENARIO_FILES:
        filepath = DATA_DIR / filename
        if not filepath.exists():
            print(f"Warning: File {filepath} not found, skipping.")
            continue

        northings, eastings, z_values = load_csv_data(filepath)
        n_pts = len(northings)
        print(f"\nProcessing '{name}' ({filename}) with {n_pts:,} points [Source {source_id}: {source_name}]...")

        # 1. Forward transformation: Stereo70 -> ETRS89
        t0 = time.perf_counter()
        lats, lons, heights = tr.st70_to_etrs89(
            northings,
            eastings,
            z_values,
            source=source_id
        )
        t_fwd = time.perf_counter() - t0
        print(f"  [+] Stereo70 -> ETRS89 (Forward): {n_pts:,} points in {t_fwd*1000:.2f} ms")

        # 2. Reverse transformation: ETRS89 -> Stereo70
        t0 = time.perf_counter()
        n_rev, e_rev, h_rev = tr.etrs89_to_st70(
            lats,
            lons,
            heights,
            source=source_id
        )
        t_rev = time.perf_counter() - t0
        print(f"  [+] ETRS89 -> Stereo70 (Reverse): {n_pts:,} points in {t_rev*1000:.2f} ms")

        total_pts_transformed += n_pts * 2

    elapsed = time.perf_counter() - start_time
    print(f"\nTotal points transformed (both ways): {total_pts_transformed:,} in {elapsed:.3f} s")

    print("\nFlushing telemetry queue and closing database...")
    logger.close()
    print("Database flushed and closed successfully!")


if __name__ == "__main__":
    main()
