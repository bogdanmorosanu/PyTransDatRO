#!/usr/bin/env python3
"""Extract statistics data from legacy TransdatOnline MySQL dump into an SQLite database.

Aggregates multiple records per calendar day by summing operations, ensuring
a strictly unique (one row per day) normalized timeline.

Usage:
    python tools/extract_legacy_stats.py [--input PATH] [--output PATH] [--quiet]
"""

import argparse
import re
import sqlite3
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Tuple


def get_default_paths() -> Tuple[Path, Path]:
    repo_root = Path(__file__).resolve().parent.parent
    # Input dump stored in logs/localhost.sql (fallback to transdatonline/db logs if missing)
    default_input = repo_root / "logs" / "localhost.sql"
    if not default_input.exists():
        fallback_input = repo_root / "transdatonline" / "db logs" / "localhost.sql"
        if fallback_input.exists():
            default_input = fallback_input
    default_output = repo_root / "logs" / "legacy_stats.db"
    return default_input, default_output


def parse_statistics_rows(sql_path: Path) -> List[Tuple[int, str, int, int, int, int]]:
    """Parses INSERT statements for the `statistics` table from the MySQL dump."""
    if not sql_path.exists():
        raise FileNotFoundError(f"Input SQL file not found: {sql_path}")

    # Read binary and decode with fallback to handle diverse MySQL encodings cleanly
    with open(sql_path, "rb") as f:
        raw_bytes = f.read()

    try:
        text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError:
        text = raw_bytes.decode("latin-1")

    # Regex matching individual row tuples: ('YYYY-MM-DD', int, int, int, int, int)
    # in INSERT INTO `statistics` statements
    pattern = re.compile(
        r"\(\s*'(\d{4}-\d{2}-\d{2})'\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)"
    )

    # Locate the statistics table dump section
    start_token = "INSERT INTO `statistics`"
    start_pos = text.find(start_token)
    if start_pos == -1:
        # Fallback without backticks
        start_token = "INSERT INTO statistics"
        start_pos = text.find(start_token)
        if start_pos == -1:
            raise ValueError("Could not find INSERT INTO `statistics` block in SQL file.")

    # Find the end of statistics section (typically marked by ALTER TABLE or next table structure)
    end_token = "-- Indexes for table `statistics`"
    end_pos = text.find(end_token, start_pos)
    if end_pos == -1:
        stats_section = text[start_pos:]
    else:
        stats_section = text[start_pos:end_pos]

    raw_rows: List[Tuple[int, str, int, int, int, int]] = []
    for match in pattern.finditer(stats_section):
        day, st70_app, st30_app, st70_srv, st30_srv, row_id = match.groups()
        raw_rows.append(
            (
                int(row_id),
                day,
                int(st70_app),
                int(st30_app),
                int(st70_srv),
                int(st30_srv),
            )
        )

    return raw_rows


def aggregate_by_day(
    raw_rows: List[Tuple[int, str, int, int, int, int]],
) -> Tuple[List[Tuple[int, str, int, int, int, int]], Dict[str, int]]:
    """Aggregates multiple records per day by summing the 4 operation counts.

    Returns:
        aggregated_rows: List of (id, day, app70, app30, srv70, srv30) sorted by day.
        duplicate_days: Mapping of day -> count of original raw rows for days with > 1 row.
    """
    daily_sums: Dict[str, List[int]] = OrderedDict()
    daily_row_counts: Dict[str, int] = OrderedDict()

    # Sort raw rows chronologically by day
    sorted_raw = sorted(raw_rows, key=lambda r: (r[1], r[0]))

    for _row_id, day, app70, app30, srv70, srv30 in sorted_raw:
        daily_row_counts[day] = daily_row_counts.get(day, 0) + 1
        if day not in daily_sums:
            daily_sums[day] = [app70, app30, srv70, srv30]
        else:
            daily_sums[day][0] += app70
            daily_sums[day][1] += app30
            daily_sums[day][2] += srv70
            daily_sums[day][3] += srv30

    aggregated_rows: List[Tuple[int, str, int, int, int, int]] = []
    for new_id, (day, counts) in enumerate(daily_sums.items(), start=1):
        aggregated_rows.append((new_id, day, counts[0], counts[1], counts[2], counts[3]))

    duplicate_days = {d: count for d, count in daily_row_counts.items() if count > 1}
    return aggregated_rows, duplicate_days


def create_and_populate_db(
    output_path: Path,
    raw_rows: List[Tuple[int, str, int, int, int, int]],
    quiet: bool = False,
) -> None:
    """Creates SQLite database, populates with daily-aggregated statistics, and asserts zero duplicates."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    aggregated_rows, duplicate_days = aggregate_by_day(raw_rows)

    # Cleanly remove existing file to ensure full overwrite
    if output_path.exists():
        try:
            output_path.unlink()
        except OSError:
            pass

    conn = sqlite3.connect(str(output_path))
    cursor = conn.cursor()

    cursor.executescript("""
        DROP TABLE IF EXISTS statistics;
        CREATE TABLE statistics (
            id INTEGER PRIMARY KEY,
            day TEXT NOT NULL UNIQUE,
            cooOpAppSt70 INTEGER NOT NULL DEFAULT 0,
            cooOpAppSt30 INTEGER NOT NULL DEFAULT 0,
            cooOpServiceSt70 INTEGER NOT NULL DEFAULT 0,
            cooOpServiceSt30 INTEGER NOT NULL DEFAULT 0
        );
        CREATE UNIQUE INDEX idx_statistics_day ON statistics(day);
    """)

    cursor.executemany(
        """
        INSERT INTO statistics (id, day, cooOpAppSt70, cooOpAppSt30, cooOpServiceSt70, cooOpServiceSt30)
        VALUES (?, ?, ?, ?, ?, ?)
    """,
        aggregated_rows,
    )

    conn.commit()

    # Verification query: Assert no duplicate dates exist (GROUP BY day HAVING COUNT(*) >= 2)
    cursor.execute("""
        SELECT day, COUNT(*)
        FROM statistics
        GROUP BY day
        HAVING COUNT(*) >= 2
    """)
    remaining_duplicates = cursor.fetchall()
    if remaining_duplicates:
        conn.close()
        raise RuntimeError(
            f"Integrity check failed: found {len(remaining_duplicates)} duplicate dates remaining in database: {remaining_duplicates}"
        )

    if not quiet:
        # Gather summary metrics
        cursor.execute("""
            SELECT
                COUNT(*),
                MIN(day),
                MAX(day),
                SUM(cooOpAppSt70),
                SUM(cooOpAppSt30),
                SUM(cooOpServiceSt70),
                SUM(cooOpServiceSt30)
            FROM statistics
        """)
        count, min_day, max_day, sum_app70, sum_app30, sum_srv70, sum_srv30 = cursor.fetchone()
        file_size_bytes = output_path.stat().st_size
        total_hits = (sum_app70 or 0) + (sum_app30 or 0) + (sum_srv70 or 0) + (sum_srv30 or 0)

        print("=" * 70)
        print(" LEGACY TRANSDATONLINE STATISTICS EXTRACTION COMPLETE")
        print("=" * 70)
        print(f" Output Database:       {output_path}")
        print(
            f" Database Size:         {file_size_bytes:,} bytes ({file_size_bytes / 1024:.2f} KB)"
        )
        print(f" Raw Rows Parsed:       {len(raw_rows):,} rows from MySQL dump")
        print(f" Aggregated Rows:       {count:,} days (strictly 1 row per calendar day)")
        print(f" Date Range:            {min_day} to {max_day}")
        print("-" * 70)
        print(f" Duplicate Days Fixed:  {len(duplicate_days)} days had multiple raw records:")
        for d, dup_count in duplicate_days.items():
            print(f"   - {d}: {dup_count} records summed")
        print("-" * 70)
        print(" Operation Volumes (Total Summed Points Transformed):")
        print(f"   - App Stereo 70:     {sum_app70:,}")
        print(f"   - App Stereo 30:     {sum_app30:,}")
        print(f"   - Service Stereo 70: {sum_srv70:,}")
        print(f"   - Service Stereo 30: {sum_srv30:,}")
        print(f"   - Total Usage Hits:  {total_hits:,}")
        print("-" * 70)
        print(" Post-Extraction Duplicate Verification:")
        print("   [OK] Query 'GROUP BY day HAVING COUNT(*) >= 2' returned 0 duplicates.")
        print("   [OK] Unique constraint on 'day' column enforced in SQLite schema.")
        print("=" * 70)

    conn.close()


def main():
    default_input, default_output = get_default_paths()

    parser = argparse.ArgumentParser(
        description="Extract statistics data from legacy TransdatOnline MySQL dump into an SQLite DB with day-level aggregation."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        default=default_input,
        help=f"Path to localhost.sql dump (default: {default_input})",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=default_output,
        help=f"Path to output SQLite database (default: {default_output})",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress summary output to stdout",
    )

    args = parser.parse_args()

    try:
        raw_rows = parse_statistics_rows(args.input)
        create_and_populate_db(args.output, raw_rows, quiet=args.quiet)
    except Exception as exc:
        print(f"Error during legacy stats extraction: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
