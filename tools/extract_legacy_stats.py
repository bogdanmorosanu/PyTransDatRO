#!/usr/bin/env python3
"""Extract statistics data from legacy TransdatOnline MySQL dump into an SQLite database.

Usage:
    python tools/extract_legacy_stats.py [--input PATH] [--output PATH] [--quiet]
"""

import argparse
import os
import re
import sqlite3
import sys
from pathlib import Path


def get_default_paths():
    repo_root = Path(__file__).resolve().parent.parent
    default_input = repo_root / "transdatonline" / "db logs" / "localhost.sql"
    default_output = repo_root / "logs" / "legacy_stats.db"
    return default_input, default_output


def parse_statistics_rows(sql_path: Path):
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

    rows = []
    for match in pattern.finditer(stats_section):
        day, st70_app, st30_app, st70_srv, st30_srv, row_id = match.groups()
        rows.append((
            int(row_id),
            day,
            int(st70_app),
            int(st30_app),
            int(st70_srv),
            int(st30_srv),
        ))

    return rows


def create_and_populate_db(output_path: Path, rows: list, quiet: bool = False):
    """Creates SQLite database and populates it with extracted statistics rows."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Cleanly remove existing file to ensure full overwrite
    if output_path.exists():
        try:
            output_path.unlink()
        except OSError as e:
            # If open or locked, try to open and truncate
            pass

    conn = sqlite3.connect(str(output_path))
    cursor = conn.cursor()

    cursor.executescript("""
        DROP TABLE IF EXISTS statistics;
        CREATE TABLE statistics (
            id INTEGER PRIMARY KEY,
            day TEXT NOT NULL,
            cooOpAppSt70 INTEGER NOT NULL DEFAULT 0,
            cooOpAppSt30 INTEGER NOT NULL DEFAULT 0,
            cooOpServiceSt70 INTEGER NOT NULL DEFAULT 0,
            cooOpServiceSt30 INTEGER NOT NULL DEFAULT 0
        );
        CREATE INDEX idx_statistics_day ON statistics(day);
    """)

    cursor.executemany("""
        INSERT INTO statistics (id, day, cooOpAppSt70, cooOpAppSt30, cooOpServiceSt70, cooOpServiceSt30)
        VALUES (?, ?, ?, ?, ?, ?)
    """, rows)

    conn.commit()

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
        print(f" Output Database: {output_path}")
        print(f" Database Size:   {file_size_bytes:,} bytes ({file_size_bytes / 1024:.2f} KB)")
        print(f" Total Rows:      {count:,}")
        print(f" Date Range:      {min_day} to {max_day}")
        print("-" * 70)
        print(" Operation Volumes:")
        print(f"   - App Stereo 70:     {sum_app70:,}")
        print(f"   - App Stereo 30:     {sum_app30:,}")
        print(f"   - Service Stereo 70: {sum_srv70:,}")
        print(f"   - Service Stereo 30: {sum_srv30:,}")
        print(f"   - Total Usage Hits:  {total_hits:,}")
        print("=" * 70)

    conn.close()


def main():
    default_input, default_output = get_default_paths()

    parser = argparse.ArgumentParser(
        description="Extract statistics data from legacy TransdatOnline MySQL dump into an SQLite DB."
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        default=default_input,
        help=f"Path to localhost.sql dump (default: {default_input})"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=default_output,
        help=f"Path to output SQLite database (default: {default_output})"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress summary output to stdout"
    )

    args = parser.parse_args()

    rows = parse_statistics_rows(args.input)
    create_and_populate_db(args.output, rows, quiet=args.quiet)


if __name__ == "__main__":
    main()
