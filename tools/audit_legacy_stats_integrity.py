#!/usr/bin/env python3
"""Audit data integrity and identify time gaps in the legacy statistics database.

Scans legacy_stats.db for:
1. Missing date intervals (days where no record exists between min and max date).
2. Zero count intervals (contiguous periods where records exist but all 4 operation counts are 0).
3. Combined zero or missing intervals (treating missing dates and all-zero days as unrecorded activity).

Outputs a comprehensive text report into the reports/ folder.

Usage:
    python tools/audit_legacy_stats_integrity.py [--db PATH] [--output PATH]
"""

import argparse
import sqlite3
import sys
from collections import defaultdict
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def get_default_paths() -> Tuple[Path, Path]:
    repo_root = Path(__file__).resolve().parent.parent
    default_db = repo_root / "logs" / "legacy_stats.db"
    default_output = repo_root / "reports" / "legacy_stats_integrity_report.txt"
    return default_db, default_output


class DailyRecord:
    __slots__ = ("day", "total_count", "app_st70", "app_st30", "srv_st70", "srv_st30", "row_count")

    def __init__(
        self,
        day: date,
        app_st70: int,
        app_st30: int,
        srv_st70: int,
        srv_st30: int,
        row_count: int = 1,
    ):
        self.day = day
        self.app_st70 = app_st70
        self.app_st30 = app_st30
        self.srv_st70 = srv_st70
        self.srv_st30 = srv_st30
        self.total_count = app_st70 + app_st30 + srv_st70 + srv_st30
        self.row_count = row_count


def load_legacy_records(db_path: Path) -> Tuple[Dict[date, DailyRecord], List[Tuple[date, int]]]:
    """Loads daily statistics aggregated by calendar day and reports duplicates."""
    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {db_path}")

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    query = """
        SELECT day, cooOpAppSt70, cooOpAppSt30, cooOpServiceSt70, cooOpServiceSt30
        FROM statistics
        ORDER BY day ASC, id ASC
    """
    rows = cursor.execute(query).fetchall()
    conn.close()

    if not rows:
        raise ValueError(f"No records found in table 'statistics' in {db_path}")

    daily_map: Dict[date, DailyRecord] = {}
    day_counts: Dict[date, int] = defaultdict(int)

    for day_str, app70, app30, srv70, srv30 in rows:
        d = datetime.strptime(day_str, "%Y-%m-%d").date()
        day_counts[d] += 1
        if d not in daily_map:
            daily_map[d] = DailyRecord(d, app70, app30, srv70, srv30, row_count=1)
        else:
            rec = daily_map[d]
            rec.app_st70 += app70
            rec.app_st30 += app30
            rec.srv_st70 += srv70
            rec.srv_st30 += srv30
            rec.total_count += app70 + app30 + srv70 + srv30
            rec.row_count += 1

    duplicates = [(d, count) for d, count in day_counts.items() if count > 1]
    duplicates.sort(key=lambda x: x[0])
    return daily_map, duplicates


class PeriodGap:
    __slots__ = ("start_date", "end_date", "days", "gap_type")

    def __init__(self, start_date: date, end_date: date, gap_type: str = "missing"):
        self.start_date = start_date
        self.end_date = end_date
        self.days = (end_date - start_date).days + 1
        self.gap_type = gap_type


def detect_periods(
    start_date: date,
    end_date: date,
    daily_map: Dict[date, DailyRecord],
    predicate_fn,
    gap_type: str,
) -> List[PeriodGap]:
    """Scans the date range continuously and clusters contiguous dates matching predicate."""
    periods: List[PeriodGap] = []
    current_start: Optional[date] = None
    current_prev: Optional[date] = None

    curr = start_date
    one_day = timedelta(days=1)

    while curr <= end_date:
        rec = daily_map.get(curr)
        match = predicate_fn(rec)
        if match:
            if current_start is None:
                current_start = curr
            current_prev = curr
        else:
            if current_start is not None and current_prev is not None:
                periods.append(PeriodGap(current_start, current_prev, gap_type))
                current_start = None
                current_prev = None
        curr += one_day

    if current_start is not None and current_prev is not None:
        periods.append(PeriodGap(current_start, current_prev, gap_type))

    return periods


def format_period_table(periods: List[PeriodGap], max_display: int = 50) -> List[str]:
    lines = []
    lines.append(f"{'#':<4} | {'Period Start':<12} | {'Period End':<12} | {'Gap Days':<10}")
    lines.append("-" * 46)
    for idx, p in enumerate(periods, start=1):
        lines.append(
            f"{idx:<4} | {p.start_date.strftime('%Y-%m-%d'):<12} | {p.end_date.strftime('%Y-%m-%d'):<12} | {p.days:<10}"
        )
        if idx >= max_display and len(periods) > max_display:
            lines.append(f"... and {len(periods) - max_display} additional smaller intervals.")
            break
    return lines


def generate_report(
    db_path: Path,
    output_path: Path,
    quiet: bool = False,
) -> str:
    daily_map, duplicates = load_legacy_records(db_path)
    all_dates = sorted(daily_map.keys())
    min_date = all_dates[0]
    max_date = all_dates[-1]
    calendar_days = (max_date - min_date).days + 1
    recorded_days = len(daily_map)
    missing_days_count = calendar_days - recorded_days

    # 1. Missing dates (where no record exists)
    missing_periods = detect_periods(
        min_date, max_date, daily_map, predicate_fn=lambda r: r is None, gap_type="missing"
    )

    # 2. Zero count dates (record exists, but all 4 fields sum to 0)
    zero_periods = detect_periods(
        min_date,
        max_date,
        daily_map,
        predicate_fn=lambda r: r is not None and r.total_count == 0,
        gap_type="zero",
    )
    total_zero_days = sum(p.days for p in zero_periods)

    # 3. Combined periods (missing record OR record exists with all 4 fields == 0)
    combined_periods = detect_periods(
        min_date,
        max_date,
        daily_map,
        predicate_fn=lambda r: r is None or r.total_count == 0,
        gap_type="missing_or_zero",
    )
    total_combined_days = sum(p.days for p in combined_periods)

    # Sorting for breakdown
    # Chronological: default
    # Length: sorted descending
    missing_sorted_by_len = sorted(missing_periods, key=lambda p: p.days, reverse=True)
    zero_sorted_by_len = sorted(zero_periods, key=lambda p: p.days, reverse=True)
    combined_sorted_by_len = sorted(combined_periods, key=lambda p: p.days, reverse=True)

    # Generate document content
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sep = "=" * 80
    subsep = "-" * 80

    out = []
    out.append(sep)
    out.append(" TRANSDATONLINE LEGACY DATABASE INTEGRITY AUDIT REPORT")
    out.append(sep)
    out.append(f" Generated:          {now_str}")
    out.append(f" Database Inspected: {db_path.resolve()}")
    out.append(f" Output Report:      {output_path.resolve()}")
    out.append(subsep)
    out.append(" EXECUTIVE SUMMARY")
    out.append(subsep)
    out.append(f" - Range of Coverage:            {min_date} to {max_date}")
    out.append(
        f" - Total Calendar Span:          {calendar_days:,} days ({calendar_days / 365.25:.2f} years)"
    )
    out.append(
        f" - Distinct Recorded Days:       {recorded_days:,} days ({recorded_days / calendar_days * 100:.2f}%)"
    )
    out.append(
        f" - Missing Calendar Days:        {missing_days_count:,} days ({missing_days_count / calendar_days * 100:.2f}%)"
    )
    out.append(f" - Recorded Days with 0 Count:   {total_zero_days:,} days")
    out.append(
        f" - Inactive / Void Days Total:   {total_combined_days:,} days ({total_combined_days / calendar_days * 100:.2f}%)"
    )
    if duplicates:
        out.append(
            f" - Days with Duplicate Rows:     {len(duplicates)} days ({sum(c for _, c in duplicates)} total rows)"
        )
    else:
        out.append(" - Days with Duplicate Rows:     None (1 row per calendar day)")
    out.append("")

    if duplicates:
        out.append(subsep)
        out.append(" DUPLICATE CALENDAR DAYS IN SOURCE DATABASE")
        out.append(subsep)
        out.append(
            " Note: These dates contain multiple rows in `statistics`. The values were summed."
        )
        for d, count in duplicates:
            rec = daily_map[d]
            out.append(
                f"  - {d}: {count} rows | AppSt70={rec.app_st70:,}, AppSt30={rec.app_st30:,}, SrvSt70={rec.srv_st70:,}, SrvSt30={rec.srv_st30:,}"
            )
        out.append("")

    # REPORT 1: MISSING DATES
    out.append(sep)
    out.append(" REPORT 1: MISSING DATES AUDIT (Gaps in Calendar Timeline)")
    out.append(sep)
    out.append(f" Total Missing Intervals: {len(missing_periods):,}")
    out.append(f" Total Days Lost:         {missing_days_count:,}")
    if missing_periods:
        out.append(
            f" Largest Gap:             {missing_sorted_by_len[0].days} days ({missing_sorted_by_len[0].start_date} to {missing_sorted_by_len[0].end_date})"
        )
    out.append("")
    out.append(" Top 20 Longest Missing Periods:")
    out.append(subsep)
    out.extend(format_period_table(missing_sorted_by_len, max_display=20))
    out.append("")
    out.append(" Complete Chronological List of Missing Periods:")
    out.append(subsep)
    out.extend(format_period_table(missing_periods, max_display=len(missing_periods)))
    out.append("")

    # REPORT 2: ZERO DATES
    out.append(sep)
    out.append(" REPORT 2: ZERO COUNT DATES AUDIT (Days Recorded with 0 Operations)")
    out.append(sep)
    out.append(f" Total Zero Intervals:    {len(zero_periods):,}")
    out.append(f" Total Zero Days:         {total_zero_days:,}")
    if zero_periods:
        out.append(
            f" Largest Zero Period:     {zero_sorted_by_len[0].days} days ({zero_sorted_by_len[0].start_date} to {zero_sorted_by_len[0].end_date})"
        )
    out.append("")
    out.append(" Chronological List of Zero Count Periods:")
    out.append(subsep)
    out.extend(format_period_table(zero_periods, max_display=len(zero_periods)))
    out.append("")

    # REPORT 3: COMBINED (MISSING + ZERO)
    out.append(sep)
    out.append(" REPORT 3: COMBINED AUDIT (Missing Dates + Zero Count Dates)")
    out.append(sep)
    out.append(" Treats both missing calendar records and all-zero records as void periods.")
    out.append(f" Total Void Intervals:    {len(combined_periods):,}")
    out.append(f" Total Void Days:         {total_combined_days:,}")
    if combined_periods:
        out.append(
            f" Largest Void Interval:   {combined_sorted_by_len[0].days} days ({combined_sorted_by_len[0].start_date} to {combined_sorted_by_len[0].end_date})"
        )
    out.append("")
    out.append(" Top 20 Longest Void Periods:")
    out.append(subsep)
    out.extend(format_period_table(combined_sorted_by_len, max_display=20))
    out.append("")
    out.append(" Complete Chronological List of Void Periods:")
    out.append(subsep)
    out.extend(format_period_table(combined_periods, max_display=len(combined_periods)))
    out.append("")
    out.append(sep)
    out.append(" END OF AUDIT REPORT")
    out.append(sep)

    report_text = "\n".join(out) + "\n"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(report_text)

    if not quiet:
        print(f"Audit report successfully written to: {output_path.resolve()}")
        print(
            f"Summary: {calendar_days:,} calendar days, {missing_days_count:,} missing days ({len(missing_periods)} gaps), {total_zero_days} zero days ({len(zero_periods)} periods), {total_combined_days} combined void days."
        )

    return report_text


def main():
    default_db, default_output = get_default_paths()
    parser = argparse.ArgumentParser(
        description="Audit legacy TransdatOnline stats database for missing dates and zero intervals."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=default_db,
        help=f"Path to SQLite database (default: {default_db})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_output,
        help=f"Path for output report text file (default: {default_output})",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress console output",
    )

    args = parser.parse_args()
    try:
        generate_report(args.db, args.output, quiet=args.quiet)
    except Exception as exc:
        print(f"Error executing audit: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
