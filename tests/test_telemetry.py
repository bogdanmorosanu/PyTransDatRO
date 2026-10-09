import sqlite3
import pytest
from pathlib import Path
from pytransdatro import TransRO
from pytransdatro.telemetry import SqliteUsageLogger


def test_telemetry_batching_and_binning(tmp_path):
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    n, e = 500000.0, 300000.0
    # Simulate 10 individual calls within the same 1km cell with source=1 (web_map)
    for _ in range(10):
        tr.st70_to_etrs89(n, e, source=1)
        n += 10.0  # 10 meters should fall in the same bin

    logger.close()  # drain queue and close DB

    db_path = log_dir / "usage.db"
    assert db_path.exists()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Verify usage table (spatial ~1 km bins)
    cursor.execute("""
        SELECT year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source, count
        FROM usage
    """)
    rows = cursor.fetchall()
    assert len(rows) == 1
    row = rows[0]
    assert row[6] == 1  # is_st70_to_etrs89
    assert row[7] == 0  # is_3d
    assert row[8] == 1  # source
    assert row[9] == 10  # count

    # Verify daily_stats table
    cursor.execute("SELECT date, source, count_calls, count_points FROM daily_stats")
    daily_rows = cursor.fetchall()
    assert len(daily_rows) == 1
    d_row = daily_rows[0]
    assert d_row[1] == 1   # source 1 (web_map)
    assert d_row[2] == 10  # 10 calls
    assert d_row[3] == 10  # 10 points

    # Verify sources lookup table
    cursor.execute("SELECT id, name FROM sources WHERE id = 1")
    src = cursor.fetchone()
    assert src == (1, "web_map")

    conn.close()


def test_telemetry_failure_isolation():
    class BrokenLogger:
        def record(self, is_st70_to_etrs89, lats, lons, is_3d, source=0):
            raise RuntimeError("Disk failed!")

    tr = TransRO(telemetry_logger=BrokenLogger())
    n, e = 500000.0, 300000.0
    lat, lon = tr.st70_to_etrs89(n, e, source=2)  # Should not crash
    assert lat is not None
    assert lon is not None
