import calendar
import datetime
import math
import sqlite3

import pytest

from pytransdatro import TransRO
from pytransdatro.exceptions import OutOfGridErr
from pytransdatro.telemetry import SqliteUsageLogger


def test_telemetry_batching_and_binning(tmp_path):
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    n, e = 500000.0, 300000.0
    # Simulate 10 individual calls within the same 1km cell with source=2 (web_map)
    for _ in range(10):
        tr.st70_to_etrs89(n, e, source=2)
        n += 10.0  # 10 meters should fall in the same bin

    logger.close()  # drain queue and close DB

    db_path = log_dir / "usage.db"
    assert db_path.exists()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Verify spatial_usage table (spatial ~1 km bins with date range)
    cursor.execute("""
        SELECT period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count
        FROM spatial_usage
    """)
    rows = cursor.fetchall()
    assert len(rows) == 1
    row = rows[0]
    now = datetime.date.today()
    last_day = calendar.monthrange(now.year, now.month)[1]
    expected_start = f"{now.year:04d}-{now.month:02d}-01"
    expected_end = f"{now.year:04d}-{now.month:02d}-{last_day:02d}"

    assert row[0] == expected_start
    assert row[1] == expected_end
    assert row[6] == 1  # operation: 1 (st70_to_etrs89)
    assert row[7] == 0  # is_3d
    assert row[8] == 2  # source
    assert row[9] == 10  # count

    # Verify computation_usage table
    cursor.execute("""
        SELECT period_start, period_end, operation, is_3d, source, count_calls, count_points
        FROM computation_usage
    """)
    comp_rows = cursor.fetchall()
    assert len(comp_rows) == 1
    c_row = comp_rows[0]
    today_str = now.strftime("%Y-%m-%d")
    assert c_row[0] == today_str
    assert c_row[1] == today_str
    assert c_row[2] == 1  # operation 1 (st70_to_etrs89)
    assert c_row[3] == 0  # is_3d 0
    assert c_row[4] == 2  # source 2 (web_map)
    assert c_row[5] == 10  # 10 calls
    assert c_row[6] == 10  # 10 points

    # Verify sources lookup table
    cursor.execute("SELECT id, name FROM sources WHERE id = 2")
    src = cursor.fetchone()
    assert src == (2, "web_map")

    cursor.execute("SELECT id, name FROM sources WHERE id = 1")
    src1 = cursor.fetchone()
    assert src1 == (1, "rest_api")

    # Verify operations lookup table
    cursor.execute("SELECT id, name FROM operations ORDER BY id")
    ops = cursor.fetchall()
    assert ops == [(1, "st70_to_etrs89"), (2, "etrs89_to_st70")]

    # Verify telemetry_config table
    cursor.execute("SELECT metric_type, window_days FROM telemetry_config ORDER BY metric_type")
    configs = cursor.fetchall()
    assert ("computation", 1) in configs
    assert ("spatial", 30) in configs

    conn.close()


def test_telemetry_failure_isolation():
    class BrokenLogger:
        def record(self, operation, lats, lons, is_3d, source=0):
            raise RuntimeError("Disk failed!")

    tr = TransRO(telemetry_logger=BrokenLogger())
    n, e = 500000.0, 300000.0
    lat, lon = tr.st70_to_etrs89(n, e, source=2)  # Should not crash
    assert lat is not None
    assert lon is not None

    # Also verify etrs89_to_st70 isolates logger failure
    lat_in, lon_in = math.radians(45.0), math.radians(25.0)
    n_res, e_res = tr.etrs89_to_st70(lat_in, lon_in, source=2)  # Should not crash
    assert n_res is not None
    assert e_res is not None


def test_telemetry_etrs89_to_st70_out_of_grid_not_recorded(tmp_path):
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    # 1. Scalar 2D out of grid (0.0, 0.0)
    with pytest.raises(OutOfGridErr):
        tr.etrs89_to_st70(0.0, 0.0, source=1)

    # 2. Scalar 3D far out of grid
    with pytest.raises(OutOfGridErr):
        tr.etrs89_to_st70(math.radians(80.0), math.radians(80.0), h=100.0, source=2)

    # 3. Batch out of grid
    with pytest.raises(OutOfGridErr):
        tr.etrs89_to_st70([0.0, 0.1], [0.0, 0.1], source=1)

    # 4. Mixed batch (one valid in-grid point, one out-of-grid point)
    with pytest.raises(OutOfGridErr):
        tr.etrs89_to_st70([math.radians(45.0), 0.0], [math.radians(25.0), 0.0], source=2)

    logger.close()

    db_path = log_dir / "usage.db"
    assert db_path.exists()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    cursor.execute("SELECT count(*) FROM spatial_usage")
    assert cursor.fetchone()[0] == 0

    cursor.execute("SELECT count(*) FROM computation_usage")
    assert cursor.fetchone()[0] == 0

    conn.close()


def test_telemetry_st70_to_etrs89_out_of_grid_not_recorded(tmp_path):
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    # 1. Scalar 2D out of grid Stereo 70
    with pytest.raises(OutOfGridErr):
        tr.st70_to_etrs89(5.0, 1.0, source=1)

    # 2. Scalar 3D out of grid Stereo 70
    with pytest.raises(OutOfGridErr):
        tr.st70_to_etrs89(0.0, 0.0, z=50.0, source=2)

    # 3. Boundary coordinates outside Stereo 70 grid coverage
    with pytest.raises(OutOfGridErr):
        tr.st70_to_etrs89(100000.0, 100000.0, source=1)

    with pytest.raises(OutOfGridErr):
        tr.st70_to_etrs89(900000.0, 900000.0, source=1)

    # 4. Batch out of grid Stereo 70
    with pytest.raises(OutOfGridErr):
        tr.st70_to_etrs89([5.0, 10.0], [1.0, 2.0], source=3)

    # 5. Mixed batch Stereo 70 (one valid in-grid point, one out-of-grid point)
    with pytest.raises(OutOfGridErr):
        tr.st70_to_etrs89([500000.0, 5.0], [300000.0, 1.0], source=2)

    logger.close()

    db_path = log_dir / "usage.db"
    assert db_path.exists()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    cursor.execute("SELECT count(*) FROM spatial_usage")
    assert cursor.fetchone()[0] == 0

    cursor.execute("SELECT count(*) FROM computation_usage")
    assert cursor.fetchone()[0] == 0

    conn.close()


def test_telemetry_in_grid_both_directions_recorded(tmp_path):
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    # In-grid ETRS89 to Stereo 70 2D (~45 deg N, ~25 deg E) -> operation 2, is_3d 0
    lat = math.radians(45.0)
    lon = math.radians(25.0)
    n_res, e_res = tr.etrs89_to_st70(lat, lon, source=2)
    assert n_res > 0 and e_res > 0

    # In-grid Stereo 70 to ETRS89 2D -> operation 1, is_3d 0
    lat_res, lon_res = tr.st70_to_etrs89(500000.0, 300000.0, source=1)
    assert lat_res > 0 and lon_res > 0

    # In-grid 3D ETRS89 to Stereo 70 -> operation 2, is_3d 1
    n_3d, e_3d, h_3d = tr.etrs89_to_st70(lat, lon, h=150.0, source=2)
    assert n_3d > 0 and e_3d > 0 and h_3d is not None

    # In-grid 3D Stereo 70 to ETRS89 -> operation 1, is_3d 1
    lat_3d, lon_3d, z_3d = tr.st70_to_etrs89(500000.0, 300000.0, z=150.0, source=1)
    assert lat_3d > 0 and lon_3d > 0 and z_3d is not None

    # In-grid batch call (2 points) -> operation 1, is_3d 0
    b_lat, b_lon = tr.st70_to_etrs89([500000.0, 500100.0], [300000.0, 300100.0], source=3)
    assert len(b_lat) == 2 and len(b_lon) == 2

    logger.close()

    db_path = log_dir / "usage.db"
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    cursor.execute("SELECT count(*) FROM spatial_usage")
    assert cursor.fetchone()[0] > 0

    # Check that both operations were recorded in spatial_usage table (1 and 2)
    cursor.execute("SELECT DISTINCT operation FROM spatial_usage ORDER BY operation")
    operations = [row[0] for row in cursor.fetchall()]
    assert operations == [1, 2]  # 1 = st70_to_etrs89, 2 = etrs89_to_st70

    # Check 3D was recorded
    cursor.execute("SELECT DISTINCT is_3d FROM spatial_usage ORDER BY is_3d")
    is_3d_values = [row[0] for row in cursor.fetchall()]
    assert is_3d_values == [0, 1]

    # Check computation_usage breakdown:
    # (operation, is_3d, source, count_calls, count_points)
    # op 1, 2D, src 1: 1 call, 1 point
    # op 1, 3D, src 1: 1 call, 1 point
    # op 2, 2D, src 2: 1 call, 1 point
    # op 2, 3D, src 2: 1 call, 1 point
    # op 1, 2D, src 3: 1 call, 2 points
    cursor.execute("""
        SELECT operation, is_3d, source, count_calls, count_points
        FROM computation_usage
        ORDER BY operation, is_3d, source
    """)
    comp_stats = cursor.fetchall()
    assert comp_stats == [
        (1, 0, 1, 1, 1),
        (1, 0, 3, 1, 2),
        (1, 1, 1, 1, 1),
        (2, 0, 2, 1, 1),
        (2, 1, 2, 1, 1),
    ]

    conn.close()


def test_telemetry_legacy_schema_migration(tmp_path):
    """Test automatic migration from legacy daily_stats and usage tables to new schema."""
    log_dir = tmp_path / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    db_path = log_dir / "usage.db"

    # Pre-create legacy database
    conn = sqlite3.connect(str(db_path))
    conn.execute("""
        CREATE TABLE daily_stats (
            date TEXT NOT NULL,
            source INTEGER NOT NULL,
            count_calls INTEGER NOT NULL DEFAULT 1,
            count_points INTEGER NOT NULL DEFAULT 1,
            PRIMARY KEY (date, source)
        );
    """)
    conn.execute("""
        CREATE TABLE usage (
            year INTEGER NOT NULL,
            month INTEGER NOT NULL,
            lat_deg INTEGER NOT NULL,
            lat_grad INTEGER NOT NULL,
            lon_deg INTEGER NOT NULL,
            lon_grad INTEGER NOT NULL,
            is_st70_to_etrs89 INTEGER NOT NULL,
            is_3d INTEGER NOT NULL,
            source INTEGER NOT NULL DEFAULT 0,
            count INTEGER NOT NULL DEFAULT 1,
            PRIMARY KEY (year, month, lat_deg, lat_grad, lon_deg, lon_grad, is_st70_to_etrs89, is_3d, source)
        );
    """)
    conn.execute("INSERT INTO daily_stats VALUES ('2025-05-15', 2, 5, 25);")
    conn.execute("INSERT INTO usage VALUES (2025, 5, 45, 50, 25, 30, 0, 1, 2, 15);")
    conn.commit()
    conn.close()

    # Now open logger, which should migrate tables automatically
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    logger.close()

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Check legacy tables were dropped
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = {r[0] for r in cursor.fetchall()}
    assert "daily_stats" not in tables
    assert "usage" not in tables
    assert "computation_usage" in tables
    assert "spatial_usage" in tables

    # Check migrated computation_usage
    cursor.execute(
        "SELECT period_start, period_end, source, count_calls, count_points FROM computation_usage"
    )
    comp_row = cursor.fetchone()
    assert comp_row == ("2025-05-15", "2025-05-15", 2, 5, 25)

    # Check migrated spatial_usage (legacy is_st70_to_etrs89=0 converted to operation=2)
    cursor.execute(
        "SELECT period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count FROM spatial_usage"
    )
    spat_row = cursor.fetchone()
    assert spat_row == ("2025-05-01", "2025-05-31", 45, 50, 25, 30, 2, 1, 2, 15)

    conn.close()
