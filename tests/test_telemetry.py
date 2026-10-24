import calendar
import datetime
import logging
import math
import sqlite3
import threading
import time
from unittest.mock import MagicMock

import pytest

import pytransdatro.exceptions as exceptions
from pytransdatro import TransRO
from pytransdatro.exceptions import OutOfGridErr
from pytransdatro.telemetry import BufferedTelemetryCollector, SqliteUsageLogger


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
        SELECT period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count, count_calls
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
    assert row[9] == 10  # count (points)
    assert row[10] == 10  # count_calls

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
    assert ("spatial", -1) in configs

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


def test_telemetry_batch_survivor_persistence(tmp_path):
    """Category A: Mixed batch persists valid survivors to SQLite with correct counts."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    valid_points = [
        [500000.0, 300000.0],
        [500100.0, 300100.0],
        [500200.0, 300200.0],
        [500300.0, 300300.0],
        [500400.0, 300400.0],
        [500500.0, 300500.0],
    ]
    invalid_points = [
        [5.0, 1.0],  # Out of grid
        [float("nan"), 300000.0],  # NaN
        [500000.0],  # Wrong dimension (<2)
        ["invalid", "coord"],  # Non-numeric
    ]
    mixed = valid_points[:3] + invalid_points[:2] + valid_points[3:] + invalid_points[2:]

    results = tr.transform_batch(mixed, from_st70=True, source=1)
    assert len(results) == 10

    logger.close()
    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()

    cur.execute("SELECT count_calls, count_points FROM computation_usage")
    comp = cur.fetchone()
    assert comp == (1, 6)

    cur.execute("SELECT SUM(count) FROM spatial_usage")
    assert cur.fetchone()[0] == 6
    conn.close()


def test_telemetry_invalid_data_rejection_no_db_writes(tmp_path):
    """Category A: NaN, Infinity, strings, and bad dimension raise errors and write 0 rows."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    for bad_val in (float("nan"), float("inf"), -float("inf"), "bad_coord"):
        with pytest.raises(exceptions.InvalidCoordinateErr):
            tr.st70_to_etrs89(bad_val, 300000.0, source=1)

    with pytest.raises(exceptions.WrongDimensionErr):
        tr.st70_to_etrs89([500000.0], [300000.0, 300100.0], source=1)

    logger.close()
    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()
    cur.execute("SELECT count(*) FROM computation_usage")
    assert cur.fetchone()[0] == 0
    cur.execute("SELECT count(*) FROM spatial_usage")
    assert cur.fetchone()[0] == 0
    conn.close()


def test_telemetry_nodata_grid_holes_not_recorded(tmp_path, monkeypatch):
    """Category A: Points hitting SPG nodata cells raise NoDataGridErr and write 0 rows."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    def mock_trans_nodata(*args, **kwargs):
        raise exceptions.NoDataGridErr(500000.0, 300000.0, tr._t_gr2d)

    monkeypatch.setattr(tr._t_gr2d, "trans", mock_trans_nodata)

    with pytest.raises(exceptions.NoDataGridErr):
        tr.st70_to_etrs89(500000.0, 300000.0, source=1)

    logger.close()
    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()
    cur.execute("SELECT count(*) FROM computation_usage")
    assert cur.fetchone()[0] == 0
    cur.execute("SELECT count(*) FROM spatial_usage")
    assert cur.fetchone()[0] == 0
    conn.close()


def test_telemetry_st70_to_etrs89_output_radians_logged(tmp_path):
    """Category B: Stereo 70 to ETRS89 logs output coordinates in radians."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    n, e = 500000.0, 300000.0
    lat_rad, lon_rad = tr.st70_to_etrs89(n, e, source=1)

    logger.close()
    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()
    cur.execute("SELECT lat_deg, lat_grad, lon_deg, lon_grad FROM spatial_usage")
    row = cur.fetchone()

    total_lat = int(round(math.degrees(lat_rad) * 100))
    total_lon = int(round(math.degrees(lon_rad) * 100))
    assert row == (total_lat // 100, total_lat % 100, total_lon // 100, total_lon % 100)
    conn.close()


def test_telemetry_etrs89_to_st70_deferred_logging(tmp_path, monkeypatch):
    """Category B: ETRS89 to Stereo 70 defers logging until after successful calculation."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    lat, lon = math.radians(45.0), math.radians(25.0)

    def mock_trans_err(*args, **kwargs):
        raise exceptions.NoDataGridErr(45.0, 25.0, tr._t_gr2d)

    monkeypatch.setattr(tr._t_gr2d, "trans", mock_trans_err)

    with pytest.raises(exceptions.NoDataGridErr):
        tr.etrs89_to_st70(lat, lon, source=1)

    logger.close()
    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()
    cur.execute("SELECT count(*) FROM spatial_usage")
    assert cur.fetchone()[0] == 0
    conn.close()


def test_telemetry_spatial_cell_bidirectional_symmetry(tmp_path):
    """Category B: Forward and inverse transformations bin to the exact same 1km cell."""
    log_dir_fwd = tmp_path / "logs_fwd"
    logger_fwd = SqliteUsageLogger(log_dir=str(log_dir_fwd))
    tr_fwd = TransRO(telemetry_logger=logger_fwd)

    n, e = 500000.0, 300000.0
    lat_rad, lon_rad = tr_fwd.st70_to_etrs89(n, e, source=1)
    logger_fwd.close()

    log_dir_inv = tmp_path / "logs_inv"
    logger_inv = SqliteUsageLogger(log_dir=str(log_dir_inv))
    tr_inv = TransRO(telemetry_logger=logger_inv)
    tr_inv.etrs89_to_st70(lat_rad, lon_rad, source=1)
    logger_inv.close()

    conn_fwd = sqlite3.connect(str(log_dir_fwd / "usage.db"))
    cell_fwd = conn_fwd.execute(
        "SELECT lat_deg, lat_grad, lon_deg, lon_grad FROM spatial_usage"
    ).fetchone()
    conn_fwd.close()

    conn_inv = sqlite3.connect(str(log_dir_inv / "usage.db"))
    cell_inv = conn_inv.execute(
        "SELECT lat_deg, lat_grad, lon_deg, lon_grad FROM spatial_usage"
    ).fetchone()
    conn_inv.close()

    assert cell_fwd is not None
    assert cell_inv is not None
    assert cell_fwd == cell_inv


def test_telemetry_angle_units_invariance(tmp_path):
    """Category B: Degrees and radians inputs normalize to the exact same spatial cell."""
    from api.schemas.point import AngleUnitEnum
    from api.services.transformer import execute_point_transformation

    lat_deg, lon_deg = 45.9997186, 24.9984476
    lat_rad, lon_rad = math.radians(lat_deg), math.radians(lon_deg)

    # 1. Via degrees
    log_dir_deg = tmp_path / "logs_deg"
    logger_deg = SqliteUsageLogger(log_dir=str(log_dir_deg))
    tr_deg = TransRO(telemetry_logger=logger_deg)
    execute_point_transformation(
        tr_deg, "ETRS89ToStereo70", [lat_deg, lon_deg], unit=AngleUnitEnum.DEGREES, source=1
    )
    logger_deg.close()

    # 2. Via radians directly
    log_dir_rad = tmp_path / "logs_rad"
    logger_rad = SqliteUsageLogger(log_dir=str(log_dir_rad))
    tr_rad = TransRO(telemetry_logger=logger_rad)
    tr_rad.etrs89_to_st70(lat_rad, lon_rad, source=1)
    logger_rad.close()

    conn_deg = sqlite3.connect(str(log_dir_deg / "usage.db"))
    cell_deg = conn_deg.execute(
        "SELECT lat_deg, lat_grad, lon_deg, lon_grad FROM spatial_usage"
    ).fetchone()
    conn_deg.close()

    conn_rad = sqlite3.connect(str(log_dir_rad / "usage.db"))
    cell_rad = conn_rad.execute(
        "SELECT lat_deg, lat_grad, lon_deg, lon_grad FROM spatial_usage"
    ).fetchone()
    conn_rad.close()

    assert cell_deg is not None
    assert cell_rad is not None
    assert cell_deg == cell_rad


def test_telemetry_non_blocking_latency_benchmark(tmp_path):
    """Category C: Telemetry record() is non-blocking with sub-0.05ms average latency."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))

    lat, lon = [math.radians(45.0)], [math.radians(25.0)]
    # Warm-up
    logger.record(1, lat, lon, is_3d=0)

    # 200 calls benchmark
    n_iterations = 200
    t0 = time.perf_counter()
    for _ in range(n_iterations):
        logger.record(1, lat, lon, is_3d=0)
    elapsed = time.perf_counter() - t0

    avg_ms = (elapsed / n_iterations) * 1000.0
    logger.close()
    assert avg_ms < 0.05


def test_telemetry_caller_isolation_during_stalled_disk_io(tmp_path, monkeypatch):
    """Category C: Worker disk I/O stalls do not block coordinate transformation callers."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    original_process = logger._process_batch

    def slow_process(*args, **kwargs):
        time.sleep(0.5)
        return original_process(*args, **kwargs)

    monkeypatch.setattr(logger, "_process_batch", slow_process)

    t0 = time.perf_counter()
    tr.st70_to_etrs89(500000.0, 300000.0, source=1)
    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    logger.close()
    assert elapsed_ms < 20.0


def test_telemetry_queue_saturation_and_dropped_events(tmp_path):
    """Category C: Queue overflow logs warning, increments dropped counter, never raises."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir), max_queue_size=5)

    for _ in range(50):
        logger.record(1, [math.radians(45.0)], [math.radians(25.0)], is_3d=0, source=1)

    diag = logger.get_diagnostics()
    assert diag["dropped_events_count"] > 0
    logger.close()


def test_telemetry_concurrent_multi_threaded_ingestion(tmp_path):
    """Category C: 10 threads concurrently recording transformations are thread-safe."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    def worker():
        for _ in range(50):
            tr.st70_to_etrs89(500000.0, 300000.0, source=2)

    threads = [threading.Thread(target=worker) for _ in range(10)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    logger.close()
    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()
    cur.execute("SELECT count_calls, count_points FROM computation_usage")
    row = cur.fetchone()
    assert row == (500, 500)
    conn.close()


def test_telemetry_batch_and_file_logger_failure_isolation():
    """Category D: Exploding logger during transform_batch does not fail transformations."""

    class ExplodingLogger:
        def record(self, *args, **kwargs):
            raise RuntimeError("Telemetry collector crashed!")

    tr = TransRO(telemetry_logger=ExplodingLogger())
    results = tr.transform_batch([[500000.0, 300000.0], [500100.0, 300100.0]], from_st70=True)
    assert len(results) == 2
    assert results[0][1] is None
    assert results[1][1] is None


def test_telemetry_sqlite_lock_contention_and_retry(tmp_path):
    """Category D: SQLite busy contention retries and succeeds once unlocked."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    db_path = log_dir / "usage.db"

    time.sleep(0.1)

    ext_conn = sqlite3.connect(str(db_path), timeout=0.1)
    ext_conn.execute("BEGIN EXCLUSIVE")

    logger.record(1, [math.radians(45.0)], [math.radians(25.0)], is_3d=0, source=1)
    time.sleep(0.15)

    ext_conn.rollback()
    ext_conn.close()

    logger.close()
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("SELECT count(*) FROM computation_usage")
    assert cur.fetchone()[0] == 1
    conn.close()


def test_telemetry_disk_write_failure_fallback_logging(tmp_path, caplog):
    """Category D: Terminal write failure logs error and fallback payload."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))

    mock_conn = MagicMock()
    mock_conn.executemany.side_effect = sqlite3.OperationalError("disk I/O error")

    with caplog.at_level(logging.ERROR):
        logger._process_batch(mock_conn, [(1, [0.785], [0.436], 0, 1)], datetime.date.today())

    assert "SQLite write failed after 3 attempts" in caplog.text
    assert "Saved unwritten batch" in caplog.text
    assert logger.fallback_path.exists()
    import json

    lines = logger.fallback_path.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 1
    record = json.loads(lines[0])
    assert "disk I/O error" in record["reason"]
    assert len(record["comp_params"]) >= 1
    logger.close()


def test_telemetry_directory_and_database_recreation(tmp_path):
    """Category E: Deleting usage.db before startup recreates schema cleanly."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    logger.close()

    db_path = log_dir / "usage.db"
    assert db_path.exists()
    db_path.unlink()
    assert not db_path.exists()

    logger2 = SqliteUsageLogger(log_dir=str(log_dir))
    logger2.record(1, [math.radians(45.0)], [math.radians(25.0)], is_3d=0, source=1)
    logger2.close()

    assert db_path.exists()
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = {r[0] for r in cur.fetchall()}
    assert "computation_usage" in tables
    assert "spatial_usage" in tables
    conn.close()


def test_telemetry_graceful_shutdown_queue_draining(tmp_path):
    """Category G: close() drains queued items completely before exiting worker thread."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir), max_queue_size=2000)

    lat, lon = [math.radians(45.0)], [math.radians(25.0)]
    for _ in range(1000):
        logger.record(1, lat, lon, is_3d=0, source=1)

    logger.close(timeout=10.0)

    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()
    cur.execute("SELECT count_points FROM computation_usage")
    row = cur.fetchone()
    assert row is not None
    assert row[0] == 1000
    conn.close()


def test_telemetry_post_close_reentrancy(tmp_path):
    """Category G: Calling record() after close() is a safe no-op."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    logger.close()

    logger.record(1, [math.radians(45.0)], [math.radians(25.0)], is_3d=0, source=1)
    assert logger._queue.empty()


def test_buffered_telemetry_collector_chunking():
    """Category C: BufferedTelemetryCollector chunks records and flushes remaining points."""
    mock_target = MagicMock()
    mock_target.get_diagnostics.return_value = {"is_alive": True}

    buf = BufferedTelemetryCollector(mock_target, chunk_size=5)

    # Record 4 points -> not yet flushed
    for i in range(4):
        buf.record(1, 0.785 + i * 0.001, 0.436 + i * 0.001, is_3d=0, source=1)

    assert mock_target.record.call_count == 0

    # Record 5th point -> flushes first chunk of 5 points
    buf.record(1, 0.789, 0.440, is_3d=0, source=1)
    assert mock_target.record.call_count == 1
    call_args = mock_target.record.call_args[0]
    assert call_args[0] == 1  # op
    assert len(call_args[1]) == 5  # lats
    assert len(call_args[2]) == 5  # lons

    # Record 2 more points and flush_all()
    buf.record(1, 0.790, 0.441, is_3d=0, source=1)
    buf.record(1, 0.791, 0.442, is_3d=0, source=1)
    buf.flush_all()

    assert mock_target.record.call_count == 2
    second_call = mock_target.record.call_args[0]
    assert len(second_call[1]) == 2

    assert buf.get_diagnostics() == {"is_alive": True}


def test_telemetry_nan_inf_direct_record_safe_isolation(tmp_path):
    """Category A/C: Direct record() calls with NaN/Inf do not crash worker thread and write 0 rows."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))

    # Directly feed NaN and Inf coordinates to record
    logger.record(1, float("nan"), float("nan"), is_3d=0, source=1)
    logger.record(1, [float("inf"), float("-inf")], [0.43, 0.44], is_3d=0, source=1)
    time.sleep(0.3)

    # Worker thread must still be alive!
    diag = logger.get_diagnostics()
    assert diag["is_alive"] is True

    # Now record a valid coordinate to prove the collector is fully functional
    lat_val, lon_val = math.radians(45.0), math.radians(25.0)
    logger.record(1, [lat_val], [lon_val], is_3d=0, source=1)
    logger.close()

    conn = sqlite3.connect(str(log_dir / "usage.db"))
    cur = conn.cursor()
    cur.execute("SELECT count_calls, count_points FROM computation_usage")
    rows = cur.fetchall()
    # Exactly 1 valid record was written, NaN/Inf calls were safely discarded without crashing
    assert rows == [(1, 1)]
    conn.close()


def test_buffered_telemetry_collector_none_target_safe():
    """Category C: BufferedTelemetryCollector with None target handles record and flush_all safely."""
    buf = BufferedTelemetryCollector(None, chunk_size=5)
    # Should not raise any exception
    buf.record(1, 0.785, 0.436, is_3d=0, source=1)
    buf.record(1, [0.786, 0.787], [0.437, 0.438], is_3d=0, source=1)
    buf.flush_all()
    assert buf.get_diagnostics() == {}


def test_buffered_telemetry_collector_batch_oversized_chunking():
    """Category C: Batches exceeding chunk_size are properly sliced into discrete chunks."""
    mock_target = MagicMock()
    buf = BufferedTelemetryCollector(mock_target, chunk_size=5)

    # Supply 12 points in a single record call
    lats = [0.785 + i * 0.001 for i in range(12)]
    lons = [0.436 + i * 0.001 for i in range(12)]
    buf.record(1, lats, lons, is_3d=0, source=2)

    # Should have emitted 2 chunks of 5 points immediately
    assert mock_target.record.call_count == 2
    c1 = mock_target.record.call_args_list[0][0]
    c2 = mock_target.record.call_args_list[1][0]
    assert len(c1[1]) == 5
    assert len(c2[1]) == 5

    # flush_all emits remaining 2 points
    buf.flush_all()
    assert mock_target.record.call_count == 3
    c3 = mock_target.record.call_args_list[2][0]
    assert len(c3[1]) == 2


def test_telemetry_worker_reconnection_failure_fallback_logged(tmp_path, caplog, monkeypatch):
    """Category D: Failed database reconnection logs error, dropped count, and unwritten payload."""
    log_dir = tmp_path / "logs"

    def failing_connect(self):
        raise sqlite3.OperationalError("Disk unmounted")

    monkeypatch.setattr(SqliteUsageLogger, "_get_db_connection", failing_connect)

    with caplog.at_level(logging.ERROR):
        logger = SqliteUsageLogger(log_dir=str(log_dir))
        time.sleep(0.1)
        logger.record(1, [math.radians(45.0)], [math.radians(25.0)], is_3d=0, source=1)
        time.sleep(0.5)

        diag = logger.get_diagnostics()
        assert diag["dropped_events_count"] >= 1
        assert "SQLite connection unavailable" in caplog.text
        assert "Saved unwritten queue batch" in caplog.text
        assert logger.fallback_path.exists()
        import json

        lines = logger.fallback_path.read_text(encoding="utf-8").strip().splitlines()
        assert len(lines) >= 1
        record = json.loads(lines[0])
        assert "Disk unmounted" in record["reason"]
        assert len(record["raw_batch"]) >= 1
        logger.close()


def test_telemetry_window_sentinel_calculations():
    """Verifies -2 (year), -1 (month), 1 (daily), and N-day window calculations."""
    logger = SqliteUsageLogger.__new__(SqliteUsageLogger)

    # 1. Calendar Year (-2)
    dt1 = datetime.date(2025, 4, 15)
    s, e = logger._calculate_period_window(dt1, -2, "2023-01-01")
    assert s == "2025-01-01"
    assert e == "2025-12-31"

    # 2. Calendar Month (-1)
    dt2 = datetime.date(2024, 2, 10)  # leap year
    s, e = logger._calculate_period_window(dt2, -1, "2023-01-01")
    assert s == "2024-02-01"
    assert e == "2024-02-29"

    # 3. Daily (<= 1)
    dt3 = datetime.date(2026, 10, 24)
    s, e = logger._calculate_period_window(dt3, 1, "2023-01-01")
    assert s == "2026-10-24"
    assert e == "2026-10-24"

    # 4. Mathematical N-day step (10 days from 2026-01-01)
    # Day 0 to 9: 2026-01-01 to 2026-01-10
    s, e = logger._calculate_period_window(datetime.date(2026, 1, 5), 10, "2026-01-01")
    assert s == "2026-01-01"
    assert e == "2026-01-10"

    # Day 15: 2026-01-11 to 2026-01-20
    s, e = logger._calculate_period_window(datetime.date(2026, 1, 16), 10, "2026-01-01")
    assert s == "2026-01-11"
    assert e == "2026-01-20"


def test_telemetry_window_config_transition(tmp_path):
    """Verifies update_window_config clamps open periods and advances cleanly with no overlap."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    conn = logger._get_db_connection()

    # Pre-populate an open monthly period row ending at end of month
    today = datetime.date.today()
    tomorrow = today + datetime.timedelta(days=1)
    month_start = today.strftime("%Y-%m-01")
    last_day = calendar.monthrange(today.year, today.month)[1]
    month_end = f"{today.year:04d}-{today.month:02d}-{last_day:02d}"

    conn.execute(
        """
        INSERT INTO spatial_usage (period_start, period_end, lat_deg, lat_grad, lon_deg, lon_grad, operation, is_3d, source, count_calls, count)
        VALUES (?, ?, 45, 10, 25, 20, 1, 0, 1, 5, 50);
        """,
        (month_start, month_end),
    )
    conn.commit()

    # Transition spatial to 7-day windows starting tomorrow
    result = SqliteUsageLogger.update_window_config(
        conn,
        metric_type="spatial",
        window_days=7,
        valid_from=tomorrow.strftime("%Y-%m-%d"),
        description="Switching to 7-day intervals",
    )

    assert result["valid_from"] == tomorrow.strftime("%Y-%m-%d")
    assert result["cutoff_date"] == today.strftime("%Y-%m-%d")

    # Verify that the open spatial_usage row's period_end was clamped to today
    cursor = conn.cursor()
    cursor.execute("SELECT period_start, period_end FROM spatial_usage")
    row = cursor.fetchone()
    assert row[0] == month_start
    assert row[1] == today.strftime("%Y-%m-%d")  # Clamped to today!

    # Verify telemetry_config has the old config closed and the new config active
    cursor.execute(
        "SELECT window_days, valid_from, valid_to FROM telemetry_config WHERE metric_type = 'spatial' ORDER BY id ASC"
    )
    configs = cursor.fetchall()
    assert len(configs) == 2
    # Old config: valid_to clamped to today
    assert configs[0][2] == today.strftime("%Y-%m-%d")
    # New config: starts tomorrow, valid_to is None
    assert configs[1][0] == 7
    assert configs[1][1] == tomorrow.strftime("%Y-%m-%d")
    assert configs[1][2] is None

    conn.close()
    logger.close()


def test_telemetry_spatial_call_and_point_counts(tmp_path):
    """Verifies that 1 call with multiple points spanning cells increments calls and points properly."""
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    # Call 1: 3 points in Cell A (45.0, 25.0) and 2 points in Cell B (46.0, 26.0)
    lats = [
        math.radians(45.001),
        math.radians(45.002),
        math.radians(45.003),
        math.radians(46.001),
        math.radians(46.002),
    ]
    lons = [
        math.radians(25.001),
        math.radians(25.002),
        math.radians(25.003),
        math.radians(26.001),
        math.radians(26.002),
    ]

    tr.etrs89_to_st70(lats, lons, source=1)

    # Call 2: 1 point in Cell A (45.0, 25.0)
    tr.etrs89_to_st70(math.radians(45.004), math.radians(25.004), source=1)

    logger.close()

    db_path = log_dir / "usage.db"
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Query spatial_usage
    cursor.execute("""
        SELECT lat_deg, lon_deg, count_calls, count
        FROM spatial_usage
        ORDER BY lat_deg ASC
    """)
    rows = cursor.fetchall()
    assert len(rows) == 2

    # Cell A (lat=45, lon=25): touched by Call 1 (3 points) and Call 2 (1 point) -> 2 calls, 4 points
    assert rows[0] == (45, 25, 2, 4)

    # Cell B (lat=46, lon=26): touched only by Call 1 (2 points) -> 1 call, 2 points
    assert rows[1] == (46, 26, 1, 2)

    # Query computation_usage: should have 2 calls and 6 points total
    cursor.execute("SELECT SUM(count_calls), SUM(count_points) FROM computation_usage")
    comp_calls, comp_points = cursor.fetchone()
    assert comp_calls == 2
    assert comp_points == 6

    # Verify points parity: sum of points in spatial_usage == sum of points in computation_usage
    cursor.execute("SELECT SUM(count) FROM spatial_usage")
    spatial_points = cursor.fetchone()[0]
    assert spatial_points == comp_points == 6

    conn.close()
