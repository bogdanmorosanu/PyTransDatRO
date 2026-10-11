import math
import sqlite3
import pytest
from pathlib import Path
from pytransdatro import TransRO
from pytransdatro.exceptions import OutOfGridErr
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

    cursor.execute("SELECT count(*) FROM usage")
    assert cursor.fetchone()[0] == 0

    cursor.execute("SELECT count(*) FROM daily_stats")
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

    cursor.execute("SELECT count(*) FROM usage")
    assert cursor.fetchone()[0] == 0

    cursor.execute("SELECT count(*) FROM daily_stats")
    assert cursor.fetchone()[0] == 0

    conn.close()


def test_telemetry_in_grid_both_directions_recorded(tmp_path):
    log_dir = tmp_path / "logs"
    logger = SqliteUsageLogger(log_dir=str(log_dir))
    tr = TransRO(telemetry_logger=logger)

    # In-grid ETRS89 to Stereo 70 2D (~45 deg N, ~25 deg E)
    lat = math.radians(45.0)
    lon = math.radians(25.0)
    n_res, e_res = tr.etrs89_to_st70(lat, lon, source=2)
    assert n_res > 0 and e_res > 0

    # In-grid Stereo 70 to ETRS89 2D
    lat_res, lon_res = tr.st70_to_etrs89(500000.0, 300000.0, source=1)
    assert lat_res > 0 and lon_res > 0

    # In-grid 3D ETRS89 to Stereo 70
    n_3d, e_3d, h_3d = tr.etrs89_to_st70(lat, lon, h=150.0, source=2)
    assert n_3d > 0 and e_3d > 0 and h_3d is not None

    # In-grid 3D Stereo 70 to ETRS89
    lat_3d, lon_3d, z_3d = tr.st70_to_etrs89(500000.0, 300000.0, z=150.0, source=1)
    assert lat_3d > 0 and lon_3d > 0 and z_3d is not None

    # In-grid batch call (2 points)
    b_lat, b_lon = tr.st70_to_etrs89([500000.0, 500100.0], [300000.0, 300100.0], source=3)
    assert len(b_lat) == 2 and len(b_lon) == 2

    logger.close()

    db_path = log_dir / "usage.db"
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    cursor.execute("SELECT count(*) FROM usage")
    assert cursor.fetchone()[0] > 0

    # Check that both directions were recorded in usage table
    cursor.execute("SELECT DISTINCT is_st70_to_etrs89 FROM usage ORDER BY is_st70_to_etrs89")
    directions = [row[0] for row in cursor.fetchall()]
    assert directions == [0, 1]  # 0 = etrs89_to_st70, 1 = st70_to_etrs89

    # Check 3D was recorded
    cursor.execute("SELECT DISTINCT is_3d FROM usage ORDER BY is_3d")
    is_3d_values = [row[0] for row in cursor.fetchall()]
    assert is_3d_values == [0, 1]

    # Check daily_stats has records for all sources:
    # source 1: 2 calls (1 2D point + 1 3D point = 2 points)
    # source 2: 2 calls (1 2D point + 1 3D point = 2 points)
    # source 3: 1 call (2 points batch = 2 points)
    cursor.execute("SELECT source, count_calls, count_points FROM daily_stats ORDER BY source")
    daily_stats = cursor.fetchall()
    assert daily_stats == [(1, 2, 2), (2, 2, 2), (3, 1, 2)]

    conn.close()
