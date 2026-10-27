"""Unit tests for legacy TransdatOnline statistics endpoints."""

import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from api.dependencies import get_legacy_db_path
from api.main import app

client = TestClient(app)

LEGACY_DB_PATH = Path("logs/legacy_stats.db")
HAS_REAL_LEGACY_DB = LEGACY_DB_PATH.is_file()


def test_legacy_stats_endpoints_when_db_missing(tmp_path):
    """Verify legacy stats endpoints return graceful empty defaults when DB is absent."""
    missing_path = tmp_path / "nonexistent.db"
    app.dependency_overrides[get_legacy_db_path] = lambda: missing_path
    try:
        # 1. Summary stats
        res_stats = client.get("/api/v1/telemetry/legacy/stats")
        assert res_stats.status_code == 200
        data_stats = res_stats.json()
        assert data_stats["total_points"] == 0
        assert data_stats["recorded_days"] == 0
        assert data_stats["from_date"] is None
        assert data_stats["to_date"] is None

        # 2. Timeline
        res_tl = client.get("/api/v1/telemetry/legacy/timeline")
        assert res_tl.status_code == 200
        data_tl = res_tl.json()
        assert data_tl["summary"]["total_points"] == 0
        assert data_tl["series"] == []
        assert data_tl["bounds"] == {}
    finally:
        app.dependency_overrides.pop(get_legacy_db_path, None)


@pytest.fixture
def synthetic_legacy_db(tmp_path):
    """Fixture providing a known SQLite legacy statistics database for hermetic CI tests."""
    db_file = tmp_path / "test_legacy_stats.db"
    conn = sqlite3.connect(str(db_file))
    conn.execute(
        """
        CREATE TABLE statistics (
            id INTEGER PRIMARY KEY,
            day TEXT NOT NULL UNIQUE,
            cooOpAppSt70 INTEGER NOT NULL DEFAULT 0,
            cooOpAppSt30 INTEGER NOT NULL DEFAULT 0,
            cooOpServiceSt70 INTEGER NOT NULL DEFAULT 0,
            cooOpServiceSt30 INTEGER NOT NULL DEFAULT 0
        )
        """
    )
    # Monthly anchor rows for Jan to Sep 2026
    for m in range(1, 10):
        day_str = f"2026-{m:02d}-15"
        conn.execute(
            "INSERT INTO statistics (day, cooOpAppSt70, cooOpAppSt30, cooOpServiceSt70, cooOpServiceSt30) VALUES (?, ?, ?, ?, ?)",
            (day_str, 100, 10, 50, 5),
        )
    # Daily rows for 2026-10-01 through 2026-10-14
    for d in range(1, 15):
        day_str = f"2026-10-{d:02d}"
        conn.execute(
            "INSERT INTO statistics (day, cooOpAppSt70, cooOpAppSt30, cooOpServiceSt70, cooOpServiceSt30) VALUES (?, ?, ?, ?, ?)",
            (day_str, 100, 10, 50, 5),
        )
    conn.commit()
    conn.close()

    app.dependency_overrides[get_legacy_db_path] = lambda: db_file
    yield db_file
    app.dependency_overrides.pop(get_legacy_db_path, None)


def test_synthetic_legacy_timeline_daily(synthetic_legacy_db):
    """Verify daily aggregation and bounds querying using synthetic database."""
    response = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-10-01", "to_date": "2026-10-14", "granularity": "daily"},
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["series"]) == 14
    assert data["bounds"]["min_date"] == "2026-01-15"
    assert data["bounds"]["max_date"] == "2026-10-14"
    first = data["series"][0]
    assert first["points"]["app_st70"] == 100
    assert first["points"]["service_st70"] == 50
    assert first["points"]["total"] == 165


def test_synthetic_legacy_timeline_monthly(synthetic_legacy_db):
    """Verify monthly aggregation across 10 months using synthetic database."""
    response = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-01-01", "to_date": "2026-10-14", "granularity": "monthly"},
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["series"]) == 10
    assert data["series"][-1]["date"] == "2026-10"


def test_synthetic_legacy_timeline_filters(synthetic_legacy_db):
    """Verify source and operation filtering using synthetic database."""
    res_app = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-10-01", "to_date": "2026-10-14", "source": "app"},
    )
    res_srv = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-10-01", "to_date": "2026-10-14", "source": "service"},
    )
    res_all = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-10-01", "to_date": "2026-10-14"},
    )
    assert res_app.status_code == 200
    assert res_srv.status_code == 200
    assert res_all.status_code == 200
    tot_app = res_app.json()["summary"]["total_points"]
    tot_srv = res_srv.json()["summary"]["total_points"]
    tot_all = res_all.json()["summary"]["total_points"]
    assert tot_app + tot_srv == tot_all


@pytest.mark.skipif(
    not HAS_REAL_LEGACY_DB, reason="Historical archive logs/legacy_stats.db not present"
)
def test_legacy_stats_endpoint():
    """Verify GET /api/v1/telemetry/legacy/stats returns accurate totals and date bounds from full archive."""
    response = client.get("/api/v1/telemetry/legacy/stats")
    assert response.status_code == 200
    data = response.json()

    assert "total_points" in data
    assert data["total_points"] == 27_745_527
    assert data["app_st70"] == 21_434_208
    assert data["app_st30"] == 561_718
    assert data["service_st70"] == 5_311_487
    assert data["service_st30"] == 438_114
    assert data["from_date"] == "2011-01-09"
    assert data["to_date"] == "2026-10-14"
    assert data["recorded_days"] == 4977


@pytest.mark.skipif(
    not HAS_REAL_LEGACY_DB, reason="Historical archive logs/legacy_stats.db not present"
)
def test_legacy_timeline_daily():
    """Verify daily aggregation and date range query against full archive."""
    response = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-10-01", "to_date": "2026-10-14", "granularity": "daily"},
    )
    assert response.status_code == 200
    data = response.json()

    assert "summary" in data
    assert "series" in data
    assert "bounds" in data
    assert data["bounds"]["min_date"] == "2011-01-09"
    assert data["bounds"]["max_date"] == "2026-10-14"

    series = data["series"]
    assert len(series) > 0
    first = series[0]
    assert "date" in first
    assert "points" in first
    assert "total" in first["points"]
    assert "app_st70" in first["points"]
    assert "service_st70" in first["points"]


@pytest.mark.skipif(
    not HAS_REAL_LEGACY_DB, reason="Historical archive logs/legacy_stats.db not present"
)
def test_legacy_timeline_monthly():
    """Verify monthly aggregation for a known time span against full archive."""
    response = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-01-01", "to_date": "2026-10-14", "granularity": "monthly"},
    )
    assert response.status_code == 200
    data = response.json()
    series = data["series"]
    assert len(series) == 10
    assert series[-1]["date"] == "2026-10"


@pytest.mark.skipif(
    not HAS_REAL_LEGACY_DB, reason="Historical archive logs/legacy_stats.db not present"
)
def test_legacy_timeline_filters():
    """Verify source and operation filtering against full archive."""
    res_app = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={
            "from_date": "2026-10-01",
            "to_date": "2026-10-14",
            "source": "app",
        },
    )
    assert res_app.status_code == 200
    data_app = res_app.json()

    res_srv = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={
            "from_date": "2026-10-01",
            "to_date": "2026-10-14",
            "source": "service",
        },
    )
    assert res_srv.status_code == 200
    data_srv = res_srv.json()

    res_all = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-10-01", "to_date": "2026-10-14"},
    )
    assert res_all.status_code == 200
    data_all = res_all.json()

    tot_app = data_app["summary"]["total_points"]
    tot_srv = data_srv["summary"]["total_points"]
    tot_all = data_all["summary"]["total_points"]
    assert tot_app + tot_srv == tot_all
