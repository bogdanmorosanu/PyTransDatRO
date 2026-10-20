"""Unit tests for legacy TransdatOnline statistics endpoints."""

from fastapi.testclient import TestClient

from api.main import app

client = TestClient(app)


def test_legacy_stats_endpoint():
    """Verify GET /api/v1/telemetry/legacy/stats returns accurate totals and date bounds."""
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


def test_legacy_timeline_daily():
    """Verify daily aggregation and date range query."""
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


def test_legacy_timeline_monthly():
    """Verify monthly aggregation for a known time span."""
    response = client.get(
        "/api/v1/telemetry/legacy/timeline",
        params={"from_date": "2026-01-01", "to_date": "2026-10-14", "granularity": "monthly"},
    )
    assert response.status_code == 200
    data = response.json()
    series = data["series"]
    assert len(series) == 10
    assert series[-1]["date"] == "2026-10"


def test_legacy_timeline_filters():
    """Verify source and operation filtering."""
    # Filter only app source
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

    # Filter only service source
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

    # Total of app + service should equal unfiltered
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
