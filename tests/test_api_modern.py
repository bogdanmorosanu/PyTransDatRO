"""Tests for modern /api/v1 endpoints, file streaming, health, and info.
"""
import io
import math
import pytest
from starlette.testclient import TestClient
from api.main import app


@pytest.fixture(scope="module")
def client():
    app.state.disable_rate_limiting = True
    with TestClient(app) as test_client:
        yield test_client
    app.state.disable_rate_limiting = False


def test_modern_point_stereo70_to_etrs89_degrees(client):
    """Verifies single point conversion with unit=degrees returns ~45.9997, 24.9984."""
    payload = {
        "op": "Stereo70ToETRS89",
        "coos": [500000.0, 500000.0, 100.0],
        "unit": "degrees"
    }
    res = client.post("/api/v1/transform/point", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["warning"] is None
    lat_deg, lon_deg, h = data["coos"]

    assert math.isclose(lat_deg, 45.9997186, abs_tol=0.0001)
    assert math.isclose(lon_deg, 24.9984476, abs_tol=0.0001)
    assert math.isclose(h, 139.6084, abs_tol=0.01)


def test_modern_point_etrs89_to_stereo70_degrees(client):
    """Verifies inverse single point conversion with unit=degrees."""
    payload = {
        "op": "ETRS89ToStereo70",
        "coos": [45.99971862804511, 24.998447635093715, 139.60844039916992],
        "unit": "degrees"
    }
    res = client.post("/api/v1/transform/point", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["warning"] is None
    n, e, h = data["coos"]
    assert math.isclose(n, 500000.0, abs_tol=0.005)
    assert math.isclose(e, 500000.0, abs_tol=0.005)
    assert math.isclose(h, 100.0, abs_tol=0.005)


def test_modern_point_out_of_grid(client):
    """Verifies out of grid returns coordinates with warning."""
    payload = {
        "op": "Stereo70ToETRS89",
        "coos": [5.0, 1.0, 1.0],
        "unit": "radians"
    }
    res = client.post("/api/v1/transform/point", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["coos"] == [5.0, 1.0, 1.0]
    assert data["warning"] == "Out of grid"


def test_modern_batch_transformation(client):
    """Verifies batch endpoint processes lists of coordinates."""
    payload = {
        "op": "Stereo70ToETRS89",
        "points": [
            [500000.0, 500000.0, 100.0],
            [505000.0, 505000.0, 120.0],
            [5.0, 1.0, 1.0]  # out of grid
        ],
        "unit": "degrees"
    }
    res = client.post("/api/v1/transform/batch", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert data["count"] == 3
    results = data["results"]

    assert results[0]["warning"] is None
    assert math.isclose(results[0]["coos"][0], 45.9997186, abs_tol=0.0001)

    assert results[1]["warning"] is None
    assert results[2]["warning"] == "Out of grid"


def test_modern_file_streaming_csv(client):
    """Verifies file streaming with CSV content, comments, and mixed coordinates."""
    csv_content = (
        "# Topo test points\n"
        "\n"
        "500000.0,500000.0,100.0\n"
        "505000.0,505000.0,120.0\n"
        "5.0,1.0,1.0\n"
    )
    files = {"file": ("test.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    data = {"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"}

    res = client.post("/api/v1/transform/file", files=files, data=data)
    assert res.status_code == 200
    assert "attachment; filename=" in res.headers.get("content-disposition", "")

    lines = res.text.strip().splitlines()
    assert len(lines) == 4
    assert lines[0] == "# Topo test points"
    assert "45.999" in lines[1]
    assert "Out of grid" in lines[3]



def test_modern_file_streaming_semicolon_auto_detect(client):
    """Verifies delimiter auto-detection with semicolon separated file."""
    csv_content = (
        "500000.0;500000.0;100.0\n"
        "505000.0;505000.0;120.0\n"
    )
    files = {"file": ("test.txt", io.BytesIO(csv_content.encode("utf-8")), "text/plain")}
    data = {"op": "Stereo70ToETRS89", "unit": "radians", "delimiter": "auto"}

    res = client.post("/api/v1/transform/file", files=files, data=data)
    assert res.status_code == 200
    lines = res.text.strip().splitlines()
    assert len(lines) == 2
    assert ";" in lines[0]
    assert "0.8028" in lines[0]


def test_health_check(client):
    """Verifies /health probe returns healthy status."""
    res = client.get("/health")
    assert res.status_code == 200
    data = res.json()
    assert data["status"] == "healthy"
    assert data["grid_loaded"] is True
    assert data["telemetry_active"] is True


def test_grid_info(client):
    """Verifies /api/v1/grid/info returns grid boundaries and metadata."""
    res = client.get("/api/v1/grid/info")
    assert res.status_code == 200
    data = res.json()
    assert "rom_grid3d" in data["grid_file"]
    bounds = data["bounds"]
    assert bounds["min_northing"] > 200000
    assert bounds["max_northing"] < 800000
    assert "EPSG:31700" in data["crs"]["projected"]


def test_telemetry_stats(client):
    """Verifies /api/v1/telemetry/stats returns summary and daily trend."""
    res = client.get("/api/v1/telemetry/stats")
    assert res.status_code == 200
    data = res.json()
    assert "summary" in data
    assert "total_calls" in data["summary"]
    assert isinstance(data["daily_trend"], list)
