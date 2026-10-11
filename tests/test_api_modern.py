"""Tests for modern /api/v1 endpoints, file streaming, health, and info.
"""
import io
import math
import pytest
from pathlib import Path
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


def test_webapp_static_mount(client):
    """Verifies the root path / serves the built webapp index.html."""
    res = client.get("/")
    assert res.status_code == 200
    assert "text/html" in res.headers.get("content-type", "")
    assert "PyTransDatRO" in res.text
    assert 'id="app-root"' in res.text


def test_modern_file_streaming_200_points_st70_to_etrs89(client):
    """Verifies streaming transformation of 200 real-world points from Stereo 70 to ETRS89."""
    file_path = Path(__file__).parent / "data" / "st70_to_etrs89_input_200.csv"
    with open(file_path, "rb") as f:
        files = {"file": ("st70_200.csv", f, "text/csv")}
        data = {"op": "Stereo70ToETRS89", "unit": "dms", "delimiter": "comma"}
        res = client.post("/api/v1/transform/file", files=files, data=data)

    assert res.status_code == 200
    lines = res.text.strip().splitlines()
    assert len(lines) == 200
    # First line must preserve point ID '1' and produce DMS without spaces
    first_line = lines[0]
    assert first_line.startswith("1,")
    assert "°" in first_line
    assert "'" in first_line
    assert '"' in first_line
    assert " # Invalid" not in first_line
    assert " # Out of grid" not in first_line


def test_modern_file_streaming_200_points_etrs89_to_st70_dms(client):
    """Verifies streaming transformation of 200 real-world points with DMS inputs from ETRS89 to Stereo 70."""
    file_path = Path(__file__).parent / "data" / "etrs89_to_st70_input_200.csv"
    with open(file_path, "rb") as f:
        files = {"file": ("etrs89_200.csv", f, "text/csv")}
        data = {"op": "ETRS89ToStereo70", "unit": "degrees", "delimiter": "comma"}
        res = client.post("/api/v1/transform/file", files=files, data=data)

    assert res.status_code == 200
    lines = res.text.strip().splitlines()
    assert len(lines) == 200
    # First line must preserve point ID '1' and produce Stereo 70 coordinates
    first_line = lines[0]
    assert first_line.startswith("1,")
    tokens = first_line.split(",")
    assert len(tokens) >= 3
    n = float(tokens[1])
    e = float(tokens[2])
    assert math.isclose(n, 230134.564, abs_tol=0.01)
    assert math.isclose(e, 434283.040, abs_tol=0.01)
    assert " # Invalid" not in first_line

