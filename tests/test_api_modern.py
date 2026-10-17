"""Tests for modern /api/v1 endpoints, file streaming, health, and info."""

import io
import math
from pathlib import Path

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
    payload = {"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0, 100.0], "unit": "degrees"}
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
        "unit": "degrees",
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
    payload = {"op": "Stereo70ToETRS89", "coos": [5.0, 1.0, 1.0], "unit": "radians"}
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
            [5.0, 1.0, 1.0],  # out of grid
        ],
        "unit": "degrees",
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
        "# Topo test points\n\n500000.0,500000.0,100.0\n505000.0,505000.0,120.0\n5.0,1.0,1.0\n"
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
    csv_content = "500000.0;500000.0;100.0\n505000.0;505000.0;120.0\n"
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
    assert "metadata" in data
    assert "params" in data
    assert "grids" in data


def test_grid_metadata_endpoint(client):
    """Verifies /api/v1/grid/metadata returns exhaustive, dynamic SPG metadata."""
    res = client.get("/api/v1/grid/metadata")
    assert res.status_code == 200
    data = res.json()

    # 1. File Info
    assert "file_info" in data
    file_info = data["file_info"]
    assert "rom_grid3d" in file_info["name"]
    assert file_info["size_bytes"] > 100000
    assert len(file_info["sha256"]) == 64

    # 2. Metadata Provenance
    assert "metadata" in data
    meta = data["metadata"]
    assert meta["created_by"] == "CNC"
    assert "license" in meta
    assert "release" in meta
    assert meta["release"]["major"] == 25

    # 3. Parameters & Helmert
    assert "params" in data
    params = data["params"]
    assert "helmert" in params
    assert "os_st70" in params["helmert"]
    assert "st70_os" in params["helmert"]
    assert params["helmert"]["os_st70"]["tE"] > 100.0
    assert "interpolation_descriptions" in params

    # 4. Layers & Diagnostics
    assert "grids" in data
    grids = data["grids"]
    assert "geodetic_shifts" in grids
    assert "geoid_heights" in grids

    # Planimetric layer
    geo_shifts = grids["geodetic_shifts"]
    assert "band_statistics" in geo_shifts
    assert len(geo_shifts["band_statistics"]) == 2
    assert geo_shifts["band_statistics"][0]["valid_nodes"] == 2689
    assert "formatted_extents" in geo_shifts
    assert geo_shifts["formatted_extents"]["easting"]["step_m"] == 11000.0

    # Vertical geoid layer
    geoid = grids["geoid_heights"]
    assert "band_statistics" in geoid
    assert geoid["band_statistics"][0]["valid_nodes"] == 43071
    assert "formatted_extents" in geoid
    assert "min_dms" in geoid["formatted_extents"]["longitude"]


def test_telemetry_stats(client):
    """Verifies /api/v1/telemetry/stats returns summary and daily trend."""
    res = client.get("/api/v1/telemetry/stats")
    assert res.status_code == 200
    data = res.json()
    assert "summary" in data
    assert "total_calls" in data["summary"]
    assert isinstance(data["daily_trend"], list)


def test_telemetry_timeline_endpoint(client):
    """Verifies /api/v1/telemetry/timeline returns structured time series with summary KPIs."""
    res = client.get("/api/v1/telemetry/timeline?granularity=daily")
    assert res.status_code == 200
    data = res.json()
    assert "summary" in data
    assert "total_calls" in data["summary"]
    assert "total_points" in data["summary"]
    assert "avg_points_per_call" in data["summary"]
    assert "series" in data
    assert isinstance(data["series"], list)

    # Weekly granularity
    res_w = client.get("/api/v1/telemetry/timeline?granularity=weekly")
    assert res_w.status_code == 200
    data_w = res_w.json()
    assert "series" in data_w

    # Monthly granularity
    res_m = client.get("/api/v1/telemetry/timeline?granularity=monthly")
    assert res_m.status_code == 200
    data_m = res_m.json()
    assert "series" in data_m

    # Filtered by source
    res_src = client.get("/api/v1/telemetry/timeline?source=1")
    assert res_src.status_code == 200


def test_client_source_header_propagation(client):
    """Verifies X-Client-Source header is recognized by the API."""
    payload = {"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0], "unit": "degrees"}
    # Request from web_app
    headers = {"X-Client-Source": "web_app"}
    res = client.post("/api/v1/transform/point", json=payload, headers=headers)
    assert res.status_code == 200


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


def test_telemetry_spatial_matrix_endpoint_points(client):
    """Verifies /api/v1/telemetry/spatial-matrix returns valid macroscopic grid and quadrants for points."""
    res = client.get("/api/v1/telemetry/spatial-matrix?metric=points")
    assert res.status_code == 200
    data = res.json()
    assert data["resolution"] == 0.1
    assert "bounds" in data
    assert "summary" in data
    assert "cells" in data

    summary = data["summary"]
    assert summary["metric"] == "points"
    assert summary["total_metric_value"] >= 0
    assert "quadrants" in summary

    quads = summary["quadrants"]
    assert quads["center"] == {"lat": 45.9, "lon": 24.9}
    for q in ["NW", "NE", "SW", "SE"]:
        assert q in quads
        assert "value" in quads[q]
        assert "percentage" in quads[q]

    if data["cells"]:
        first_cell = data["cells"][0]
        assert "lat" in first_cell
        assert "lon" in first_cell
        assert "hits" in first_cell
        assert first_cell["hits"] > 0


def test_telemetry_spatial_matrix_endpoint_calls(client):
    """Verifies /api/v1/telemetry/spatial-matrix works with metric=calls."""
    res = client.get("/api/v1/telemetry/spatial-matrix?metric=calls")
    assert res.status_code == 200
    data = res.json()
    assert data["summary"]["metric"] == "calls"


def test_telemetry_spatial_matrix_filters(client):
    """Verifies spatial-matrix filters (dates, sources, operation, 3d)."""
    res = client.get(
        "/api/v1/telemetry/spatial-matrix?from_date=2026-01-01&to_date=2026-12-31&sources=1,2&operation=1&is_3d=1"
    )
    assert res.status_code == 200
    data = res.json()
    assert data["summary"]["from_date"] == "2026-01-01"
    assert data["summary"]["to_date"] == "2026-12-31"


def test_telemetry_spatial_matrix_invalid_metric(client):
    """Verifies error on invalid metric parameter."""
    res = client.get("/api/v1/telemetry/spatial-matrix?metric=invalid")
    assert res.status_code == 422


def test_grid_comparison_endpoint(client):
    """Verifies /api/v1/grid/comparison returns official ANCPI benchmark metrics."""
    res = client.get("/api/v1/grid/comparison")
    assert res.status_code == 200
    data = res.json()
    assert data["total_points"] >= 9000
    assert data["overall_status"] == "PASS"
    assert "etrs89_to_st70" in data
    assert "st70_to_etrs89" in data
    assert data["etrs89_to_st70"]["status"] == "PASS"
    assert data["st70_to_etrs89"]["status"] == "PASS"
    assert "max_abs_error" in data["etrs89_to_st70"]
    assert "max_abs_error" in data["st70_to_etrs89"]


def test_uniform_validation_error_response(client):
    """Verifies that RequestValidationError returns uniform 422 JSON with code='VALIDATION_ERROR'."""
    res = client.post("/api/v1/transform/point", json={"op": "InvalidOp", "coos": [500000.0]})
    assert res.status_code == 422
    data = res.json()
    assert "detail" in data
    assert data.get("code") == "VALIDATION_ERROR"


def test_error_response_schema():
    """Verifies that ErrorResponse Pydantic schema serializes as expected."""
    from api.schemas.error import ErrorResponse

    err = ErrorResponse(detail="Boundary error", code="OUT_OF_BOUNDS")
    assert err.detail == "Boundary error"
    assert err.code == "OUT_OF_BOUNDS"
    dumped = err.model_dump()
    assert dumped == {"detail": "Boundary error", "code": "OUT_OF_BOUNDS"}


def test_modern_point_wrong_dimension_inline_warning(client):
    """Verifies /api/v1/transform/point returns HTTP 200 with 'Wrong dimension' for 1D and 4D inputs."""
    res_1d = client.post(
        "/api/v1/transform/point",
        json={"op": "Stereo70ToETRS89", "coos": [500000.0], "unit": "degrees"},
    )
    assert res_1d.status_code == 200
    assert res_1d.json() == {"coos": [500000.0], "warning": "Wrong dimension"}

    res_4d = client.post(
        "/api/v1/transform/point",
        json={
            "op": "Stereo70ToETRS89",
            "coos": [500000.0, 500000.0, 100.0, 200.0],
            "unit": "degrees",
        },
    )
    assert res_4d.status_code == 200
    assert res_4d.json() == {
        "coos": [500000.0, 500000.0, 100.0, 200.0],
        "warning": "Wrong dimension",
    }


def test_modern_point_and_batch_reject_dms_unit(client):
    """Verifies that unit='dms' is rejected with 422 on numeric /point and /batch endpoints."""
    res_pt = client.post(
        "/api/v1/transform/point",
        json={"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0], "unit": "dms"},
    )
    assert res_pt.status_code == 422

    res_batch = client.post(
        "/api/v1/transform/batch",
        json={"op": "Stereo70ToETRS89", "points": [[500000.0, 500000.0]], "unit": "dms"},
    )
    assert res_batch.status_code == 422


def test_modern_batch_4d_wrong_dimension(client):
    """Verifies /api/v1/transform/batch marks 4D coordinates with 'Wrong dimension' and echoes input."""
    res = client.post(
        "/api/v1/transform/batch",
        json={
            "op": "Stereo70ToETRS89",
            "points": [[500000.0], [500000.0, 500000.0, 100.0, 999.0]],
            "unit": "degrees",
        },
    )
    assert res.status_code == 200
    data = res.json()
    assert data["results"][0] == {"coos": [500000.0], "warning": "Wrong dimension"}
    assert data["results"][1] == {
        "coos": [500000.0, 500000.0, 100.0, 999.0],
        "warning": "Wrong dimension",
    }


def test_modern_file_delimiter_and_column_validation(client):
    """Verifies /api/v1/transform/file rejects invalid delimiter with 422 and 5+ columns with '# Wrong dimension'."""
    csv_content = (
        "P1,500000.0,500000.0,100.0,extra_col\n"
        "P2,45°65'00\"N,25°00'00\"E,100.0\n"
        "P3,500000.0,500000.0,100.0\n"
    )
    # 1. Invalid delimiter -> 422
    files_bad = {"file": ("test.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    res_bad_delim = client.post(
        "/api/v1/transform/file",
        files=files_bad,
        data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "pipe"},
    )
    assert res_bad_delim.status_code == 422

    # 2. 5 columns -> # Wrong dimension
    files_ok = {"file": ("test.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    res_st70 = client.post(
        "/api/v1/transform/file",
        files=files_ok,
        data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"},
    )
    assert res_st70.status_code == 200
    lines = res_st70.text.strip().splitlines()
    assert lines[0].endswith("# Wrong dimension")
    assert lines[1].endswith("# Invalid coordinate data")
    assert lines[2].startswith("P3,45.999")

    # 3. Invalid DMS minutes >= 60 in ETRS89ToStereo70 -> # Invalid coordinate data
    dms_csv = "P1,45°65'00\"N,25°00'00\"E,100.0\n"
    files_dms = {"file": ("dms.csv", io.BytesIO(dms_csv.encode("utf-8")), "text/csv")}
    res_dms = client.post(
        "/api/v1/transform/file",
        files=files_dms,
        data={"op": "ETRS89ToStereo70", "unit": "dms", "delimiter": "comma"},
    )
    assert res_dms.status_code == 200
    assert res_dms.text.strip().endswith("# Invalid coordinate data")


def test_modern_batch_per_point_isolation_option_b(client):
    """Verifies /api/v1/transform/batch isolates null, non-numeric, wrong-dim, and out-of-grid points."""
    res = client.post(
        "/api/v1/transform/batch",
        json={
            "op": "Stereo70ToETRS89",
            "points": [
                [500000.0, 500000.0],
                [500000.0, None],
                ["bad", 500000.0],
                [500000.0],
                [5.0, 1.0],
            ],
            "unit": "degrees",
        },
    )
    assert res.status_code == 200
    results = res.json()["results"]
    assert len(results) == 5
    assert results[0]["warning"] is None and len(results[0]["coos"]) == 2
    assert results[1] == {"coos": [], "warning": "Invalid coordinate data"}
    assert results[2] == {"coos": [], "warning": "Invalid coordinate data"}
    assert results[3] == {"coos": [500000.0], "warning": "Wrong dimension"}
    assert results[4] == {"coos": [5.0, 1.0], "warning": "Out of grid"}


def test_modern_file_unquoted_dms_and_empty_middle_column(client):
    """Verifies unquoted DMS trailing quote preservation and empty middle column rejection."""
    csv_content = (
        'P1,46°00\'00",25°00\'00",100.0\n"P2","46°00\'00""","25°00\'00""","100.0"\nP3,46.0,,100.0\n'
    )
    files = {"file": ("dms_edge.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")}
    res = client.post(
        "/api/v1/transform/file",
        files=files,
        data={"op": "ETRS89ToStereo70", "unit": "degrees", "delimiter": "comma"},
    )
    assert res.status_code == 200
    lines = res.text.strip().splitlines()
    assert len(lines) == 3
    # P1 (unquoted DMS ending with ") must succeed and produce ~500031.2654, ~500120.2236
    assert lines[0].startswith("P1,5000") and "# Invalid" not in lines[0]
    # P3 (empty middle column) must be rejected as Invalid coordinate data, not collapsed to 2D
    assert lines[2].endswith("# Invalid coordinate data")
