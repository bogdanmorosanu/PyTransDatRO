"""Tests for modern /api/v1 endpoints, file streaming, health, and info."""

import io
import math
from pathlib import Path

import pytest
from starlette.testclient import TestClient

import pytransdatro
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
    assert data["engine_version"] == pytransdatro.__version__
    assert app.version == pytransdatro.__version__
    assert data["telemetry_active"] is True
    assert data["telemetry_status"] == "healthy"


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


def test_modern_point_complete_engine_to_api_matrix(client):
    """Verifies /api/v1/transform/point across all engine exceptions and both directions."""
    # 1. Stereo70ToETRS89: Valid 2D & 3D in radians
    res_st70_2d = client.post(
        "/api/v1/transform/point",
        json={"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0], "unit": "radians"},
    )
    assert res_st70_2d.status_code == 200
    d_2d = res_st70_2d.json()
    assert d_2d["warning"] is None and len(d_2d["coos"]) == 2
    assert math.isclose(d_2d["coos"][0], 0.8028465450500996, abs_tol=1e-12)

    # 2. Stereo70ToETRS89: OutOfGrid1DErr (inside 2D grid <=783000, outside 1D geoid >776510)
    res_st70_out1d = client.post(
        "/api/v1/transform/point",
        json={
            "op": "Stereo70ToETRS89",
            "coos": [778000.0, 500000.0, 100.0],
            "unit": "degrees",
        },
    )
    assert res_st70_out1d.status_code == 200
    assert res_st70_out1d.json() == {
        "coos": [778000.0, 500000.0, 100.0],
        "warning": "Out of grid",
    }

    # 3. Stereo70ToETRS89: Pre-check passes in Stereo70, final Grid1D fails (<43.40°N)
    res_st70_wedge = client.post(
        "/api/v1/transform/point",
        json={
            "op": "Stereo70ToETRS89",
            "coos": [205000.0, 260000.0, 100.0],
            "unit": "degrees",
        },
    )
    assert res_st70_wedge.status_code == 200
    assert res_st70_wedge.json() == {
        "coos": [205000.0, 260000.0, 100.0],
        "warning": "Out of grid",
    }

    # 4. Stereo70ToETRS89: NoDataGridErr (2D and 3D)
    for coos in [[549364.0, 147553.0], [549364.0, 147553.0, 100.0]]:
        res_nd = client.post(
            "/api/v1/transform/point",
            json={"op": "Stereo70ToETRS89", "coos": coos, "unit": "degrees"},
        )
        assert res_nd.status_code == 200
        assert res_nd.json() == {"coos": coos, "warning": "No data on grid"}

    # 5. ETRS89ToStereo70: Valid 2D in radians
    res_etrs_2d = client.post(
        "/api/v1/transform/point",
        json={
            "op": "ETRS89ToStereo70",
            "coos": [0.8028465450500996, 0.43630521911977493],
            "unit": "radians",
        },
    )
    assert res_etrs_2d.status_code == 200
    assert res_etrs_2d.json()["warning"] is None
    assert math.isclose(res_etrs_2d.json()["coos"][0], 500000.0, abs_tol=0.005)

    # 6. ETRS89ToStereo70: OutOfGrid2DErr in degrees echoes input degrees
    res_etrs_oog = client.post(
        "/api/v1/transform/point",
        json={"op": "ETRS89ToStereo70", "coos": [10.0, 10.0], "unit": "degrees"},
    )
    assert res_etrs_oog.status_code == 200
    assert res_etrs_oog.json() == {"coos": [10.0, 10.0], "warning": "Out of grid"}

    # 7. ETRS89ToStereo70: OutOfGrid1DErr (48.45°N, 25.00°E is inside Grid2D but >48.40°N Grid1D)
    res_etrs_out1d = client.post(
        "/api/v1/transform/point",
        json={"op": "ETRS89ToStereo70", "coos": [48.45, 25.00, 100.0], "unit": "degrees"},
    )
    assert res_etrs_out1d.status_code == 200
    assert res_etrs_out1d.json() == {
        "coos": [48.45, 25.00, 100.0],
        "warning": "Out of grid",
    }

    # 8. ETRS89ToStereo70: Curved-border SW wedge (passes 6-point precheck, fails final Grid2D)
    res_etrs_wedge = client.post(
        "/api/v1/transform/point",
        json={"op": "ETRS89ToStereo70", "coos": [43.50, 19.90], "unit": "degrees"},
    )
    assert res_etrs_wedge.status_code == 200
    assert res_etrs_wedge.json() == {"coos": [43.50, 19.90], "warning": "Out of grid"}

    # 9. ETRS89ToStereo70: Northern parallel bulge (48.50°N, 25.00°E)
    res_etrs_bulge = client.post(
        "/api/v1/transform/point",
        json={"op": "ETRS89ToStereo70", "coos": [48.50, 25.00], "unit": "degrees"},
    )
    assert res_etrs_bulge.status_code == 200
    assert res_etrs_bulge.json() == {"coos": [48.50, 25.00], "warning": "Out of grid"}

    # 10. ETRS89ToStereo70: NoDataGridErr (2D and 3D in degrees)
    lat_nd_deg = math.degrees(0.8087707898617313)
    lon_nd_deg = math.degrees(0.35742127559921274)
    for coos in [[lat_nd_deg, lon_nd_deg], [lat_nd_deg, lon_nd_deg, 100.0]]:
        res_etrs_nd = client.post(
            "/api/v1/transform/point",
            json={"op": "ETRS89ToStereo70", "coos": coos, "unit": "degrees"},
        )
        assert res_etrs_nd.status_code == 200
        assert res_etrs_nd.json() == {"coos": coos, "warning": "No data on grid"}


def test_modern_point_invalid_coordinate_elements_and_0d(client):
    """Verifies /api/v1/transform/point rejects null, bool, non-numeric string, NaN/Inf, and 0D cleanly."""
    from api.schemas.point import AngleUnitEnum, OperationEnum
    from api.services.transformer import execute_point_transformation, format_modern_response

    # 0D empty array -> Wrong dimension
    res_0d = client.post(
        "/api/v1/transform/point",
        json={"op": "Stereo70ToETRS89", "coos": [], "unit": "degrees"},
    )
    assert res_0d.status_code == 200
    assert res_0d.json() == {"coos": [], "warning": "Wrong dimension"}

    # Null, bool, and non-numeric string -> Invalid coordinate data
    for bad_coos in [
        [None, 500000.0],
        [True, 500000.0],
        [False, 500000.0],
        ["bad", 500000.0],
        ["NaN", 500000.0],
        ["Infinity", 500000.0],
    ]:
        res_bad = client.post(
            "/api/v1/transform/point",
            json={"op": "Stereo70ToETRS89", "coos": bad_coos, "unit": "degrees"},
        )
        assert res_bad.status_code == 200
        assert res_bad.json() == {"coos": [], "warning": "Invalid coordinate data"}

    # Direct float NaN and Inf check on execute_point_transformation
    for non_finite in [[float("nan"), 500000.0], [500000.0, float("inf")]]:
        out_coos, warn = execute_point_transformation(
            app.state.trans,
            OperationEnum.STEREO70_TO_ETRS89.value,
            non_finite,
            unit=AngleUnitEnum.DEGREES,
        )
        assert format_modern_response(out_coos, warn).model_dump() == {
            "coos": [],
            "warning": "Invalid coordinate data",
        }


def test_modern_batch_complete_matrix_both_directions(client):
    """Verifies /api/v1/transform/batch across all warning types in both directions and max_batch_points."""
    from api.config import Settings, get_settings

    # 1. ETRS89ToStereo70 batch covering valid, out of grid (2D, 1D, wedge), no data, wrong dim, invalid
    lat_nd_deg = math.degrees(0.8087707898617313)
    lon_nd_deg = math.degrees(0.35742127559921274)
    res_etrs = client.post(
        "/api/v1/transform/batch",
        json={
            "op": "ETRS89ToStereo70",
            "unit": "degrees",
            "points": [
                [45.99971862804511, 24.998447635093715],  # 0: valid 2D
                [45.99971862804511, 24.998447635093715, 139.60844039916992],  # 1: valid 3D
                [10.0, 10.0],  # 2: out of grid 2D
                [48.45, 25.00, 100.0],  # 3: out of grid 1D geoid
                [43.50, 19.90],  # 4: SW curved-border wedge
                [lat_nd_deg, lon_nd_deg],  # 5: no data on grid 2D
                [lat_nd_deg, lon_nd_deg, 100.0],  # 6: no data on grid 3D
                [],  # 7: wrong dimension 0D
                [True, 25.0],  # 8: invalid (bool)
                ["bad", "25.0"],  # 9: invalid (non-numeric string)
            ],
        },
    )
    assert res_etrs.status_code == 200
    results = res_etrs.json()["results"]
    assert len(results) == 10
    assert results[0]["warning"] is None and len(results[0]["coos"]) == 2
    assert results[1]["warning"] is None and len(results[1]["coos"]) == 3
    assert results[2] == {"coos": [10.0, 10.0], "warning": "Out of grid"}
    assert results[3] == {"coos": [48.45, 25.00, 100.0], "warning": "Out of grid"}
    assert results[4] == {"coos": [43.50, 19.90], "warning": "Out of grid"}
    assert results[5] == {"coos": [lat_nd_deg, lon_nd_deg], "warning": "No data on grid"}
    assert results[6] == {"coos": [lat_nd_deg, lon_nd_deg, 100.0], "warning": "No data on grid"}
    assert results[7] == {"coos": [], "warning": "Wrong dimension"}
    assert results[8] == {"coos": [], "warning": "Invalid coordinate data"}
    assert results[9] == {"coos": [], "warning": "Invalid coordinate data"}

    # 2. Exceeding max_batch_points returns 422
    app.dependency_overrides[get_settings] = lambda: Settings(max_batch_points=2)
    try:
        res_limit = client.post(
            "/api/v1/transform/batch",
            json={
                "op": "Stereo70ToETRS89",
                "points": [[500000.0, 500000.0], [500000.0, 500000.0], [500000.0, 500000.0]],
            },
        )
        assert res_limit.status_code == 422
        assert "exceeds maximum limit" in res_limit.json()["detail"].lower()
    finally:
        app.dependency_overrides.pop(get_settings, None)


def test_modern_file_complete_matrix_and_non_finite_guards(client):
    """Verifies /api/v1/transform/file for No data, 1D/wedge Out of grid, and 3-token -Inf/Inf/NaN guards."""
    # 1. Stereo70ToETRS89 file with No data, OutOfGrid1D, Wedge, and 3-token -Inf/Inf/NaN
    st70_csv = (
        "P_OK,500000.0,500000.0,100.0\n"
        "P_ND2D,549364.0,147553.0\n"
        "P_ND3D,549364.0,147553.0,100.0\n"
        "P_OOG1D,778000.0,500000.0,100.0\n"
        "P_WEDGE,205000.0,260000.0,100.0\n"
        "-Inf,500000.0,500000.0\n"
        "NaN,500000.0,500000.0\n"
        "500000.0,500000.0,Inf\n"
    )
    res_st70 = client.post(
        "/api/v1/transform/file",
        files={"file": ("st70.csv", io.BytesIO(st70_csv.encode("utf-8")), "text/csv")},
        data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"},
    )
    assert res_st70.status_code == 200
    st70_lines = res_st70.text.strip().splitlines()
    assert len(st70_lines) == 8
    assert st70_lines[0].startswith("P_OK,45.999")
    assert st70_lines[1].endswith("# No data on grid")
    assert st70_lines[2].endswith("# No data on grid")
    assert st70_lines[3].endswith("# Out of grid")
    assert st70_lines[4].endswith("# Out of grid")
    assert st70_lines[5] == "-Inf,500000.0,500000.0 # Invalid coordinate data"
    assert st70_lines[6] == "NaN,500000.0,500000.0 # Invalid coordinate data"
    assert st70_lines[7] == "500000.0,500000.0,Inf # Invalid coordinate data"

    # 2. ETRS89ToStereo70 file with No data, OutOfGrid1D, SW wedge, North bulge, and 3-token Inf
    lat_nd_deg = math.degrees(0.8087707898617313)
    lon_nd_deg = math.degrees(0.35742127559921274)
    etrs_csv = (
        f"P_ND,{lat_nd_deg},{lon_nd_deg},100.0\n"
        "P_OOG1D,48.45,25.00,100.0\n"
        "P_WEDGE,43.50,19.90\n"
        "P_BULGE,48.50,25.00\n"
        "Inf,46.0,25.0\n"
        "-Inf,46.0,25.0\n"
    )
    res_etrs = client.post(
        "/api/v1/transform/file",
        files={"file": ("etrs.csv", io.BytesIO(etrs_csv.encode("utf-8")), "text/csv")},
        data={"op": "ETRS89ToStereo70", "unit": "degrees", "delimiter": "comma"},
    )
    assert res_etrs.status_code == 200
    etrs_lines = res_etrs.text.strip().splitlines()
    assert len(etrs_lines) == 6
    assert etrs_lines[0].endswith("# No data on grid")
    assert etrs_lines[1].endswith("# Out of grid")
    assert etrs_lines[2].endswith("# Out of grid")
    assert etrs_lines[3].endswith("# Out of grid")
    assert etrs_lines[4] == "Inf,46.0,25.0 # Invalid coordinate data"
    assert etrs_lines[5] == "-Inf,46.0,25.0 # Invalid coordinate data"

    # 3. Tab and space delimiters
    tab_csv = "P1\t500000.0\t500000.0\t100.0\n"
    res_tab = client.post(
        "/api/v1/transform/file",
        files={"file": ("tab.txt", io.BytesIO(tab_csv.encode("utf-8")), "text/plain")},
        data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "tab"},
    )
    assert res_tab.status_code == 200
    assert res_tab.text.strip().startswith("P1\t45.999")

    space_csv = "P1 500000.0 500000.0 100.0\n"
    res_space = client.post(
        "/api/v1/transform/file",
        files={"file": ("space.txt", io.BytesIO(space_csv.encode("utf-8")), "text/plain")},
        data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "space"},
    )
    assert res_space.status_code == 200
    assert res_space.text.strip().startswith("P1 45.999")


def test_modern_file_upload_size_and_row_limits(client):
    """Verifies HTTP 413 on upload size limit and inline truncation when max_file_points is exceeded."""
    from api.config import Settings, get_settings

    csv_content = "500000.0,500000.0,100.0\n505000.0,505000.0,120.0\n510000.0,510000.0,130.0\n"

    # 1. Exceeding max_file_upload_bytes returns 413 Payload Too Large
    app.dependency_overrides[get_settings] = lambda: Settings(max_file_upload_bytes=10)
    try:
        res_413 = client.post(
            "/api/v1/transform/file",
            files={"file": ("big.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")},
            data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"},
        )
        assert res_413.status_code == 413
        assert "exceeds maximum allowed size" in res_413.json()["detail"].lower()
    finally:
        app.dependency_overrides.pop(get_settings, None)

    # 2. Exceeding max_file_points truncates stream with '# Limit exceeded'
    app.dependency_overrides[get_settings] = lambda: Settings(max_file_points=2)
    try:
        res_rows = client.post(
            "/api/v1/transform/file",
            files={"file": ("rows.csv", io.BytesIO(csv_content.encode("utf-8")), "text/csv")},
            data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"},
        )
        assert res_rows.status_code == 200
        lines = res_rows.text.strip().splitlines()
        assert len(lines) == 3
        assert "45.999" in lines[0]
        assert lines[2] == "# Limit exceeded: maximum 2 points allowed"
    finally:
        app.dependency_overrides.pop(get_settings, None)


def test_modern_file_skip_header_option(client):
    """Verifies skip_header=True passes the first line verbatim without error comment or warning."""
    csv_with_header = (
        "PointID,StereoN,StereoE,Height\nP1,500000.0,500000.0,100.0\nP2,500000.0,500000.0,120.0\n"
    )

    # 1. With skip_header=False (default): header line is parsed as coordinates and flagged
    files_no_skip = {
        "file": ("header.csv", io.BytesIO(csv_with_header.encode("utf-8")), "text/csv")
    }
    res_no_skip = client.post(
        "/api/v1/transform/file",
        files=files_no_skip,
        data={
            "op": "Stereo70ToETRS89",
            "unit": "degrees",
            "delimiter": "comma",
            "skip_header": "false",
        },
    )
    assert res_no_skip.status_code == 200
    lines_no_skip = res_no_skip.text.strip().splitlines()
    assert len(lines_no_skip) == 3
    assert lines_no_skip[0].endswith("# Invalid coordinate data")
    assert lines_no_skip[1].startswith("P1,45.999")

    # 2. With skip_header=True: header line is preserved verbatim
    files_skip = {"file": ("header.csv", io.BytesIO(csv_with_header.encode("utf-8")), "text/csv")}
    res_skip = client.post(
        "/api/v1/transform/file",
        files=files_skip,
        data={
            "op": "Stereo70ToETRS89",
            "unit": "degrees",
            "delimiter": "comma",
            "skip_header": "true",
        },
    )
    assert res_skip.status_code == 200
    lines_skip = res_skip.text.strip().splitlines()
    assert len(lines_skip) == 3
    assert lines_skip[0] == "PointID,StereoN,StereoE,Height"
    assert lines_skip[1].startswith("P1,45.999")
    assert lines_skip[2].startswith("P2,45.999")


def test_health_check_degraded_and_telemetry_flags(client):
    """Verifies /health reports degraded when trans is None, and tracks telemetry health."""
    from unittest.mock import MagicMock

    orig_trans = app.state.trans
    orig_telemetry = getattr(app.state, "telemetry", None)
    try:
        # 1. Telemetry disabled
        app.state.trans = None
        app.state.telemetry = None
        res = client.get("/health")
        assert res.status_code == 200
        data = res.json()
        assert data["status"] == "degraded"
        assert data["grid_loaded"] is False
        assert data["telemetry_active"] is False
        assert data["telemetry_status"] == "disabled"

        # 2. Worker thread dead
        app.state.trans = orig_trans
        mock_dead = MagicMock()
        mock_dead.get_diagnostics.return_value = {"is_alive": False, "last_error": None}
        app.state.telemetry = mock_dead
        res_dead = client.get("/health")
        assert res_dead.status_code == 200
        d_dead = res_dead.json()
        assert d_dead["status"] == "degraded"
        assert d_dead["telemetry_status"] == "degraded"
        assert d_dead["telemetry_active"] is False
        assert "Worker thread is not running" in d_dead["telemetry_error"]

        # 3. Telemetry error state
        mock_err = MagicMock()
        mock_err.get_diagnostics.return_value = {
            "is_alive": True,
            "last_error": "database disk image is malformed",
        }
        app.state.telemetry = mock_err
        res_err = client.get("/health")
        assert res_err.status_code == 200
        d_err = res_err.json()
        assert d_err["status"] == "degraded"
        assert d_err["telemetry_status"] == "degraded"
        assert d_err["telemetry_active"] is True
        assert d_err["telemetry_error"] == "database disk image is malformed"
    finally:
        app.state.trans = orig_trans
        app.state.telemetry = orig_telemetry


def test_client_source_resolution_all_origins(client):
    """Verifies X-Client-Source maps to 1 (rest_api), 2 (web_map), and 3 (desktop_cli)."""
    from unittest.mock import MagicMock

    orig_collector = app.state.trans._collector
    mock_collector = MagicMock()
    app.state.trans._collector = mock_collector

    try:
        source_expectations = [
            ("", 1),
            ("rest_api", 1),
            ("web_map", 2),
            ("web_app", 2),
            ("desktop_cli", 3),
            ("cli", 3),
        ]
        for header_val, expected_src in source_expectations:
            mock_collector.reset_mock()
            headers = {"X-Client-Source": header_val} if header_val else {}
            res = client.post(
                "/api/v1/transform/point",
                json={"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0], "unit": "degrees"},
                headers=headers,
            )
            assert res.status_code == 200
            assert mock_collector.record.called
            call_kwargs = mock_collector.record.call_args[1]
            assert call_kwargs.get("source") == expected_src
    finally:
        app.state.trans._collector = orig_collector


def test_streaming_file_upload_telemetry_chunking(client):
    """Verifies streaming file upload chunks telemetry points into batches of 500 lines without saturation."""
    from unittest.mock import MagicMock

    orig_collector = app.state.trans._collector
    mock_collector = MagicMock()
    app.state.trans._collector = mock_collector

    try:
        # Create a 1,200 line CSV of valid Stereo 70 coordinates
        lines = ["#ID,N,E"]
        for i in range(1200):
            lines.append(f"P{i:04d},500000.0,300000.0")
        csv_data = "\n".join(lines).encode("utf-8")

        res = client.post(
            "/api/v1/transform/file",
            data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"},
            files={"file": ("coords.csv", io.BytesIO(csv_data), "text/csv")},
            headers={"X-Client-Source": "desktop_cli"},
        )
        assert res.status_code == 200
        output_lines = [ln for ln in res.text.splitlines() if ln.strip() and not ln.startswith("#")]
        assert len(output_lines) == 1200

        # Telemetry should be chunked: 1200 points in chunks of 500 -> 3 calls (500, 500, 200)
        assert mock_collector.record.call_count == 3
        c1 = mock_collector.record.call_args_list[0]
        c2 = mock_collector.record.call_args_list[1]
        c3 = mock_collector.record.call_args_list[2]

        assert len(c1[0][1]) == 500
        assert len(c2[0][1]) == 500
        assert len(c3[0][1]) == 200

        # Verify source routing preserved on chunked calls
        assert c1[1].get("source") == 3
        assert c2[1].get("source") == 3
        assert c3[1].get("source") == 3
    finally:
        app.state.trans._collector = orig_collector


def test_telemetry_endpoints_with_seeded_sqlite_db(client, tmp_path):
    """Verifies /telemetry/stats, /telemetry/timeline, and /telemetry/spatial-matrix against a seeded DB."""
    from api.dependencies import get_db_path
    from pytransdatro.telemetry import SqliteUsageLogger

    # 1. Missing DB path returns empty responses gracefully
    missing_db = tmp_path / "nonexistent" / "usage.db"
    app.dependency_overrides[get_db_path] = lambda: missing_db
    try:
        assert client.get("/api/v1/telemetry/stats").json()["summary"]["total_calls"] == 0
        assert client.get("/api/v1/telemetry/timeline").json()["summary"]["total_calls"] == 0
        assert (
            client.get("/api/v1/telemetry/spatial-matrix").json()["summary"]["total_metric_value"]
            == 0
        )
    finally:
        app.dependency_overrides.pop(get_db_path, None)

    # 2. Seed a temporary SQLite database with known usage and 4-quadrant spatial data
    logger_inst = SqliteUsageLogger(log_dir=str(tmp_path))
    # Flush initial schema creation and close worker
    logger_inst.close()

    seeded_db = tmp_path / "usage.db"
    import sqlite3

    conn = sqlite3.connect(str(seeded_db))
    cur = conn.cursor()
    cur.executemany(
        """
        INSERT INTO computation_usage (period_start, period_end, source, operation, is_3d, count_calls, count_points)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            ("2026-03-01", "2026-03-01", 1, 1, 1, 4, 40),  # rest_api, st70->etrs, 3D
            ("2026-03-02", "2026-03-02", 2, 2, 0, 6, 60),  # web_map, etrs->st70, 2D
            ("2026-03-10", "2026-03-10", 3, 1, 0, 2, 20),  # desktop_cli, st70->etrs, 2D
        ],
    )
    # Seed 4 quadrants around center (lat=45.9, lon=24.9):
    # NW: lat=47.1 (47, 10), lon=23.2 (23, 20) -> 10 points
    # NE: lat=47.1 (47, 10), lon=26.5 (26, 50) -> 20 points
    # SW: lat=44.5 (44, 50), lon=23.2 (23, 20) -> 30 points
    # SE: lat=44.5 (44, 50), lon=26.5 (26, 50) -> 40 points
    cur.executemany(
        """
        INSERT INTO spatial_usage (period_start, period_end, source, operation, is_3d, lat_deg, lat_grad, lon_deg, lon_grad, count)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            ("2026-03-01", "2026-03-01", 1, 1, 1, 47, 10, 23, 20, 10),  # NW
            ("2026-03-01", "2026-03-01", 1, 1, 1, 47, 10, 26, 50, 20),  # NE
            ("2026-03-02", "2026-03-02", 2, 2, 0, 44, 50, 23, 20, 30),  # SW
            ("2026-03-02", "2026-03-02", 2, 2, 0, 44, 50, 26, 50, 40),  # SE
        ],
    )
    conn.commit()
    conn.close()

    app.dependency_overrides[get_db_path] = lambda: seeded_db
    try:
        # Verify /telemetry/stats
        res_stats = client.get("/api/v1/telemetry/stats")
        assert res_stats.status_code == 200
        d_stats = res_stats.json()
        assert d_stats["summary"]["total_calls"] == 12
        assert d_stats["summary"]["total_points"] == 120
        assert len(d_stats["daily_trend"]) == 3

        # Verify /telemetry/timeline with filters and granularities
        for gran in ["daily", "weekly", "monthly"]:
            res_tl = client.get(
                f"/api/v1/telemetry/timeline?from_date=2026-03-01&to_date=2026-03-31&granularity={gran}"
            )
            assert res_tl.status_code == 200
            d_tl = res_tl.json()
            assert d_tl["summary"]["total_calls"] == 12
            assert d_tl["summary"]["total_points"] == 120
            assert d_tl["summary"]["avg_points_per_call"] == 10.0

        # Filter timeline by operation=ETRS89ToStereo70 and is_3d=0
        res_tl_filt = client.get(
            "/api/v1/telemetry/timeline?operation=ETRS89ToStereo70&is_3d=0&source=2"
        )
        assert res_tl_filt.status_code == 200
        assert res_tl_filt.json()["summary"]["total_points"] == 60

        # Verify /telemetry/spatial-matrix quadrant math (total=100 -> 10%, 20%, 30%, 40%)
        res_sm = client.get("/api/v1/telemetry/spatial-matrix?metric=points")
        assert res_sm.status_code == 200
        d_sm = res_sm.json()
        assert d_sm["summary"]["total_metric_value"] == 100
        assert d_sm["summary"]["active_cells_count"] == 4
        assert d_sm["summary"]["max_cell_hits"] == 40
        quads = d_sm["summary"]["quadrants"]
        assert quads["NW"] == {"value": 10, "percentage": 10.0}
        assert quads["NE"] == {"value": 20, "percentage": 20.0}
        assert quads["SW"] == {"value": 30, "percentage": 30.0}
        assert quads["SE"] == {"value": 40, "percentage": 40.0}

        # Verify metric=calls on spatial-matrix (1 row per cell -> 4 total calls, 25% each)
        res_sm_calls = client.get("/api/v1/telemetry/spatial-matrix?metric=calls")
        assert res_sm_calls.status_code == 200
        assert res_sm_calls.json()["summary"]["total_metric_value"] == 4
        assert res_sm_calls.json()["summary"]["quadrants"]["NW"]["percentage"] == 25.0
    finally:
        app.dependency_overrides.pop(get_db_path, None)


def test_telemetry_effective_date_snapping_and_parity(client):
    """Verifies effective date snapping on monthly/multi-day buckets and points parity."""
    # 1. Query a 7-day mid-month window in March 2025 (where spatial is monthly: 2025-03-01 to 2025-03-31)
    res_sm = client.get("/api/v1/telemetry/spatial-matrix?from_date=2025-03-10&to_date=2025-03-17")
    assert res_sm.status_code == 200
    d_sm = res_sm.json()
    assert d_sm["summary"]["from_date"] == "2025-03-10"
    assert d_sm["summary"]["to_date"] == "2025-03-17"
    assert d_sm["summary"]["effective_from_date"] == "2025-03-01"
    assert d_sm["summary"]["effective_to_date"] == "2025-03-31"

    # 2. Query timeline for the same exact full month
    res_tl = client.get("/api/v1/telemetry/timeline?from_date=2025-03-01&to_date=2025-03-31")
    assert res_tl.status_code == 200
    d_tl = res_tl.json()
    assert d_tl["summary"]["effective_from_date"] == "2025-03-01"
    assert d_tl["summary"]["effective_to_date"] == "2025-03-31"

    # 3. Verify points parity: for the full month of March 2025, total_points must match
    month_sm = client.get(
        "/api/v1/telemetry/spatial-matrix?from_date=2025-03-01&to_date=2025-03-31&metric=points"
    )
    assert month_sm.status_code == 200
    month_sm_data = month_sm.json()

    assert month_sm_data["summary"]["total_metric_value"] == d_tl["summary"]["total_points"]
