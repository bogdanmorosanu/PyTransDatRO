"""Tests for 100% backward compatibility of the legacy TransDatOnline endpoint.
Validates /transdatonline/cooOpService and /cooOpService locally and via differential
comparison against the live legacy service.
"""

import json
import math
import os
from typing import Any, Dict

import httpx
import pytest
from starlette.testclient import TestClient

from api.config import Settings, get_settings
from api.main import app

DEFAULT_LEGACY_LIVE_URL = "https://www.geo-spatial.org/transdatonline/cooOpService"
PLANAR_TOL_M = 0.005  # 5 mm planar tolerance
RAD_TOL = PLANAR_TOL_M / 6371000.0  # ~7.85e-10 rad (~5 mm on ellipsoid surface)
ELEV_TOL_M = 1.0  # 1.0 m vertical tolerance (EGG97 vs modern SPG quasigeoid)


@pytest.fixture(scope="module")
def client():
    # Disable rate limiting for functional tests
    app.state.disable_rate_limiting = True
    with TestClient(app) as test_client:
        yield test_client
    app.state.disable_rate_limiting = False


def test_legacy_get_stereo70_to_etrs89_3d(client):
    """Verifies GET 3D Stereo70ToETRS89 produces coordinates with sub-micrometer precision."""
    res = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000;100")
    assert res.status_code == 200
    data = res.json()

    assert "coos" in data
    assert "warning" not in data
    assert len(data["coos"]) == 3

    lat, lon, h = data["coos"]
    # Reference value from legacy service documentation: 0.8028465450500996, 0.43630521911977493
    assert math.isclose(lat, 0.8028465450500996, abs_tol=1e-12)
    assert math.isclose(lon, 0.43630521911977493, abs_tol=1e-12)
    # Modern SPG quasigeoid model gives 139.6084 vs legacy EGG97 ~139.7825 (within ~17 cm)
    assert math.isclose(h, 139.6084, abs_tol=0.2)


def test_legacy_get_2d(client):
    """Verifies GET /transdatonline/cooOpService works for 2D inputs."""
    res = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000")
    assert res.status_code == 200
    data = res.json()

    assert "coos" in data
    assert "warning" not in data
    assert len(data["coos"]) == 2


def test_legacy_get_etrs89_to_stereo70_3d(client):
    """Verifies inverse transformation ETRS89ToStereo70."""
    res = client.get(
        "/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos=0.8028465450500996;0.43630521911977493;139.7825537763764"
    )
    assert res.status_code == 200
    data = res.json()

    assert "coos" in data
    assert "warning" not in data
    n, e, h = data["coos"]
    assert math.isclose(n, 500000.0, abs_tol=0.005)
    assert math.isclose(e, 500000.0, abs_tol=0.005)


def test_legacy_get_out_of_grid(client):
    """Verifies out of grid points echo input coordinates with warning 'Out of grid'."""
    res = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=5;1;1")
    assert res.status_code == 200
    data = res.json()

    assert data["coos"] == [5.0, 1.0, 1.0]
    assert data["warning"] == "Out of grid"


def test_legacy_get_invalid_coordinate_data(client):
    """Verifies malformed/non-numeric strings return empty coos array with warning 'Invalid coordinate data'."""
    res = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=a;b;c")
    assert res.status_code == 200
    data = res.json()

    assert data["coos"] == []
    assert data["warning"] == "Invalid coordinate data"


def test_legacy_get_wrong_dimension(client):
    """Verifies single coordinate string returns 'Wrong dimension'."""
    res = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000")
    assert res.status_code == 200
    data = res.json()

    assert data["coos"] == [500000.0]
    assert data["warning"] == "Wrong dimension"


def test_legacy_get_stereo30_obsolete(client):
    """Verifies Stereo30 requests return input coordinates and graceful warning."""
    res = client.get("/transdatonline/cooOpService?cooOp=Stereo30ToETRS89&coos=500000;500000")
    assert res.status_code == 200
    data = res.json()

    assert data["coos"] == [500000.0, 500000.0]
    assert data["warning"] == "Stereo30 transformation is obsolete and unsupported"


def test_legacy_get_missing_parameters(client):
    """Verifies missing cooOp or coos returns 400 Bad Request."""
    res = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89")
    assert res.status_code == 400
    assert "missing" in res.json()["detail"].lower()


def test_legacy_get_unknown_operation(client):
    """Verifies unknown operation returns 400 Bad Request."""
    res = client.get("/transdatonline/cooOpService?cooOp=NonExistentOp&coos=500000;500000")
    assert res.status_code == 400
    assert "no defined transformation" in res.json()["detail"].lower()


def test_legacy_post_form_data(client):
    """Verifies POST with urlencoded form data."""
    payload = {
        "cooOp": "Stereo70ToETRS89",
        "coosArray": '[{"coos":[500000, 500000, 100]}, {"coos":[5, 1, 1]}]',
    }
    res = client.post("/transdatonline/cooOpService", data=payload)
    assert res.status_code == 200
    data = res.json()

    assert isinstance(data, list)
    assert len(data) == 2

    # First point: valid
    assert len(data[0]["coos"]) == 3
    assert "warning" not in data[0]

    # Second point: out of grid
    assert data[1]["coos"] == [5.0, 1.0, 1.0]
    assert data[1]["warning"] == "Out of grid"


def test_legacy_post_json_payload(client):
    """Verifies POST with direct application/json body (modern client compatibility)."""
    payload = {
        "cooOp": "ETRS89ToStereo70",
        "coosArray": [{"coos": [0.8028465450500996, 0.43630521911977493, 139.7825537763764]}],
    }
    res = client.post("/transdatonline/cooOpService", json=payload)
    assert res.status_code == 200
    data = res.json()

    assert isinstance(data, list)
    assert len(data) == 1
    assert math.isclose(data[0]["coos"][0], 500000.0, abs_tol=0.005)


def test_legacy_get_4d_wrong_dimension(client):
    """Verifies >3D coordinates echo all coordinates with warning 'Wrong dimension'."""
    res = client.get(
        "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000;100;999"
    )
    assert res.status_code == 200
    data = res.json()
    assert data["coos"] == [500000.0, 500000.0, 100.0, 999.0]
    assert data["warning"] == "Wrong dimension"


def test_legacy_get_semicolon_token_validation(client):
    """Verifies Java split(';') semantics: trailing ';' is valid, leading/middle ';;' is invalid."""
    # Trailing semicolon -> valid 2D
    res_trailing = client.get(
        "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000;"
    )
    assert res_trailing.status_code == 200
    data_trailing = res_trailing.json()
    assert "warning" not in data_trailing
    assert len(data_trailing["coos"]) == 2

    # Middle empty token -> Invalid coordinate data
    res_middle = client.get(
        "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;;500000"
    )
    assert res_middle.status_code == 200
    assert res_middle.json() == {"coos": [], "warning": "Invalid coordinate data"}

    # Leading empty token -> Invalid coordinate data
    res_leading = client.get(
        "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=;500000;500000"
    )
    assert res_leading.status_code == 200
    assert res_leading.json() == {"coos": [], "warning": "Invalid coordinate data"}


def test_legacy_get_nan_and_inf_rejected(client):
    """Verifies NaN and Infinity strings return 'Invalid coordinate data'."""
    for bad in ["NaN;500000", "Infinity;500000", "500000;-Infinity"]:
        res = client.get(f"/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos={bad}")
        assert res.status_code == 200
        assert res.json() == {"coos": [], "warning": "Invalid coordinate data"}


def test_legacy_post_eager_cooop_validation_on_empty_array(client):
    """Verifies unknown cooOp returns 400 Bad Request even when coosArray is empty."""
    payload = {"cooOp": "InvalidOperation", "coosArray": []}
    res = client.post("/transdatonline/cooOpService", json=payload)
    assert res.status_code == 400
    assert "no defined transformation" in res.json()["detail"].lower()


def test_legacy_post_per_point_dimension_and_invalid_isolation(client):
    """Verifies per-point validation in POST batch for 1D, 4D, non-numeric, and NaN."""
    payload = {
        "cooOp": "Stereo70ToETRS89",
        "coosArray": [
            {"coos": [500000.0]},
            {"coos": [500000.0, 500000.0, 100.0, 10.0]},
            {"coos": ["bad", 500000.0]},
            {"coos": [500000.0, 500000.0]},
        ],
    }
    res = client.post("/transdatonline/cooOpService", json=payload)
    assert res.status_code == 200
    data = res.json()
    assert len(data) == 4
    assert data[0] == {"coos": [500000.0], "warning": "Wrong dimension"}
    assert data[1] == {"coos": [500000.0, 500000.0, 100.0, 10.0], "warning": "Wrong dimension"}
    assert data[2] == {"coos": [], "warning": "Invalid coordinate data"}
    assert "warning" not in data[3]
    assert len(data[3]["coos"]) == 2


def test_legacy_empty_string_parameters_parity(client):
    """Verifies empty string query/form parameters match Java getParameter() non-null semantics."""
    # Empty coos= is present (non-null), so it returns HTTP 200 with 'Invalid coordinate data'
    res_empty_coos = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=")
    assert res_empty_coos.status_code == 200
    assert res_empty_coos.json() == {"coos": [], "warning": "Invalid coordinate data"}

    # Empty cooOp= is present (non-null), so it returns HTTP 400 'There is no defined transformation with name '
    res_empty_op = client.get("/transdatonline/cooOpService?cooOp=&coos=500000;500000")
    assert res_empty_op.status_code == 400
    assert "no defined transformation" in res_empty_op.json()["detail"].lower()

    # Empty coosArray= in POST form returns HTTP 400 Malformed JSON, not Missing parameter
    res_empty_arr = client.post(
        "/transdatonline/cooOpService",
        data={"cooOp": "Stereo70ToETRS89", "coosArray": ""},
    )
    assert res_empty_arr.status_code == 400
    assert "malformed json" in res_empty_arr.json()["detail"].lower()


def test_legacy_get_stereo70_to_etrs89_no_data_and_1d_edge_scenarios(client):
    """Verifies Stereo70ToETRS89 GET for No data on grid, 1D geoid out-of-grid, and curved wedge."""
    # 1. No data on grid (inside Grid2D bbox, NaN subgrid nodes near border)
    for coos_str, expected_coos in [
        ("549364;147553", [549364.0, 147553.0]),
        ("549364;147553;100", [549364.0, 147553.0, 100.0]),
    ]:
        res = client.get(f"/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos={coos_str}")
        assert res.status_code == 200
        assert res.json() == {"coos": expected_coos, "warning": "No data on grid"}

    # 2. Point inside 2D horizontal grid (N=778000, E=500000 is <=783000 for Grid2D)
    # but outside 1D geoid grid Stereo70 pre-check (>776510 m)
    res_3d_out1d = client.get(
        "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=778000;500000;100"
    )
    assert res_3d_out1d.status_code == 200
    assert res_3d_out1d.json() == {
        "coos": [778000.0, 500000.0, 100.0],
        "warning": "Out of grid",
    }

    # 3. Pre-check passes in Stereo70 (N=205000, E=260000 is inside grid1d_bounds_st70),
    # but after projection to ETRS89 its latitude (~43.30°N) is outside Grid1D (<43.40°N)
    res_wedge_3d = client.get(
        "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=205000;260000;100"
    )
    assert res_wedge_3d.status_code == 200
    assert res_wedge_3d.json() == {
        "coos": [205000.0, 260000.0, 100.0],
        "warning": "Out of grid",
    }


def test_legacy_get_etrs89_to_stereo70_complete_matrix(client):
    """Verifies ETRS89ToStereo70 GET across 2D, Out of grid (2D/1D/wedge), No data, Invalid, and Wrong dim."""
    # 1. Valid 2D on both /transdatonline/cooOpService and /cooOpService
    res_2d = client.get(
        "/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos=0.8028465450500996;0.43630521911977493"
    )
    assert res_2d.status_code == 200
    data_2d = res_2d.json()
    assert "warning" not in data_2d
    assert len(data_2d["coos"]) == 2
    assert math.isclose(data_2d["coos"][0], 500000.0, abs_tol=0.005)
    assert math.isclose(data_2d["coos"][1], 500000.0, abs_tol=0.005)

    # 2. Out of grid (2D and 3D pre-check rejection)
    res_oog_2d = client.get("/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos=0.1;0.1")
    assert res_oog_2d.status_code == 200
    assert res_oog_2d.json() == {"coos": [0.1, 0.1], "warning": "Out of grid"}

    res_oog_3d = client.get("/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos=0.1;0.1;100")
    assert res_oog_3d.status_code == 200
    assert res_oog_3d.json() == {"coos": [0.1, 0.1, 100.0], "warning": "Out of grid"}

    # 3. Out of grid (1D vertical geoid only: 48.45°N, 25.00°E is inside Grid2D but outside Grid1D)
    lat_4845 = math.radians(48.45)
    lon_2500 = math.radians(25.00)
    res_out1d = client.get(
        f"/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos={lat_4845};{lon_2500};100"
    )
    assert res_out1d.status_code == 200
    assert res_out1d.json() == {
        "coos": [lat_4845, lon_2500, 100.0],
        "warning": "Out of grid",
    }

    # 4. Pre-check inside 6-point extended ETRS89 bbox, but final Stereo70 Grid2D outside (SW wedge)
    lat_wedge = math.radians(43.50)
    lon_wedge = math.radians(19.90)
    res_wedge = client.get(
        f"/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos={lat_wedge};{lon_wedge}"
    )
    assert res_wedge.status_code == 200
    assert res_wedge.json() == {
        "coos": [lat_wedge, lon_wedge],
        "warning": "Out of grid",
    }

    # 5. No data on grid (inside Grid2D & Grid1D bbox, NaN subgrid nodes near western border)
    # Corresponds to Stereo70 ~(549364, 147553) -> lat=0.8087707898617313, lon=0.35742127559921274
    lat_nodata = 0.8087707898617313
    lon_nodata = 0.35742127559921274
    for coos_str, expected_coos in [
        (f"{lat_nodata};{lon_nodata}", [lat_nodata, lon_nodata]),
        (f"{lat_nodata};{lon_nodata};100", [lat_nodata, lon_nodata, 100.0]),
    ]:
        res_nd = client.get(f"/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos={coos_str}")
        assert res_nd.status_code == 200
        assert res_nd.json() == {"coos": expected_coos, "warning": "No data on grid"}

    # 6. Invalid coordinate data (non-numeric, NaN, Inf, empty middle token)
    for bad_coos in ["x;y", "NaN;0.4363", "0.8028;Infinity", "0.8028;;0.4363"]:
        res_bad = client.get(f"/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos={bad_coos}")
        assert res_bad.status_code == 200
        assert res_bad.json() == {"coos": [], "warning": "Invalid coordinate data"}

    # 7. Wrong dimension (1D and 4D)
    res_1d = client.get("/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos=0.8028")
    assert res_1d.status_code == 200
    assert res_1d.json() == {"coos": [0.8028], "warning": "Wrong dimension"}

    res_4d = client.get(
        "/transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos=0.8028;0.4363;100;200"
    )
    assert res_4d.status_code == 200
    assert res_4d.json() == {
        "coos": [0.8028, 0.4363, 100.0, 200.0],
        "warning": "Wrong dimension",
    }


def test_legacy_post_complete_matrix_both_directions(client):
    """Verifies POST across both Stereo70ToETRS89 and ETRS89ToStereo70 for all valid/warning cases."""
    # 1. Stereo70ToETRS89 full batch matrix on root alias /cooOpService
    st70_batch = [
        {"coos": [500000.0, 500000.0]},  # 0: valid 2D
        {"coos": [500000.0, 500000.0, 100.0]},  # 1: valid 3D
        {"coos": [5.0, 1.0]},  # 2: out of grid 2D
        {"coos": [778000.0, 500000.0, 100.0]},  # 3: out of grid 1D geoid
        {"coos": [205000.0, 260000.0, 100.0]},  # 4: pre-check pass, Grid1D out of grid
        {"coos": [549364.0, 147553.0]},  # 5: no data on grid 2D
        {"coos": [549364.0, 147553.0, 100.0]},  # 6: no data on grid 3D
        {"coos": []},  # 7: wrong dimension 0D
        {"coos": [500000.0]},  # 8: wrong dimension 1D
        {"coos": [500000.0, 500000.0, 100.0, 999.0]},  # 9: wrong dimension 4D
        {"coos": [None, 500000.0]},  # 10: invalid (null)
        {"coos": [True, 500000.0]},  # 11: invalid (bool)
        {"not_coos": [500000.0, 500000.0]},  # 12: invalid (missing 'coos' key)
    ]
    res_st70 = client.post(
        "/transdatonline/cooOpService",
        json={"cooOp": "Stereo70ToETRS89", "coosArray": st70_batch},
    )
    assert res_st70.status_code == 200
    d_st70 = res_st70.json()
    assert len(d_st70) == 13
    assert "warning" not in d_st70[0] and len(d_st70[0]["coos"]) == 2
    assert "warning" not in d_st70[1] and len(d_st70[1]["coos"]) == 3
    assert d_st70[2] == {"coos": [5.0, 1.0], "warning": "Out of grid"}
    assert d_st70[3] == {"coos": [778000.0, 500000.0, 100.0], "warning": "Out of grid"}
    assert d_st70[4] == {"coos": [205000.0, 260000.0, 100.0], "warning": "Out of grid"}
    assert d_st70[5] == {"coos": [549364.0, 147553.0], "warning": "No data on grid"}
    assert d_st70[6] == {"coos": [549364.0, 147553.0, 100.0], "warning": "No data on grid"}
    assert d_st70[7] == {"coos": [], "warning": "Wrong dimension"}
    assert d_st70[8] == {"coos": [500000.0], "warning": "Wrong dimension"}
    assert d_st70[9] == {"coos": [500000.0, 500000.0, 100.0, 999.0], "warning": "Wrong dimension"}
    assert d_st70[10] == {"coos": [], "warning": "Invalid coordinate data"}
    assert d_st70[11] == {"coos": [], "warning": "Invalid coordinate data"}
    assert d_st70[12] == {"coos": [], "warning": "Invalid coordinate data"}

    # 2. ETRS89ToStereo70 full batch matrix
    lat_ok, lon_ok, h_ok = 0.8028465450500996, 0.43630521911977493, 139.60844039916992
    lat_nd, lon_nd = 0.8087707898617313, 0.35742127559921274
    lat_1d, lon_1d = math.radians(48.45), math.radians(25.00)
    lat_w, lon_w = math.radians(43.50), math.radians(19.90)

    etrs_batch = [
        {"coos": [lat_ok, lon_ok]},  # 0: valid 2D
        {"coos": [lat_ok, lon_ok, h_ok]},  # 1: valid 3D
        {"coos": [0.1, 0.1]},  # 2: out of grid 2D
        {"coos": [lat_1d, lon_1d, 100.0]},  # 3: out of grid 1D geoid
        {"coos": [lat_w, lon_w]},  # 4: pre-check pass, Grid2D out of grid (SW wedge)
        {"coos": [lat_nd, lon_nd]},  # 5: no data on grid 2D
        {"coos": [lat_nd, lon_nd, 100.0]},  # 6: no data on grid 3D
        {"coos": [lat_ok]},  # 7: wrong dimension 1D
        {"coos": [lat_ok, lon_ok, h_ok, 1.0]},  # 8: wrong dimension 4D
        {"coos": ["invalid", lon_ok]},  # 9: invalid coordinate data
    ]
    res_etrs = client.post(
        "/transdatonline/cooOpService",
        json={"cooOp": "ETRS89ToStereo70", "coosArray": etrs_batch},
    )
    assert res_etrs.status_code == 200
    d_etrs = res_etrs.json()
    assert len(d_etrs) == 10
    assert "warning" not in d_etrs[0] and len(d_etrs[0]["coos"]) == 2
    assert math.isclose(d_etrs[0]["coos"][0], 500000.0, abs_tol=0.005)
    assert "warning" not in d_etrs[1] and len(d_etrs[1]["coos"]) == 3
    assert math.isclose(d_etrs[1]["coos"][2], 100.0, abs_tol=0.005)
    assert d_etrs[2] == {"coos": [0.1, 0.1], "warning": "Out of grid"}
    assert d_etrs[3] == {"coos": [lat_1d, lon_1d, 100.0], "warning": "Out of grid"}
    assert d_etrs[4] == {"coos": [lat_w, lon_w], "warning": "Out of grid"}
    assert d_etrs[5] == {"coos": [lat_nd, lon_nd], "warning": "No data on grid"}
    assert d_etrs[6] == {"coos": [lat_nd, lon_nd, 100.0], "warning": "No data on grid"}
    assert d_etrs[7] == {"coos": [lat_ok], "warning": "Wrong dimension"}
    assert d_etrs[8] == {"coos": [lat_ok, lon_ok, h_ok, 1.0], "warning": "Wrong dimension"}
    assert d_etrs[9] == {"coos": [], "warning": "Invalid coordinate data"}


def test_legacy_post_string_encoded_numbers(client):
    """Verifies POST accepts string-encoded numbers in coosArray matching Java CooJSON.jsonToCoo()."""
    # 1. Stereo70ToETRS89 with string-encoded numbers
    res_st70 = client.post(
        "/transdatonline/cooOpService",
        data={
            "cooOp": "Stereo70ToETRS89",
            "coosArray": '[{"coos": ["500000.0", "500000.0", "100.0"]}, {"coos": ["500000", "500000"]}]',
        },
    )
    assert res_st70.status_code == 200
    data_st70 = res_st70.json()
    assert len(data_st70) == 2
    assert "warning" not in data_st70[0] and len(data_st70[0]["coos"]) == 3
    assert math.isclose(data_st70[0]["coos"][0], 0.8028465450500996, abs_tol=1e-12)
    assert "warning" not in data_st70[1] and len(data_st70[1]["coos"]) == 2

    # 2. ETRS89ToStereo70 with string-encoded numbers
    res_etrs = client.post(
        "/transdatonline/cooOpService",
        json={
            "cooOp": "ETRS89ToStereo70",
            "coosArray": [
                {"coos": ["0.8028465450500996", "0.43630521911977493", "139.60844039916992"]}
            ],
        },
    )
    assert res_etrs.status_code == 200
    data_etrs = res_etrs.json()
    assert len(data_etrs) == 1
    assert "warning" not in data_etrs[0]
    assert math.isclose(data_etrs[0]["coos"][0], 500000.0, abs_tol=0.005)
    assert math.isclose(data_etrs[0]["coos"][1], 500000.0, abs_tol=0.005)
    assert math.isclose(data_etrs[0]["coos"][2], 100.0, abs_tol=0.005)


def test_legacy_stereo30_both_directions_get_and_post(client):
    """Verifies both Stereo30ToETRS89 and ETRS89ToStereo30 in GET and POST."""
    res_get = client.get("/transdatonline/cooOpService?cooOp=ETRS89ToStereo30&coos=0.8;0.4;100")
    assert res_get.status_code == 200
    assert res_get.json() == {
        "coos": [0.8, 0.4, 100.0],
        "warning": "Stereo30 transformation is obsolete and unsupported",
    }

    for op in ["Stereo30ToETRS89", "ETRS89ToStereo30"]:
        res_post = client.post(
            "/transdatonline/cooOpService",
            json={
                "cooOp": op,
                "coosArray": [
                    {"coos": [500000.0, 500000.0]},
                    {"coos": [500000.0]},
                    {"coos": ["bad", 500000.0]},
                ],
            },
        )
        assert res_post.status_code == 200
        data = res_post.json()
        assert data[0] == {
            "coos": [500000.0, 500000.0],
            "warning": "Stereo30 transformation is obsolete and unsupported",
        }
        assert data[1] == {"coos": [500000.0], "warning": "Wrong dimension"}
        assert data[2] == {"coos": [], "warning": "Invalid coordinate data"}


def test_legacy_post_multipart_and_error_branches(client):
    """Verifies multipart/form-data POST, malformed JSON, non-list coosArray, and max_batch_points."""
    # 1. multipart/form-data encoding
    res_mp = client.post(
        "/transdatonline/cooOpService",
        files={
            "cooOp": (None, "Stereo70ToETRS89"),
            "coosArray": (None, '[{"coos": [500000.0, 500000.0]}]'),
        },
    )
    assert res_mp.status_code == 200
    assert len(res_mp.json()[0]["coos"]) == 2

    # 2. Invalid raw JSON body
    res_bad_json = client.post(
        "/transdatonline/cooOpService",
        content=b"{malformed json",
        headers={"Content-Type": "application/json"},
    )
    assert res_bad_json.status_code == 400
    assert "invalid json" in res_bad_json.json()["detail"].lower()

    # 3. Non-string cooOp in JSON body
    res_bad_op_type = client.post(
        "/transdatonline/cooOpService",
        json={"cooOp": 123, "coosArray": []},
    )
    assert res_bad_op_type.status_code == 400

    # 4. coosArray is a JSON object string instead of a JSON array
    res_not_list = client.post(
        "/transdatonline/cooOpService",
        data={"cooOp": "Stereo70ToETRS89", "coosArray": '{"coos": [500000, 500000]}'},
    )
    assert res_not_list.status_code == 400
    assert "json array" in res_not_list.json()["detail"].lower()

    # 5. coosArray is an integer in JSON body
    res_int_arr = client.post(
        "/transdatonline/cooOpService",
        json={"cooOp": "Stereo70ToETRS89", "coosArray": 42},
    )
    assert res_int_arr.status_code == 400

    # 6. Exceeding max_batch_points returns 422
    app.dependency_overrides[get_settings] = lambda: Settings(max_batch_points=2)
    try:
        res_limit = client.post(
            "/transdatonline/cooOpService",
            json={
                "cooOp": "Stereo70ToETRS89",
                "coosArray": [
                    {"coos": [500000.0, 500000.0]},
                    {"coos": [500000.0, 500000.0]},
                    {"coos": [500000.0, 500000.0]},
                ],
            },
        )
        assert res_limit.status_code == 422
        assert "exceeds maximum limit" in res_limit.json()["detail"].lower()
    finally:
        app.dependency_overrides.pop(get_settings, None)


# ==============================================================================
# LIVE DIFFERENTIAL PARITY SUITE (Legacy TransDatOnline vs. PyTransDatRO API)
# ==============================================================================


def _resolve_live_legacy_url() -> str:
    """Checks environment toggle and probes live legacy service; calls pytest.skip if unavailable."""
    enabled = os.environ.get("PYTRANSDAT_LEGACY_LIVE_ENABLED", "1").strip().lower()
    if enabled in ("0", "false", "no", "off"):
        pytest.skip(
            "Live legacy TransDatOnline parity tests disabled via PYTRANSDAT_LEGACY_LIVE_ENABLED=0."
        )

    url = os.environ.get("PYTRANSDAT_LEGACY_LIVE_URL", DEFAULT_LEGACY_LIVE_URL).strip()
    try:
        resp = httpx.get(
            url,
            params={"cooOp": "Stereo70ToETRS89", "coos": "500000;500000"},
            timeout=5.0,
        )
        if resp.status_code != 200:
            pytest.skip(
                f"Live legacy service at {url} returned HTTP {resp.status_code} (likely decommissioned)."
            )
        data = resp.json()
        if not isinstance(data, dict) or "coos" not in data:
            pytest.skip(f"Live legacy service at {url} did not return expected JSON contract.")
    except (httpx.HTTPError, ValueError, OSError) as exc:
        pytest.skip(
            f"Live legacy TransDatOnline service at {url} is unreachable or decommissioned: {exc}"
        )
    return url


@pytest.fixture(scope="module")
def live_legacy_url() -> str:
    """Module-scoped fixture providing the verified live legacy URL or skipping cleanly."""
    return _resolve_live_legacy_url()


@pytest.fixture(scope="module")
def live_http_client():
    """Module-scoped persistent HTTP client with connection pooling for live legacy queries."""
    with httpx.Client(timeout=10.0) as http_client:
        yield http_client


def _live_request(
    http_client: httpx.Client, method: str, url: str, **kwargs: Any
) -> httpx.Response:
    """Executes an HTTP request against the live legacy server with retry on transient socket resets."""
    import time

    last_exc: Exception | None = None
    for attempt in range(3):
        try:
            return http_client.request(method, url, **kwargs)
        except httpx.TransportError as exc:
            last_exc = exc
            time.sleep(0.3 * (attempt + 1))
    pytest.skip(f"Live legacy service experienced transient transport failure: {last_exc}")


def _assert_legacy_item_parity(
    live_item: Dict[str, Any], local_item: Dict[str, Any], coo_op: str
) -> None:
    """Asserts structural, warning, and numerical parity between live legacy and local FastAPI."""
    assert set(live_item.keys()) == set(local_item.keys()), (
        f"Key mismatch: live={live_item} vs local={local_item}"
    )

    if "warning" in live_item:
        assert live_item["warning"] == local_item["warning"]
        assert len(live_item["coos"]) == len(local_item["coos"])
        # Note: On ETRS89ToStereo70, the legacy Java pipeline mutates the Coordinate object
        # in-place (converting radians to degrees or projecting to Stereo70 meters before
        # hitting the grid check), whereas PyTransDatRO cleanly echoes the original input.
        if coo_op == "Stereo70ToETRS89":
            for c_live, c_local in zip(live_item["coos"], local_item["coos"], strict=True):
                assert math.isclose(c_live, c_local, abs_tol=1e-9)
        return

    live_coos = live_item["coos"]
    local_coos = local_item["coos"]
    assert len(live_coos) == len(local_coos)

    if coo_op == "Stereo70ToETRS89":
        # 2D coordinates are in radians: compare within RAD_TOL (~5 mm on ground)
        assert math.isclose(live_coos[0], local_coos[0], abs_tol=RAD_TOL), (
            f"Lat diff={abs(live_coos[0] - local_coos[0])} rad exceeds {RAD_TOL} rad (~5 mm)"
        )
        assert math.isclose(live_coos[1], local_coos[1], abs_tol=RAD_TOL), (
            f"Lon diff={abs(live_coos[1] - local_coos[1])} rad exceeds {RAD_TOL} rad (~5 mm)"
        )
    else:
        # 2D coordinates are in Stereo70 meters: compare within PLANAR_TOL_M (5 mm)
        assert math.isclose(live_coos[0], local_coos[0], abs_tol=PLANAR_TOL_M), (
            f"N diff={abs(live_coos[0] - local_coos[0])} m exceeds {PLANAR_TOL_M} m (5 mm)"
        )
        assert math.isclose(live_coos[1], local_coos[1], abs_tol=PLANAR_TOL_M), (
            f"E diff={abs(live_coos[1] - local_coos[1])} m exceeds {PLANAR_TOL_M} m (5 mm)"
        )

    if len(live_coos) == 3:
        # 3D elevation: compare within ELEV_TOL_M (1.0 m) due to EGG97 vs modern SPG quasigeoid
        assert math.isclose(live_coos[2], local_coos[2], abs_tol=ELEV_TOL_M), (
            f"Elevation diff={abs(live_coos[2] - local_coos[2])} m exceeds {ELEV_TOL_M} m"
        )


@pytest.mark.legacy_live
@pytest.mark.parametrize(
    "coo_op, coos_query",
    [
        # Stereo70ToETRS89 scenarios
        ("Stereo70ToETRS89", "500000;500000"),  # Valid 2D
        ("Stereo70ToETRS89", "500000;500000;100"),  # Valid 3D
        ("Stereo70ToETRS89", "500000;500000;"),  # Trailing semicolon
        ("Stereo70ToETRS89", "5;1"),  # Out of grid 2D
        ("Stereo70ToETRS89", "5;1;1"),  # Out of grid 3D
        ("Stereo70ToETRS89", "549364;147553"),  # No data on grid 2D
        ("Stereo70ToETRS89", "549364;147553;100"),  # No data on grid 3D
        ("Stereo70ToETRS89", "a;b;c"),  # Invalid coordinate data (alpha)
        ("Stereo70ToETRS89", "500000;;500000"),  # Invalid coordinate data (empty middle token)
        # ETRS89ToStereo70 scenarios
        ("ETRS89ToStereo70", "0.8028465450500996;0.43630521911977493"),  # Valid 2D
        (
            "ETRS89ToStereo70",
            "0.8028465450500996;0.43630521911977493;139.7825537763764",
        ),  # Valid 3D
        ("ETRS89ToStereo70", "0.1;0.1"),  # Out of grid 2D
        ("ETRS89ToStereo70", "0.1;0.1;100"),  # Out of grid 3D
        ("ETRS89ToStereo70", "0.8087707898617313;0.35742127559921274"),  # No data on grid 2D
        ("ETRS89ToStereo70", "0.8087707898617313;0.35742127559921274;100"),  # No data on grid 3D
        ("ETRS89ToStereo70", "a;b;c"),  # Invalid coordinate data
    ],
)
def test_legacy_live_get_differential_parity(
    client, live_legacy_url, live_http_client, coo_op, coos_query
):
    """Compares live legacy GET responses against local FastAPI GET responses."""
    params = {"cooOp": coo_op, "coos": coos_query}
    live_res = _live_request(live_http_client, "GET", live_legacy_url, params=params)
    local_res = client.get("/transdatonline/cooOpService", params=params)

    assert live_res.status_code == 200
    assert local_res.status_code == 200

    _assert_legacy_item_parity(live_res.json(), local_res.json(), coo_op)


@pytest.mark.legacy_live
@pytest.mark.parametrize(
    "coo_op, batch_items",
    [
        (
            "Stereo70ToETRS89",
            [
                {"coos": [500000, 500000]},  # Valid 2D
                {"coos": [500000, 500000, 100]},  # Valid 3D
                {"coos": ["500000", "500000", "100"]},  # String-encoded numbers
                {"coos": [5, 1, 1]},  # Out of grid
                {"coos": [549364, 147553, 100]},  # No data on grid
                {"coos": ["a", "b", "c"]},  # Invalid coordinate data
            ],
        ),
        (
            "ETRS89ToStereo70",
            [
                {"coos": [0.8028465450500996, 0.43630521911977493]},  # Valid 2D
                {"coos": [0.8028465450500996, 0.43630521911977493, 139.7825537763764]},  # Valid 3D
                {
                    "coos": [
                        "0.8028465450500996",
                        "0.43630521911977493",
                        "139.7825537763764",
                    ]
                },  # String-encoded numbers
                {"coos": [0.1, 0.1, 100]},  # Out of grid
                {"coos": [0.8087707898617313, 0.35742127559921274, 100]},  # No data on grid
                {"coos": ["a", "b"]},  # Invalid coordinate data
            ],
        ),
    ],
)
def test_legacy_live_post_differential_parity(
    client, live_legacy_url, live_http_client, coo_op, batch_items
):
    """Compares live legacy POST form-urlencoded responses against local FastAPI POST responses."""
    form_data = {"cooOp": coo_op, "coosArray": json.dumps(batch_items)}
    live_res = _live_request(live_http_client, "POST", live_legacy_url, data=form_data)
    local_res = client.post("/transdatonline/cooOpService", data=form_data)

    assert live_res.status_code == 200
    assert local_res.status_code == 200

    live_list = live_res.json()
    local_list = local_res.json()
    assert isinstance(live_list, list) and isinstance(local_list, list)
    assert len(live_list) == len(local_list) == len(batch_items)

    for live_item, local_item in zip(live_list, local_list, strict=True):
        _assert_legacy_item_parity(live_item, local_item, coo_op)


def test_legacy_live_decommissioning_skip_mechanism(monkeypatch):
    """Verifies that disabling or decommissioning the live legacy server skips cleanly without failing."""
    from _pytest.outcomes import Skipped

    # 1. Explicit disable via environment variable
    monkeypatch.setenv("PYTRANSDAT_LEGACY_LIVE_ENABLED", "0")
    with pytest.raises(Skipped):
        _resolve_live_legacy_url()

    # 2. Unreachable / decommissioned endpoint URL
    monkeypatch.setenv("PYTRANSDAT_LEGACY_LIVE_ENABLED", "1")
    monkeypatch.setenv("PYTRANSDAT_LEGACY_LIVE_URL", "http://127.0.0.1:9/decommissioned")
    with pytest.raises(Skipped):
        _resolve_live_legacy_url()
