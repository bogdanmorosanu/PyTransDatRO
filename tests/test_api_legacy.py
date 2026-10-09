"""Tests for 100% backward compatibility of the legacy TransDatOnline endpoint.
Validates /transdatonline/cooOpService and /cooOpService.
"""
import math
import pytest
from starlette.testclient import TestClient
from api.main import app


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


def test_legacy_get_root_alias_2d(client):
    """Verifies root alias /cooOpService works for 2D inputs."""
    res = client.get("/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000")
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
        "coosArray": '[{"coos":[500000, 500000, 100]}, {"coos":[5, 1, 1]}]'
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
        "coosArray": [
            {"coos": [0.8028465450500996, 0.43630521911977493, 139.7825537763764]}
        ]
    }
    res = client.post("/transdatonline/cooOpService", json=payload)
    assert res.status_code == 200
    data = res.json()

    assert isinstance(data, list)
    assert len(data) == 1
    assert math.isclose(data[0]["coos"][0], 500000.0, abs_tol=0.005)
