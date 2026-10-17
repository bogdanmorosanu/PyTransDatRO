"""Tests for sliding-window rate limiting middleware."""

from starlette.testclient import TestClient

from api.main import app


def test_rate_limiter_point_route():
    """Verifies that point requests exceeding the 60 req/min limit return HTTP 429."""
    app.state.disable_rate_limiting = False
    with TestClient(app) as client:
        # Send 60 requests (allowed)
        for i in range(60):
            res = client.get(
                "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000"
            )
            assert res.status_code == 200, f"Request {i + 1} failed unexpectedly"

        # 61st request must trigger 429 Too Many Requests
        res_blocked = client.get(
            "/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000"
        )
        assert res_blocked.status_code == 429
        assert "rate limit exceeded" in res_blocked.json()["detail"].lower()
        assert "retry-after" in res_blocked.headers


def test_rate_limiter_batch_route():
    """Verifies that batch requests exceeding the 10 req/min limit return HTTP 429."""
    app.state.disable_rate_limiting = False
    with TestClient(app) as client:
        payload = {
            "op": "Stereo70ToETRS89",
            "points": [[500000.0, 500000.0, 100.0]],
            "unit": "radians",
        }
        # Send 10 requests (allowed)
        for i in range(10):
            res = client.post("/api/v1/transform/batch", json=payload)
            assert res.status_code == 200, f"Batch request {i + 1} failed"

        # 11th request must trigger 429
        res_blocked = client.post("/api/v1/transform/batch", json=payload)
        assert res_blocked.status_code == 429
        assert "Batch rate limit exceeded" in res_blocked.json()["detail"]
        assert "retry-after" in res_blocked.headers


def test_rate_limiter_modern_point_per_ip_isolation_and_unmetered_routes():
    """Verifies /api/v1/transform/point rate limit, X-Forwarded-For IP isolation, and unmetered routes."""
    app.state.disable_rate_limiting = False
    ip_a = {"X-Forwarded-For": "198.51.100.10"}
    ip_b = {"X-Forwarded-For": "198.51.100.20"}
    payload = {"op": "Stereo70ToETRS89", "coos": [500000.0, 500000.0], "unit": "degrees"}

    with TestClient(app) as client:
        for _ in range(60):
            res = client.post("/api/v1/transform/point", json=payload, headers=ip_a)
            assert res.status_code == 200

        # 61st request from IP A is rate-limited
        res_blocked = client.post("/api/v1/transform/point", json=payload, headers=ip_a)
        assert res_blocked.status_code == 429

        # IP B is completely unaffected
        res_ip_b = client.post("/api/v1/transform/point", json=payload, headers=ip_b)
        assert res_ip_b.status_code == 200

        # Unmetered routes (/health, /api/v1/grid/info) remain accessible for IP A
        assert client.get("/health", headers=ip_a).status_code == 200
        assert client.get("/api/v1/grid/info", headers=ip_a).status_code == 200


def test_rate_limiter_legacy_post_and_file_upload_share_batch_bucket():
    """Verifies POST /transdatonline/cooOpService and POST /api/v1/transform/file share the 10 req/min bucket."""
    import io

    app.state.disable_rate_limiting = False
    headers = {"X-Forwarded-For": "198.51.100.30"}
    with TestClient(app) as client:
        # 5 Legacy POST requests
        for _ in range(5):
            res = client.post(
                "/transdatonline/cooOpService",
                json={"cooOp": "Stereo70ToETRS89", "coosArray": [{"coos": [500000.0, 500000.0]}]},
                headers=headers,
            )
            assert res.status_code == 200

        # 5 Modern File Upload requests (total = 10 batch requests)
        for _ in range(5):
            res = client.post(
                "/api/v1/transform/file",
                files={"file": ("t.csv", io.BytesIO(b"500000.0,500000.0\n"), "text/csv")},
                data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"},
                headers=headers,
            )
            assert res.status_code == 200

        # 11th request on either endpoint triggers 429
        res_legacy_blocked = client.post(
            "/transdatonline/cooOpService",
            json={"cooOp": "Stereo70ToETRS89", "coosArray": [{"coos": [500000.0, 500000.0]}]},
            headers=headers,
        )
        assert res_legacy_blocked.status_code == 429

        res_file_blocked = client.post(
            "/api/v1/transform/file",
            files={"file": ("t.csv", io.BytesIO(b"500000.0,500000.0\n"), "text/csv")},
            data={"op": "Stereo70ToETRS89", "unit": "degrees", "delimiter": "comma"},
            headers=headers,
        )
        assert res_file_blocked.status_code == 429


def test_rate_limiter_sliding_window_expiration(monkeypatch):
    """Verifies that advancing time beyond the 60s window evicts old hits and restores access."""
    app.state.disable_rate_limiting = False
    headers = {"X-Forwarded-For": "198.51.100.40"}
    payload = {"op": "Stereo70ToETRS89", "points": [[500000.0, 500000.0]], "unit": "degrees"}

    current_time = 1_700_000_000.0
    monkeypatch.setattr("api.middleware.rate_limit.time.time", lambda: current_time)

    with TestClient(app) as client:
        for _ in range(10):
            assert (
                client.post("/api/v1/transform/batch", json=payload, headers=headers).status_code
                == 200
            )
        assert (
            client.post("/api/v1/transform/batch", json=payload, headers=headers).status_code == 429
        )

        # Advance clock by 61 seconds past the 60-second window
        current_time += 61.0
        res_recovered = client.post("/api/v1/transform/batch", json=payload, headers=headers)
        assert res_recovered.status_code == 200
