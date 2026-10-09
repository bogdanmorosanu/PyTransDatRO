"""Tests for sliding-window rate limiting middleware.
"""
from starlette.testclient import TestClient
from api.main import app


def test_rate_limiter_point_route():
    """Verifies that point requests exceeding the 60 req/min limit return HTTP 429."""
    app.state.disable_rate_limiting = False
    with TestClient(app) as client:
        # Send 60 requests (allowed)
        for i in range(60):
            res = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000")
            assert res.status_code == 200, f"Request {i+1} failed unexpectedly"

        # 61st request must trigger 429 Too Many Requests
        res_blocked = client.get("/transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000")
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
            "unit": "radians"
        }
        # Send 10 requests (allowed)
        for i in range(10):
            res = client.post("/api/v1/transform/batch", json=payload)
            assert res.status_code == 200, f"Batch request {i+1} failed"

        # 11th request must trigger 429
        res_blocked = client.post("/api/v1/transform/batch", json=payload)
        assert res_blocked.status_code == 429
        assert "Batch rate limit exceeded" in res_blocked.json()["detail"]
        assert "retry-after" in res_blocked.headers
