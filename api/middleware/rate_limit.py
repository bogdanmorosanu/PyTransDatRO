"""In-memory sliding-window rate limiting middleware for FastAPI.
Zero external dependencies, standard library only.
"""
import time
from collections import defaultdict, deque
from typing import Dict, Deque
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse, Response


class SlidingWindowRateLimiter(BaseHTTPMiddleware):
    """Sliding-window in-memory rate limiter.
    
    Default rules:
    - Point transformations (GET legacy, POST /point): 60 requests / minute / IP
    - Batch and file uploads (POST legacy, POST /batch, POST /file): 10 requests / minute / IP
    - Health, grid info, documentation: unmetered
    """
    def __init__(
        self,
        app,
        point_limit: int = 60,
        batch_limit: int = 10,
        window_seconds: float = 60.0
    ):
        super().__init__(app)
        self.point_limit = point_limit
        self.batch_limit = batch_limit
        self.window_seconds = window_seconds

        # IP -> deque of request timestamps
        self._point_hits: Dict[str, Deque[float]] = defaultdict(deque)
        self._batch_hits: Dict[str, Deque[float]] = defaultdict(deque)

    def _get_client_ip(self, request: Request) -> str:
        forwarded = request.headers.get("x-forwarded-for")
        if forwarded:
            return forwarded.split(",")[0].strip()
        if request.client and request.client.host:
            return request.client.host
        return "127.0.0.1"

    def _check_rate_limit(self, hits_map: Dict[str, Deque[float]], client_ip: str, limit: int) -> int:
        """Checks and records request. Returns 0 if allowed, or retry_after in seconds if rate-limited."""
        now = time.time()
        cutoff = now - self.window_seconds
        client_hits = hits_map[client_ip]

        # Evict timestamps outside the window
        while client_hits and client_hits[0] < cutoff:
            client_hits.popleft()

        if len(client_hits) >= limit:
            oldest = client_hits[0]
            retry_after = max(1, int(self.window_seconds - (now - oldest)))
            return retry_after

        client_hits.append(now)
        return 0

    async def dispatch(self, request: Request, call_next) -> Response:
        # Allow disabling rate limiting during specific test suites if flag set
        if getattr(request.app.state, "disable_rate_limiting", False):
            return await call_next(request)

        path = request.url.path.lower()
        method = request.method.upper()

        # Classify route
        is_point_route = (
            ("/api/v1/transform/point" in path) or
            ("cooopservice" in path and method == "GET")
        )
        is_batch_route = (
            ("/api/v1/transform/batch" in path) or
            ("/api/v1/transform/file" in path) or
            ("cooopservice" in path and method == "POST")
        )

        client_ip = self._get_client_ip(request)

        if is_batch_route:
            retry_after = self._check_rate_limit(self._batch_hits, client_ip, self.batch_limit)
            if retry_after > 0:
                return JSONResponse(
                    status_code=429,
                    content={"detail": "Too Many Requests. Batch rate limit exceeded (10 requests/minute)."},
                    headers={"Retry-After": str(retry_after)}
                )

        elif is_point_route:
            retry_after = self._check_rate_limit(self._point_hits, client_ip, self.point_limit)
            if retry_after > 0:
                return JSONResponse(
                    status_code=429,
                    content={"detail": "Too Many Requests. Point rate limit exceeded (60 requests/minute)."},
                    headers={"Retry-After": str(retry_after)}
                )

        return await call_next(request)
