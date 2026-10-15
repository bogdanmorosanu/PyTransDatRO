"""FastAPI application factory and entry point for PyTransDatRO Web Service."""

from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from api.config import get_settings
from api.middleware.rate_limit import SlidingWindowRateLimiter
from api.routes import batch, feedback, file, health, info, legacy, point
from pytransdatro.telemetry import SqliteUsageLogger
from pytransdatro.trans_ro import TransRO


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan manager: pre-warms TransRO grid singleton and telemetry."""
    settings = get_settings()
    # Initialize background SQLite WAL telemetry logger
    telemetry = SqliteUsageLogger(log_dir=settings.log_dir)
    app.state.telemetry = telemetry

    # Load SPG grid into RAM once on worker boot
    trans = TransRO(telemetry_logger=telemetry)
    app.state.trans = trans

    yield

    # Flush pending queue and gracefully close telemetry worker on exit
    if hasattr(app.state, "telemetry") and app.state.telemetry:
        app.state.telemetry.close()


def create_app() -> FastAPI:
    """Factory creating and configuring the FastAPI instance."""
    settings = get_settings()

    app = FastAPI(
        title="PyTransDatRO Web API",
        description=(
            "High-performance REST API for Romanian geodetic coordinate transformations "
            "(Stereo 70 <-> ETRS89), featuring 100% backward compatibility with the legacy "
            "TransDatOnline service, modern streaming endpoints, and zero-overhead telemetry."
        ),
        version="1.0.0",
        lifespan=lifespan,
        docs_url="/docs",
        redoc_url="/redoc",
    )

    # Enable configurable CORS for web map applications and GIS plugins
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.parsed_cors_origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # In-memory sliding-window rate limiting middleware configured from settings
    app.add_middleware(
        SlidingWindowRateLimiter,
        point_limit=settings.rate_limit_point,
        batch_limit=settings.rate_limit_batch,
        window_seconds=settings.rate_limit_window_seconds,
    )

    # Register routers
    app.include_router(legacy.router)
    app.include_router(point.router)
    app.include_router(batch.router)
    app.include_router(file.router)
    app.include_router(info.router)
    app.include_router(health.router)
    app.include_router(feedback.router)

    # Mount static web application assets if built
    webapp_dist = Path(__file__).resolve().parent.parent / "webapp" / "dist"
    if webapp_dist.is_dir():
        app.mount("/", StaticFiles(directory=str(webapp_dist), html=True), name="webapp")

    return app


app = create_app()
