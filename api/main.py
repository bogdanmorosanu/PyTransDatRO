import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from starlette.exceptions import HTTPException as StarletteHTTPException

from api.config import configure_logging, get_settings
from api.middleware.rate_limit import SlidingWindowRateLimiter
from api.routes import batch, feedback, file, health, info, legacy, point
from pytransdatro.telemetry import SqliteUsageLogger
from pytransdatro.trans_ro import TransRO

logger = logging.getLogger("api.main")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan manager: pre-warms TransRO grid singleton and telemetry."""
    settings = get_settings()
    configure_logging(settings)
    logger.info("Initializing PyTransDatRO Web Service (Log Level: %s)", settings.log_level)

    # Initialize background SQLite WAL telemetry logger
    telemetry = SqliteUsageLogger(log_dir=settings.log_dir)
    app.state.telemetry = telemetry

    # Load SPG grid into RAM once on worker boot
    logger.info("Loading SPG grid into RAM...")
    trans = TransRO(telemetry_logger=telemetry)
    app.state.trans = trans
    logger.info("TransRO transformation engine pre-warmed and ready.")

    yield

    # Flush pending queue and gracefully close telemetry worker on exit
    logger.info("Shutting down PyTransDatRO Web Service worker.")
    if hasattr(app.state, "telemetry") and app.state.telemetry:
        app.state.telemetry.close()


def create_app() -> FastAPI:
    """Factory creating and configuring the FastAPI instance."""
    settings = get_settings()
    configure_logging(settings)

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

    # Global uniform exception handlers
    @app.exception_handler(StarletteHTTPException)
    async def http_exception_handler(request: Request, exc: StarletteHTTPException) -> JSONResponse:
        return JSONResponse(
            status_code=exc.status_code,
            content={
                "detail": exc.detail,
                "code": getattr(exc, "code", None),
            },
            headers=getattr(exc, "headers", None),
        )

    @app.exception_handler(RequestValidationError)
    async def validation_exception_handler(
        request: Request, exc: RequestValidationError
    ) -> JSONResponse:
        return JSONResponse(
            status_code=422,
            content={
                "detail": exc.errors(),
                "code": "VALIDATION_ERROR",
            },
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
    webapp_dist = settings.webapp_dist_dir
    if webapp_dist.is_dir():
        app.mount("/", StaticFiles(directory=str(webapp_dist), html=True), name="webapp")

    return app


app = create_app()
