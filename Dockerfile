# ==============================================================================
# PyTransDatRO - Backend REST API & Core Geodetic Engine Runtime
# ==============================================================================
FROM python:3.11-slim AS runtime
WORKDIR /app

# Prevent Python from writing .pyc files and buffer stdout/stderr
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

# Install Python runtime dependencies
COPY requirements.txt requirements-api.txt ./
RUN pip install --no-cache-dir -r requirements.txt -r requirements-api.txt

# Copy backend service and pure Python core transformation library (including binary grids)
COPY pytransdatro/ /app/pytransdatro/
COPY api/ /app/api/

# Ensure persistent telemetry and logging directory exists
RUN mkdir -p /app/logs

# Expose internal API port
EXPOSE 8000

# Container healthcheck: verifies grid pre-warming and telemetry initialization
HEALTHCHECK --interval=15s --timeout=5s --start-period=5s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')" || exit 1

# Launch production Uvicorn ASGI server
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
