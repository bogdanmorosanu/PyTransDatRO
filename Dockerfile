# ==============================================================================
# Stage 1: Build Frontend Web Application (TypeScript + Vite)
# ==============================================================================
FROM node:20-alpine AS frontend-builder
WORKDIR /build

# Install dependencies with lockfile consistency
COPY webapp/package.json webapp/package-lock.json ./
RUN npm ci

# Copy web application source and build production bundle
COPY webapp/ ./
RUN npm run build

# ==============================================================================
# Stage 2: Python Backend & Unified Web Service Runtime
# ==============================================================================
FROM python:3.11-slim AS runtime
WORKDIR /app

# Ensure Python doesn't buffer standard output and doesn't write .pyc files
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

# Install Python dependencies
COPY requirements.txt requirements-api.txt ./
RUN pip install --no-cache-dir -r requirements.txt -r requirements-api.txt

# Copy backend service and core transformation library (including binary grids)
COPY pytransdatro/ /app/pytransdatro/
COPY api/ /app/api/

# Copy compiled webapp static distribution from Stage 1 into webapp/dist
# (FastAPI mounts /webapp/dist automatically at "/" when present)
COPY --from=frontend-builder /build/dist /app/webapp/dist

# Ensure persistent telemetry directory exists
RUN mkdir -p /app/logs

# Expose FastAPI port
EXPOSE 8000

# Container healthcheck using standard library urllib
HEALTHCHECK --interval=30s --timeout=5s --start-period=5s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000/health')" || exit 1

# Launch production Uvicorn server
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
