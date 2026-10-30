<#
.SYNOPSIS
    Runs the Playwright End-to-End (E2E) UI test suite for PyTransDatRO.

.DESCRIPTION
    Builds the web application into dist-e2e, boots an isolated FastAPI test server
    using the conda python environment, and executes the Playwright test runner.
    Supports interactive UI mode with -Ui or -Interactive.

.EXAMPLE
    .\test-e2e.ps1
    Runs all 19 E2E browser tests in headless mode.

.EXAMPLE
    .\test-e2e.ps1 -Ui
    Opens Playwright's visual interactive UI runner with time-travel debugging.
#>

[CmdletBinding()]
param(
    [string]$PythonPath = "",
    [switch]$Ui,
    [switch]$Interactive
)

$rootDir = $PSScriptRoot
$webappDir = Join-Path $rootDir "webapp"

function Find-ProjectPython {
    param([string]$Explicit)
    if ($Explicit -and (Test-Path $Explicit)) { return $Explicit }
    if ($env:PYTHON -and (Test-Path $env:PYTHON)) { return $env:PYTHON }
    if ($env:CONDA_PREFIX) {
        $condaPy = Join-Path $env:CONDA_PREFIX "python.exe"
        if (Test-Path $condaPy) { return $condaPy }
    }
    $candidates = @(
        (Join-Path $env:USERPROFILE "anaconda3\envs\pytransdat\python.exe"),
        (Join-Path $env:USERPROFILE "miniconda3\envs\pytransdat\python.exe"),
        "D:\anaconda3\envs\pytransdat\python.exe",
        (Join-Path $rootDir ".venv\Scripts\python.exe")
    )
    foreach ($cand in $candidates) {
        if ($cand -and (Test-Path $cand)) { return $cand }
    }
    $cmd = Get-Command python -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    return "python"
}

$pyExe = Find-ProjectPython -Explicit $PythonPath
$env:PYTHON = $pyExe

Write-Host "=================================================" -ForegroundColor Cyan
Write-Host "      PyTransDatRO E2E Test Runner (Playwright)  " -ForegroundColor Cyan
Write-Host "=================================================" -ForegroundColor Cyan

# 1. Validate Python environment
Write-Host "Using Python: $pyExe" -ForegroundColor Green

# 2. Execute Playwright tests
Push-Location $webappDir
try {
    if ($Ui -or $Interactive) {
        Write-Host "Launching Playwright Interactive UI Mode..." -ForegroundColor Yellow
        npm run test:e2e:ui
    } else {
        Write-Host "Running Playwright E2E Suite (Headless)..." -ForegroundColor Green
        npm run test:e2e
    }
} finally {
    Pop-Location
}
