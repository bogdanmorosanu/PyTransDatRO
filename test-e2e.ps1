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
    [string]$PythonPath = "D:\anaconda3\envs\pytransdat\python.exe",
    [switch]$Ui,
    [switch]$Interactive
)

$rootDir = $PSScriptRoot
$webappDir = Join-Path $rootDir "webapp"

Write-Host "=================================================" -ForegroundColor Cyan
Write-Host "      PyTransDatRO E2E Test Runner (Playwright)  " -ForegroundColor Cyan
Write-Host "=================================================" -ForegroundColor Cyan

# 1. Validate Python environment
if (Test-Path $PythonPath) {
    $env:PYTHON = $PythonPath
    Write-Host "Using Python: $PythonPath" -ForegroundColor Green
} else {
    Write-Warning "Conda python not found at $PythonPath, falling back to system 'python'..."
}

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
