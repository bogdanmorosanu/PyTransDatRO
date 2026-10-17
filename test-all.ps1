<#
.SYNOPSIS
    Runs the complete PyTransDatRO verification suite (Core, REST API, Web Unit, E2E).

.DESCRIPTION
    Executes in sequence:
      1. Python unit, integration, and API tests (pytest)
      2. Frontend unit and parsing tests (vitest)
      3. Browser End-to-End flows (playwright)
    Displays a concise consolidated status report at the end.

.EXAMPLE
    .\test-all.ps1
#>

[CmdletBinding()]
param(
    [string]$PythonPath = "D:\anaconda3\envs\pytransdat\python.exe",
    [switch]$SkipE2E
)

$rootDir = $PSScriptRoot
$webappDir = Join-Path $rootDir "webapp"

Write-Host "=================================================" -ForegroundColor Cyan
Write-Host "       PyTransDatRO Full Verification Suite      " -ForegroundColor Cyan
Write-Host "=================================================" -ForegroundColor Cyan

$pyExe = if (Test-Path $PythonPath) { $PythonPath } else { "python" }
$env:PYTHON = $pyExe

$results = [ordered]@{}

# 1. Python pytest suite
Write-Host "`n[1/3] Running Python Core & REST API Suite (pytest)..." -ForegroundColor Yellow
& $pyExe -m pytest -m "not legacy_live" -q
if ($LASTEXITCODE -eq 0) {
    $results["Python (Core & API)"] = "PASSED"
} else {
    $results["Python (Core & API)"] = "FAILED"
}

# 2. Frontend Vitest unit suite
Write-Host "`n[2/3] Running Frontend Unit & Parser Suite (vitest)..." -ForegroundColor Yellow
npm --prefix $webappDir test
if ($LASTEXITCODE -eq 0) {
    $results["Frontend Unit (Vitest)"] = "PASSED"
} else {
    $results["Frontend Unit (Vitest)"] = "FAILED"
}

# 3. Playwright E2E browser suite
if (-not $SkipE2E) {
    Write-Host "`n[3/3] Running Browser End-to-End Suite (playwright)..." -ForegroundColor Yellow
    npm --prefix $webappDir run test:e2e
    if ($LASTEXITCODE -eq 0) {
        $results["Frontend E2E (Playwright)"] = "PASSED"
    } else {
        $results["Frontend E2E (Playwright)"] = "FAILED"
    }
} else {
    $results["Frontend E2E (Playwright)"] = "SKIPPED"
}

# Consolidated summary
Write-Host "`n=================================================" -ForegroundColor Cyan
Write-Host "                 SUMMARY REPORT                  " -ForegroundColor Cyan
Write-Host "=================================================" -ForegroundColor Cyan

$allPassed = $true
foreach ($key in $results.Keys) {
    $status = $results[$key]
    if ($status -eq "PASSED") {
        Write-Host "  [✓] $key : $status" -ForegroundColor Green
    } elseif ($status -eq "SKIPPED") {
        Write-Host "  [-] $key : $status" -ForegroundColor Yellow
    } else {
        Write-Host "  [✕] $key : $status" -ForegroundColor Red
        $allPassed = $false
    }
}
Write-Host "=================================================" -ForegroundColor Cyan

if (-not $allPassed) {
    exit 1
}
