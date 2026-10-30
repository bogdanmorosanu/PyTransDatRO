<#
.SYNOPSIS
    Starts the PyTransDatRO FastAPI backend and Vite frontend development servers.

.DESCRIPTION
    Launches Uvicorn (REST API) using the project conda environment and Vite (Web App),
    opens the application in the default web browser, and shuts down both child processes
    cleanly when Enter or Ctrl+C is pressed.
#>

[CmdletBinding()]
param(
    [string]$PythonPath = "",
    [int]$ApiPort = 8000,
    [switch]$NoBrowser
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
Write-Host "       PyTransDatRO Development Runner           " -ForegroundColor Cyan
Write-Host "=================================================" -ForegroundColor Cyan

# 1. Validate Python environment
Write-Host "Using Python: $pyExe" -ForegroundColor Green

# 2. Check for port conflicts on API port
$existingApi = Get-NetTCPConnection -LocalPort $ApiPort -ErrorAction SilentlyContinue | Select-Object -ExpandProperty OwningProcess -Unique
if ($existingApi) {
    Write-Warning "Port $ApiPort is currently occupied by process ID(s): $($existingApi -join ', ')."
    Write-Host "Please close the process or choose another port." -ForegroundColor Yellow
}

# 3. Start Backend (FastAPI / Uvicorn)
Write-Host "[1/2] Starting REST API on http://127.0.0.1:$ApiPort..." -ForegroundColor Green
$backendArgs = "-m uvicorn api.main:app --host 127.0.0.1 --port $ApiPort --reload"
$backendProc = Start-Process -FilePath $pyExe `
    -ArgumentList $backendArgs `
    -WorkingDirectory $rootDir `
    -PassThru

# 4. Start Frontend (Vite)
Write-Host "[2/2] Starting Web App (Vite)..." -ForegroundColor Green
$frontendProc = Start-Process -FilePath "cmd.exe" `
    -ArgumentList "/c npm run dev" `
    -WorkingDirectory $webappDir `
    -PassThru

# Brief pause to let Vite bind and display the port
Start-Sleep -Seconds 2

if (-not $NoBrowser) {
    Write-Host "Opening web app in default browser..." -ForegroundColor Magenta
    Start-Process "http://localhost:5173"
}

Write-Host ""
Write-Host "Development servers are active:" -ForegroundColor Cyan
Write-Host "  - REST API:     http://127.0.0.1:$ApiPort" -ForegroundColor White
Write-Host "  - Swagger Docs: http://127.0.0.1:$ApiPort/docs" -ForegroundColor White
Write-Host "  - Web App:      http://localhost:5173" -ForegroundColor White
Write-Host ""
Write-Host "Press [Enter] to stop all servers and exit..." -ForegroundColor Yellow

try {
    [void][System.Console]::ReadLine()
}
finally {
    Write-Host "`nStopping development servers..." -ForegroundColor Red
    
    if ($backendProc -and -not $backendProc.HasExited) {
        Stop-Process -Id $backendProc.Id -Force -ErrorAction SilentlyContinue
        # Kill any child python processes spawned by reload
        Get-CimInstance Win32_Process | Where-Object { $_.ParentProcessId -eq $backendProc.Id } | ForEach-Object {
            Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
        }
    }
    
    if ($frontendProc -and -not $frontendProc.HasExited) {
        # Kill cmd and node child processes
        Get-CimInstance Win32_Process | Where-Object { $_.ParentProcessId -eq $frontendProc.Id } | ForEach-Object {
            Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
        }
        Stop-Process -Id $frontendProc.Id -Force -ErrorAction SilentlyContinue
    }

    Write-Host "All servers stopped cleanly. Goodbye!" -ForegroundColor Green
}
