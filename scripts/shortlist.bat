@echo off
REM Filter an existing results.json down to a shortlist. No network,
REM no re-run, instant.
setlocal enabledelayedexpansion
cd /d "%~dp0.."

set "PY=python"
if exist ".venv\Scripts\python.exe" set "PY=.venv\Scripts\python.exe"
if exist "venv\Scripts\python.exe" set "PY=venv\Scripts\python.exe"

set "RESULTS="
set /p RESULTS="Path to results.json (e.g. results\PAH_up\results.json): "
if "!RESULTS!"=="" (
    echo No path entered.
    pause
    exit /b 1
)
if not exist "!RESULTS!" (
    echo File not found: !RESULTS!
    pause
    exit /b 1
)

echo.
echo   strict     what you would actually order next
echo   mechanism  routes you understand well enough to design around
echo   empirical  only what has been measured
echo   specific   nothing that moves half the transcriptome
echo   novel      what nobody has tried on this gene
echo.
set "PRESET=strict"
set /p PRESET="Preset [strict]: "
if "!PRESET!"=="" set "PRESET=strict"

REM Checked here rather than left to argparse, because the preset also
REM becomes part of the output filename below - an unchecked value ends
REM up in a path.
set "OK="
for %%P in (strict mechanism empirical specific novel) do (
    if /i "!PRESET!"=="%%P" set "OK=1"
)
if not defined OK (
    echo "!PRESET!" is not one of the presets listed above.
    pause
    exit /b 1
)

"%PY%" -m generoutes.shortlist "!RESULTS!" --preset !PRESET! ^
    --out "shortlist_!PRESET!.md"
if errorlevel 1 (
    echo.
    echo The shortlist failed - see the error above.
    pause
    exit /b 1
)

echo.
echo Written to shortlist_!PRESET!.md
pause
