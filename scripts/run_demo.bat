@echo off
REM Double-click this to check the install works. No internet needed:
REM the demo runs on synthetic fixtures and makes no request at all.
setlocal
cd /d "%~dp0.."

REM A virtual environment in the checkout wins. Without this the script
REM runs against whatever "python" is first on PATH, which on a machine
REM with several installs is rarely the one the package was installed
REM into - and the failure then reads as a broken package rather than
REM the wrong interpreter.
set "PY=python"
if exist ".venv\Scripts\python.exe" set "PY=.venv\Scripts\python.exe"
if exist "venv\Scripts\python.exe" set "PY=venv\Scripts\python.exe"

"%PY%" -c "import generoutes" 2>nul
if errorlevel 1 (
    echo The package is not importable by %PY%.
    echo.
    echo From this folder:
    echo     python -m pip install -e .
    pause
    exit /b 1
)

"%PY%" -m generoutes --demo --gene DEMOGENE --out results\demo
if errorlevel 1 (
    echo.
    echo The demo failed. That is an installation problem rather than a
    echo data one: it reads no files and makes no requests.
    pause
    exit /b 1
)

echo.
echo Done. Open results\demo\report.md
pause
