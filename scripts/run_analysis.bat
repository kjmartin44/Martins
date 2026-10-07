@echo off
REM Real analysis on a real gene. Needs an internet connection.
REM For every method at once, including docking, use run_everything.bat.
setlocal enabledelayedexpansion
cd /d "%~dp0.."

REM A virtual environment in the checkout wins, so the script does not
REM run against whatever "python" happens to be first on PATH.
set "PY=python"
if exist ".venv\Scripts\python.exe" set "PY=.venv\Scripts\python.exe"
if exist "venv\Scripts\python.exe" set "PY=venv\Scripts\python.exe"

REM --- motif file: download once, reuse forever -------------------------
if not exist data\jaspar2024.txt (
    echo Downloading JASPAR motif file ^(~2 MB, one time^)...
    REM -f so an HTTP error exits non-zero rather than saving the error
    REM page under the name of the file you wanted.
    curl.exe -fL --retry 3 -o data\jaspar2024.txt "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    if errorlevel 1 (
        echo.
        echo Download failed. Check your connection, or fetch it manually
        echo from https://jaspar.elixir.no/downloads/ and save it here as
        echo data\jaspar2024.txt
        if exist data\jaspar2024.txt del data\jaspar2024.txt
        pause
        exit /b 1
    )
)

REM --- is it the file, or something else of about the right size? ------
REM Checked by CONTENT rather than by length. A captive portal page or a
REM maintenance notice can be any size at all, and the old check - "at
REM least 100 KB" - passes anything long enough to look plausible. Every
REM JASPAR matrix header begins with ">MA".
findstr /m /c:">MA" data\jaspar2024.txt >nul 2>&1
if errorlevel 1 (
    echo.
    echo data\jaspar2024.txt is not a JASPAR PFM file - the download
    echo returned something else. Delete it and try again.
    pause
    exit /b 1
)

REM --- ask what to run --------------------------------------------------
set "GENE="
set /p GENE="Gene symbol (e.g. PAH, SMN2, LDLR): "
if "!GENE!"=="" (
    echo No gene entered.
    pause
    exit /b 1
)
REM Letters, digits, hyphen and underscore. Anything else is a typo, or
REM an attempt to put a second command on the line - and the value is
REM used in an output path below.
echo !GENE!| findstr /r /c:"^[A-Za-z0-9_-][A-Za-z0-9_-]*$" >nul
if errorlevel 1 (
    echo "!GENE!" is not a gene symbol.
    pause
    exit /b 1
)

set "DIR=up"
set /p DIR="Direction - type up or down [up]: "
if "!DIR!"=="" set "DIR=up"
if /i not "!DIR!"=="up" if /i not "!DIR!"=="down" (
    echo Direction must be up or down, not "!DIR!".
    pause
    exit /b 1
)

echo.
echo Running !GENE! ^(!DIR!regulation^). This hits Ensembl, ChEMBL and
echo DGIdb, so a first run takes a minute or two. Later runs are cached.
echo Checkpoints print as it goes.
echo.

"%PY%" -m promoter_tf_drugs --gene !GENE! --pfm-file data\jaspar2024.txt ^
    --direction !DIR! --out "results\!GENE!_!DIR!" -v

if errorlevel 1 (
    echo.
    echo The run failed - see the error above.
    pause
    exit /b 1
)

echo.
echo Done. Report written to results\!GENE!_!DIR!\report.md
echo.
echo Read the evidence table at the top before the ranking underneath
echo it: a source that failed is recorded as failed, never as zero.
pause
