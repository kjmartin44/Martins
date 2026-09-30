@echo off
REM Real analysis on a real gene. Needs an internet connection.
setlocal
cd /d "%~dp0.."

REM --- motif file: download once, reuse forever -------------------------
if not exist data\jaspar2024.txt (
    echo Downloading JASPAR motif file ^(~2 MB, one time^)...
    curl.exe -L -o data\jaspar2024.txt "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
    if errorlevel 1 (
        echo.
        echo Download failed. Check your internet connection, or download the file
        echo manually from https://jaspar.elixir.no/downloads/ and save it here as
        echo data\jaspar2024.txt
        pause
        exit /b 1
    )
)

REM --- sanity check: a real file is ~2 MB, an error page is tiny --------
for %%F in (data\jaspar2024.txt) do set SIZE=%%~zF
if %SIZE% LSS 100000 (
    echo.
    echo data\jaspar2024.txt is only %SIZE% bytes - the download did not work.
    echo Delete it and try again, or fetch it manually.
    pause
    exit /b 1
)

REM --- ask what to run --------------------------------------------------
set "GENE="
set /p GENE="Gene symbol (e.g. PAH, SMN2, LDLR): "
if "%GENE%"=="" (
    echo No gene entered.
    pause
    exit /b 1
)

set "DIR=up"
set /p DIR="Direction - type up or down [up]: "
if "%DIR%"=="" set "DIR=up"

echo.
echo Running %GENE% (%DIR%regulation). This hits Ensembl, ChEMBL and DGIdb,
echo so the first run takes a minute or two. Later runs are cached.
echo.

python -m promoter_tf_drugs --gene %GENE% --pfm-file data\jaspar2024.txt ^
    --direction %DIR% --out "results\%GENE%_%DIR%" -v

if errorlevel 1 (
    echo.
    echo The run failed - see the error above.
    pause
    exit /b 1
)

echo.
echo Done. Report written to results\%GENE%_%DIR%\report.md
pause
