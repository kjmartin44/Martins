@echo off
REM Full analysis: every method enabled. Needs internet + rdkit + scikit-learn.
setlocal enabledelayedexpansion
cd /d "%~dp0"

echo === Checking dependencies ===
python -c "import numpy" 2>nul || (echo Installing numpy... & pip install numpy)
python -c "import rdkit, sklearn" 2>nul || (echo Installing rdkit + scikit-learn... & pip install rdkit scikit-learn)

REM --- motif file (required) -------------------------------------------
if not exist jaspar2024.txt (
    echo === Downloading JASPAR motifs ^(~2 MB^) ===
    curl.exe -L -o jaspar2024.txt "https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
)
for %%F in (jaspar2024.txt) do if %%~zF LSS 100000 (
    echo JASPAR download failed - delete jaspar2024.txt and retry.
    pause & exit /b 1
)

REM --- regulon file (optional, improves specificity scoring) ------------
if not exist trrust.tsv (
    echo === Downloading TRRUST regulons ^(optional^) ===
    curl.exe -L -o trrust.tsv "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
)
set REGULON=
for %%F in (trrust.tsv) do if %%~zF GTR 100000 set REGULON=--regulon-file trrust.tsv
if "%REGULON%"=="" echo   ^(TRRUST unavailable - using bundled tiers^)

REM --- inputs ----------------------------------------------------------
set "GENE="
set /p GENE="Gene symbol: "
if "%GENE%"=="" (echo No gene entered. & pause & exit /b 1)

set "DIR=up"
set /p DIR="Direction - up or down [up]: "
if "%DIR%"=="" set "DIR=up"

set "PROLIF=unknown"
set /p PROLIF="Tissue - cycling / post-mitotic / mixed / unknown [unknown]: "
if "%PROLIF%"=="" set "PROLIF=unknown"

set "TISS="
set /p CNS="Is the target tissue in the CNS? y/N: "
if /i "%CNS%"=="y" set "TISS=--tissue cns"

echo.
echo === Running every method for %GENE% (%DIR%regulation) ===
echo This takes 10-30 minutes on a first run. Results are cached afterwards.
echo.

python -m promoter_tf_drugs --gene %GENE% --direction %DIR% ^
    --pfm-file jaspar2024.txt %REGULON% ^
    --chip-api --upstream 5000 --downstream 1000 --all-tss ^
    --tf-weight 0.7 --proliferation %PROLIF% --cell-type-check %TISS% ^
    --predict-direction --upstream-network --upstream-depth 3 ^
    --lincs --genetics --screen --design-screen ^
    --out "results\%GENE%_%DIR%" -v

if errorlevel 1 (
    echo.
    echo Run failed - see the error above.
    echo.
    echo Common causes:
    echo   * Ensembl temporarily down ^(HTTP 500^) - wait a few minutes and rerun.
    echo     Everything already fetched is cached, so a rerun resumes quickly.
    echo   * Gene symbol not recognised - check it at https://www.genenames.org
    echo.
    pause & exit /b 1
)

echo.
echo === Done ===
echo Report:  results\%GENE%_%DIR%\report.md
echo Table:   results\%GENE%_%DIR%\tf_summary.tsv
echo Browser: results\%GENE%_%DIR%\sites.bed
pause
