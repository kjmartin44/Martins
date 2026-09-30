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

REM --- scan window and threshold ---------------------------------------
REM A 6 kb window at p=1e-4 across ~880 matrices expects ~740 chance hits.
REM Real promoters do not contain 900 genuine sites, so that combination
REM returns a list that is ~83%% noise no matter which gene you run - the
REM ranking below it is then a ranking of chance matches. A core-promoter
REM window at a stricter threshold is what makes the scan interpretable.
REM Widen deliberately (and tighten --pvalue further if you do), not by
REM default.
set "WINDOW=--upstream 1000 --downstream 200 --pvalue 1e-5"

REM --- ChIP evidence ---------------------------------------------------
REM This is the step that turns "could bind" into "was observed bound", and
REM it is worth more than every scoring tweak combined. A local ReMap BED is
REM parsed directly; the API is used only as a fallback because a changed
REM response shape silently yields no peaks at all.
set "CHIP=--chip-api"
if exist remap2022.bed.gz set "CHIP=--chip-bed remap2022.bed.gz"
if "%CHIP%"=="--chip-api" echo   ^(no local remap2022.bed.gz - falling back to the ReMap API^)

REM --- tissue expression ------------------------------------------------
REM Open Targets' expression endpoint changes shape periodically; when it
REM does, every factor reports "unknown" and the co-expression filter
REM silently does nothing. A local HPA/GTEx table removes the dependency.
set "EXPR="
if exist rna_tissue_consensus.tsv set "EXPR=--expression-file rna_tissue_consensus.tsv"
if "%EXPR%"=="" echo   ^(no local rna_tissue_consensus.tsv - using the Open Targets API^)

REM --tf-weight: 0.7 puts most of each grade on the promoter call, which is
REM only justified once ChIP evidence is actually present. Without it the
REM weight sits back at 0.5 so a motif guess cannot carry a grade.
set "TFW=--tf-weight 0.5"
if not "%CHIP%"=="--chip-api" set "TFW=--tf-weight 0.7"

python -m promoter_tf_drugs --gene %GENE% --direction %DIR% ^
    --pfm-file jaspar2024.txt %REGULON% ^
    %CHIP% %WINDOW% --all-tss ^
    %TFW% --proliferation %PROLIF% --cell-type-check %TISS% %EXPR% ^
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
