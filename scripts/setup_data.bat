@echo off
REM Fetch every reference file the analysis reads from data\ and that git
REM does not carry: motifs, ChIP peaks, annotation, expression, LINCS,
REM DepMap, PRISM, GWAS, gnomAD and the rest. See scripts\setup_data.py.
REM
REM Some tables are BUILT rather than downloaded - the named transcript
REM regulators are assembled from TargetScan and ENCODE eCLIP by
REM scripts\build_transcript_sites.py, which this script calls for you.
REM
REM Two are neither fetched nor built, and the run says so rather than
REM reporting the question as answered:
REM   data\drug_protein_abundance.tsv   decryptE or ProTargetMiner, placed
REM                                     by hand - see docs\methods.md
REM   data\transcript_expression.tsv    optional; without it GTEx is
REM                                     queried live instead
REM
REM Usage:  setup_data.bat              (shows what is missing, asks, fetches)
REM         setup_data.bat --yes        (no question)
REM         setup_data.bat --list       (inventory only)
REM         setup_data.bat --small      (skip the files over about 500 MB)
REM         setup_data.bat --with-optional
REM         setup_data.bat --only remap,collectri
REM
REM Files already in data\ are left alone, so this is safe to rerun after a
REM failure or to add a file you deleted.
setlocal
cd /d "%~dp0.."

set "PY=python"
if exist ".venv\Scripts\python.exe" set "PY=.venv\Scripts\python.exe"
if exist "venv\Scripts\python.exe" set "PY=venv\Scripts\python.exe"
echo Interpreter: %PY%

"%PY%" -c "import numpy" 2>nul || ("%PY%" -m pip install numpy)
if not exist data mkdir data

"%PY%" scripts\setup_data.py %*
set "RC=%errorlevel%"
if not "%RC%"=="0" (
    echo.
    echo Some files were not fetched. Rerun this script to retry only those.
)
if not defined SETUP_NOPAUSE pause
exit /b %RC%
