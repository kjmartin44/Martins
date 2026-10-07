@echo off
REM The ReMap 2022 peak file: measured TF occupancy, and the single
REM largest gain available to a run. About 1.4 GB.
REM
REM Checked three ways, because a download that fails quietly is worse
REM than one that fails loudly: curl is told to fail on an HTTP error
REM rather than save the error page, its exit code is read, and the
REM result is opened as gzip to confirm it is the file rather than an
REM apology page that happens to be large.
setlocal
cd /d "%~dp0.."

set "PY=python"
if exist ".venv\Scripts\python.exe" set "PY=.venv\Scripts\python.exe"
if exist "venv\Scripts\python.exe" set "PY=venv\Scripts\python.exe"

set "URL=https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz"
set "DEST=data\remap2022.bed.gz"
set "TMP=data\remap2022.bed.gz.part"

if not exist data mkdir data
if exist "%DEST%" (
    echo %DEST% already exists. Delete it to download again.
    pause
    exit /b 0
)

echo Downloading ReMap 2022, about 1.4 GB. This takes a while.
echo.
REM -f makes curl exit non-zero on 4xx/5xx rather than saving the body.
REM Downloaded to a .part file so a failure cannot leave something that
REM looks like the real thing in its place.
curl.exe -fL --retry 3 -o "%TMP%" "%URL%"
if errorlevel 1 (
    echo.
    echo Download failed. Nothing was written to %DEST%.
    echo Check it in a browser: https://remap.univ-amu.fr/
    if exist "%TMP%" del "%TMP%"
    pause
    exit /b 1
)

REM Content, not size. A captive portal or an error page can be any
REM length at all; gzip is the one thing an HTML page is not.
"%PY%" -c "import gzip,sys; gzip.open(sys.argv[1],'rb').read(1024)" "%TMP%"
if errorlevel 1 (
    echo.
    echo What arrived is not a gzip file - most likely an error page.
    echo Left at %TMP% so you can look at it.
    pause
    exit /b 1
)

move /y "%TMP%" "%DEST%" >nul
echo.
echo Done: %DEST%
echo Use it with  --chip-bed %DEST%
pause
