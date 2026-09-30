@echo off
REM Filter an existing results.json down to a shortlist. No network, instant.
setlocal
cd /d "%~dp0"

set "RESULTS="
set /p RESULTS="Path to results.json (e.g. results\PAH_up\results.json): "
if not exist "%RESULTS%" (echo File not found. & pause & exit /b 1)

echo.
echo   strict     what you would actually order next
echo   mechanism  routes you understand well enough to design around
echo   empirical  only what has been measured
echo   specific   nothing that moves half the transcriptome
echo   novel      what nobody has tried on this gene
echo.
set "PRESET=strict"
set /p PRESET="Preset [strict]: "
if "%PRESET%"=="" set "PRESET=strict"

python -m promoter_tf_drugs.shortlist "%RESULTS%" --preset %PRESET% --out shortlist_%PRESET%.md
echo.
echo Written to shortlist_%PRESET%.md
pause
