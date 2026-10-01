@echo off
REM Double-click this to check the install works. No internet needed.
cd /d "%~dp0.."
python -m promoter_tf_drugs --demo --gene DEMOGENE --out results\demo
echo.
echo Done. Open results\demo\report.md
pause
