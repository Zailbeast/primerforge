@echo off
REM PrimerForge launcher - double-click to start the app.
cd /d "%~dp0"
echo Starting PrimerForge...
python -m pip install -q -r requirements.txt 2>nul
python app.py %*
if errorlevel 1 (
  echo.
  echo PrimerForge exited with an error. Press any key to close.
  pause >nul
)
