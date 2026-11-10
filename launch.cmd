@echo off
set PYTHONUTF8=1
cd /d "%~dp0"
if defined AMOEBA_PYTHON (
  "%AMOEBA_PYTHON%" gui.py
) else (
  python gui.py
)
if errorlevel 1 pause
