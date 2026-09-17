@echo off
REM PrimerForge network server - logins on, reachable from other computers on port 8080.
REM Don't run this alongside run.bat: both would start the same BLAT servers and job queue.
cd /d "%~dp0"
set PRIMERFORGE_AUTH=1
set PRIMERFORGE_HOST=0.0.0.0
if not defined PRIMERFORGE_PORT set PRIMERFORGE_PORT=8080
python serve.py >> data\server.log 2>&1
