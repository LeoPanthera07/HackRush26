@echo off
cd /d "%~dp0"
echo Running VMD psfgen for 1GGY_A_732_YB...
vmd -dispdev text -e psfgen.tcl
pause
