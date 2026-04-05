@echo off
cd /d "%~dp0"
echo Running VMD psfgen for 1H87_A_291_GD...
vmd -dispdev text -e psfgen.tcl
pause
