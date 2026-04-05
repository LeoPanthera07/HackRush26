@echo off
cd /d "%~dp0"
echo Running VMD psfgen for 1M3D_L_2004_LU...
vmd -dispdev text -e psfgen.tcl
pause
