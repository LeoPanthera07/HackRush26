@echo off
cd /d "%~dp0"
echo Running VMD psfgen for 1M3D_C_2001_LU...
vmd -dispdev text -e psfgen.tcl
pause
