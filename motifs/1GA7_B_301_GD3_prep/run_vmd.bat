@echo off
:: HackRush '26 | VMD psfgen runner
:: Make sure VMD is installed and in PATH
:: Or set VMD_EXE to full path below, e.g.: C:\Program Files\VMD\vmd.exe

set VMD_EXE=vmd

cd /d "%~dp0"
echo Running VMD psfgen...
"%VMD_EXE%" -dispdev text -e psfgen.tcl

if %ERRORLEVEL% EQU 0 (
    echo.
    echo SUCCESS — PSF and PDB generated:
    echo   1GA7_B_301_GD3_final.psf
    echo   1GA7_B_301_GD3_final.pdb
    echo.
    echo Load in VMD:
    echo   vmd 1GA7_B_301_GD3_final.pdb -psf 1GA7_B_301_GD3_final.psf
) else (
    echo.
    echo ERROR — check psfgen.tcl and prep.log for issues.
)
pause
