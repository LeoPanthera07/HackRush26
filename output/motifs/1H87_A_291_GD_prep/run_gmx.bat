@echo off
rem GROMACS runner for 1H87_A_291_GD
cd /d "%~dp0"

echo [1/4] Running pdb2gmx...
gmx pdb2gmx -f protein_gmx.pdb -o 1H87_A_291_GD_protein.gro -p 1H87_A_291_GD_protein.top -i posre.itp -ff charmm36 -water tip3p -ter -ignh

echo [2/4] Setting box...
gmx editconf -f 1H87_A_291_GD_protein.gro -o 1H87_A_291_GD_box.gro -c -d 1.0 -bt cubic

echo [3/4] Energy minimisation...
gmx grompp -f em.mdp -c 1H87_A_291_GD_box.gro -p topol.top -o em.tpr
gmx mdrun  -v -deffnm em

echo [4/4] Done. Check em.log for convergence.
pause
