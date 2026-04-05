# psfgen.tcl - 1H87_A_291_GD
# VMD CHARMM36 PSF/PDB generation
# Run: vmd -dispdev text -e psfgen.tcl
#
# Place CHARMM36 topology files in ./charmm36/ before running.
# Required: top_all36_prot.rtf, toppar/toppar_water_ions.str

package require psfgen

set CHARMMDIR "./charmm36"
topology $CHARMMDIR/top_all36_prot.rtf
topology $CHARMMDIR/toppar/toppar_water_ions.str
topology ./ln_custom.rtf

# --- Protein segments ---

# --- Water segment ---
segment SOLV {
    auto none
    pdb water_charmm.pdb
}
coordpdb water_charmm.pdb SOLV

# --- Lanthanide GD ---
segment LN01 {
    auto none
    pdb ions_charmm.pdb
}
coordpdb ions_charmm.pdb LN01

guesscoord
writepsf 1H87_A_291_GD_vmd.psf
writepdb 1H87_A_291_GD_vmd.pdb
puts "\nDone: 1H87_A_291_GD_vmd.psf and 1H87_A_291_GD_vmd.pdb written."
