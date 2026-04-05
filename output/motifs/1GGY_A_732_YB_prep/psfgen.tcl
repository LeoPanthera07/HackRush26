# psfgen.tcl - 1GGY_A_732_YB
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
segment PROA {
    pdb protein_A_charmm.pdb
    first NTER
    last  CTER
}
coordpdb protein_A_charmm.pdb PROA


# --- Water segment ---
segment SOLV {
    auto none
    pdb water_charmm.pdb
}
coordpdb water_charmm.pdb SOLV

# --- Lanthanide YB ---
segment LN01 {
    auto none
    pdb ions_charmm.pdb
}
coordpdb ions_charmm.pdb LN01

guesscoord
writepsf 1GGY_A_732_YB_vmd.psf
writepdb 1GGY_A_732_YB_vmd.pdb
puts "\nDone: 1GGY_A_732_YB_vmd.psf and 1GGY_A_732_YB_vmd.pdb written."
