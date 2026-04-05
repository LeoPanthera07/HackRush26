# psfgen.tcl - 1M3D_F_2002_LU
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
segment PROC {
    pdb protein_C_charmm.pdb
    first NTER
    last  CTER
}
coordpdb protein_C_charmm.pdb PROC

segment PROF {
    pdb protein_F_charmm.pdb
    first NTER
    last  CTER
}
coordpdb protein_F_charmm.pdb PROF


# --- Water segment ---
segment SOLV {
    auto none
    pdb water_charmm.pdb
}
coordpdb water_charmm.pdb SOLV

# --- Lanthanide LU ---
segment LN01 {
    auto none
    pdb ions_charmm.pdb
}
coordpdb ions_charmm.pdb LN01

guesscoord
writepsf 1M3D_F_2002_LU_vmd.psf
writepdb 1M3D_F_2002_LU_vmd.pdb
puts "\nDone: 1M3D_F_2002_LU_vmd.psf and 1M3D_F_2002_LU_vmd.pdb written."
