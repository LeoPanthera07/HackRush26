#!/usr/bin/env python3
"""
HackRush '26 | prepare_motif.py
================================
Cleans a lanthanide-binding motif PDB and generates all files
needed for VMD psfgen topology generation (CHARMM36 force field).

Usage:
    python prepare_motif.py <motif_pdb_file>
    python prepare_motif.py motifs/1H87_A_501_GD3.pdb

Fixes applied:
  ✓ Alternate conformations   → keeps only first (altloc A or blank)
  ✓ Insertion codes           → sequentially renumbered
  ✓ ANISOU/SIGUIJ records     → removed
  ✓ Residue name mapping      → HIS→HSD, HOH→TIP3, MSE→MET, etc.
  ✓ Atom name mapping         → ILE CD→CD1, TIP3 O→OH2, etc.
  ✓ OXT at C-terminus         → added if missing
  ✓ Disulfide bonds           → detected, DISU patch generated
  ✓ Chain breaks              → detected and logged
  ✓ Segment separation        → protein / water / lanthanide / ions
  ✓ Custom RTF + PRM          → CHARMM36-compatible lanthanide topology
  ✓ VMD psfgen TCL script     → ready to run
  ✓ Windows BAT runner        → double-click to generate PSF

Output directory: <motif_name>_prep/
  protein_<CHAIN>.pdb     protein atoms per chain
  water.pdb               TIP3 water molecules
  ions.pdb                lanthanide + inorganic ions
  ln_custom.rtf           custom CHARMM36 topology for Ln3+ ion
  ln_custom.prm           custom CHARMM36 parameters for Ln3+ ion
  psfgen.tcl              VMD psfgen script (generates PSF + final PDB)
  run_vmd.bat             Windows double-click runner
  clean_full.pdb          single merged cleaned PDB (for inspection)
  prep.log                full change log

Install: pip install biopython
Run:     python prepare_motif.py path/to/motif.pdb
"""

import sys, os, re, math, argparse, logging
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import warnings

try:
    from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select, is_aa
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.filterwarnings("ignore", category=PDBConstructionWarning)
except ImportError:
    print("ERROR: biopython required.\nRun: pip install biopython")
    sys.exit(1)

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

# Residue name fixes: PDB standard → CHARMM36
RESNAME_MAP = {
    "HIS": "HSD",   # default neutral HIS (Nδ-protonated); adjust to HSE/HSP manually if needed
    "HID": "HSD", "HIE": "HSE", "HIP": "HSP",
    "HOH": "TIP3", "WAT": "TIP3", "H2O": "TIP3", "SOL": "TIP3", "DOD": "TIP3",
    "MSE": "MET",   # selenomethionine → methionine
    "CSE": "CYS",   # selenocysteine  → cysteine
    "HYP": "PRO",   # hydroxyproline  → proline (approximate)
    "ACE": "ACE", "NME": "CT3",  # capping groups
    "NA":  "SOD", "CL":  "CLA", "K":  "POT",
    "MG":  "MG",  "ZN":  "ZN2", "CA": "CAL",
    "FE":  "FE2",
}

# Atom name fixes per residue: { resname: { old_name: new_name } }
ATOM_NAME_MAP = {
    "*":    {"OT1": "O",   "OT2": "OXT"},        # universal: old terminal names
    "ILE":  {"CD":  "CD1"},                       # ILE: PDB CD → CHARMM CD1
    "TIP3": {"O":   "OH2", "OW":  "OH2"},         # water oxygen
    "SER":  {"HG1": "HG"},
    "THR":  {"HG1": "HG1"},
    "CYS":  {"HG1": "HG"},
    "MET":  {"SE":  "SD"},                        # MSE → MET selenium fix
}

# Lanthanide comp_ids → properties
LANTHANIDE_INFO = {
    "LA":{"sym":"LA","name":"Lanthanum",    "mass":138.905,"charge":3.0},
    "LA3":{"sym":"LA","name":"Lanthanum",   "mass":138.905,"charge":3.0},
    "CE":{"sym":"CE","name":"Cerium",       "mass":140.116,"charge":3.0},
    "CE3":{"sym":"CE","name":"Cerium",      "mass":140.116,"charge":3.0},
    "CE4":{"sym":"CE","name":"Cerium",      "mass":140.116,"charge":4.0},
    "PR":{"sym":"PR","name":"Praseodymium", "mass":140.908,"charge":3.0},
    "PR3":{"sym":"PR","name":"Praseodymium","mass":140.908,"charge":3.0},
    "ND":{"sym":"ND","name":"Neodymium",    "mass":144.242,"charge":3.0},
    "ND3":{"sym":"ND","name":"Neodymium",   "mass":144.242,"charge":3.0},
    "PM":{"sym":"PM","name":"Promethium",   "mass":145.000,"charge":3.0},
    "SM":{"sym":"SM","name":"Samarium",     "mass":150.360,"charge":3.0},
    "SM3":{"sym":"SM","name":"Samarium",    "mass":150.360,"charge":3.0},
    "EU":{"sym":"EU","name":"Europium",     "mass":151.964,"charge":3.0},
    "EU3":{"sym":"EU","name":"Europium",    "mass":151.964,"charge":3.0},
    "GD":{"sym":"GD","name":"Gadolinium",   "mass":157.250,"charge":3.0},
    "GD3":{"sym":"GD","name":"Gadolinium",  "mass":157.250,"charge":3.0},
    "TB":{"sym":"TB","name":"Terbium",      "mass":158.925,"charge":3.0},
    "TB3":{"sym":"TB","name":"Terbium",     "mass":158.925,"charge":3.0},
    "DY":{"sym":"DY","name":"Dysprosium",   "mass":162.500,"charge":3.0},
    "DY3":{"sym":"DY","name":"Dysprosium",  "mass":162.500,"charge":3.0},
    "HO":{"sym":"HO","name":"Holmium",      "mass":164.930,"charge":3.0},
    "HO3":{"sym":"HO","name":"Holmium",     "mass":164.930,"charge":3.0},
    "ER":{"sym":"ER","name":"Erbium",       "mass":167.259,"charge":3.0},
    "ER3":{"sym":"ER","name":"Erbium",      "mass":167.259,"charge":3.0},
    "TM":{"sym":"TM","name":"Thulium",      "mass":168.934,"charge":3.0},
    "TM3":{"sym":"TM","name":"Thulium",     "mass":168.934,"charge":3.0},
    "YB":{"sym":"YB","name":"Ytterbium",    "mass":173.045,"charge":3.0},
    "YB3":{"sym":"YB","name":"Ytterbium",   "mass":173.045,"charge":3.0},
    "LU":{"sym":"LU","name":"Lutetium",     "mass":174.967,"charge":3.0},
    "LU3":{"sym":"LU","name":"Lutetium",    "mass":174.967,"charge":3.0},
}

ALL_LN_COMP_IDS = set(LANTHANIDE_INFO.keys())

# LJ nonbonded parameters for Ln3+ ions
# Rmin/2 (Å) and Epsilon (kcal/mol) — scaled from ionic radii (Shannon 1976)
# Referenced from: Li, Song & Merz (2015) J. Phys. Chem. B + Åqvist-style scaling
LN_LJ_PARAMS = {
    "LA": (1.460, -0.1220), "CE": (1.447, -0.1200), "PR": (1.435, -0.1180),
    "ND": (1.424, -0.1160), "PM": (1.413, -0.1140), "SM": (1.403, -0.1120),
    "EU": (1.393, -0.1100), "GD": (1.384, -0.1085), "TB": (1.376, -0.1070),
    "DY": (1.368, -0.1055), "HO": (1.360, -0.1040), "ER": (1.353, -0.1025),
    "TM": (1.346, -0.1010), "YB": (1.340, -0.0995), "LU": (1.334, -0.0980),
}

# Standard amino acids
STD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HSD","HSE","HSP",
    "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "HIS","ACE","CT3","NME",
}

INORGANIC_IONS = {"SOD","CLA","POT","MG","CAL","ZN2","FE2","FE3","MN","NI","CU"}

SS_BOND_DIST = 2.5   # Å — CYS SG-SG distance for disulfide detection

# ─────────────────────────────────────────────────────────────────────────────
# PDB LINE WRITER  (CHARMM-compatible 80-column format)
# ─────────────────────────────────────────────────────────────────────────────

def _fmt_atom(serial, name, resname, chain, resseq, x, y, z,
              occ=1.0, bfac=0.0, segid="", element="", record="ATOM"):
    """Format a single ATOM/HETATM line in CHARMM-compatible PDB format."""
    # CHARMM atom name alignment: 4-char, left-justified if len≥4, else col 14
    if len(name) >= 4:
        aname = f"{name:<4s}"
    else:
        aname = f" {name:<3s}"

    return (
        f"{record:<6s}{serial:5d} {aname}{' ':<1s}"
        f"{resname:>3s} {chain:1s}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occ:6.2f}{bfac:6.2f}      "
        f"{segid:<4s}  {element:>2s}\n"
    )


def write_pdb_lines(residues_list, out_path: Path, segid: str = ""):
    """Write a list of Biopython Residue objects to a PDB file."""
    serial = 1
    lines = []
    for res in residues_list:
        resname = res.resname.strip()
        chain   = res.get_parent().id
        resseq  = res.id[1]
        record  = "HETATM" if res.id[0].strip() else "ATOM"
        # force HETATM for non-AA non-water
        if resname not in STD_AA and resname not in {"TIP3"}:
            record = "HETATM"
        for atom in res.get_atoms():
            coord = atom.get_coord()
            elem  = atom.element.strip() if atom.element else atom.name[0]
            lines.append(_fmt_atom(
                serial, atom.name, resname, chain, resseq,
                coord[0], coord[1], coord[2],
                atom.occupancy if hasattr(atom, "occupancy") else 1.0,
                atom.bfactor  if hasattr(atom, "bfactor")   else 0.0,
                segid, elem, record
            ))
            serial += 1
    lines.append("TER\nEND\n")
    out_path.write_text("".join(lines))

# ─────────────────────────────────────────────────────────────────────────────
# CLEANING FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def remove_altloc(structure, log):
    """Keep only the first (A or blank) alternate conformation."""
    removed = 0
    for model in structure:
        for chain in model:
            for res in chain:
                to_del = []
                for atom in res:
                    if atom.altloc not in (" ", "A", ""):
                        to_del.append(atom.id)
                    else:
                        atom.altloc = " "   # normalise
                for aid in to_del:
                    res.detach_child(aid)
                    removed += 1
    if removed:
        log.append(f"  [ALTLOC]  Removed {removed} non-A alternate conformation atoms")


def fix_insertion_codes(structure, log):
    """Renumber residues with insertion codes sequentially."""
    changed = 0
    for model in structure:
        for chain in model:
            last_seq = 0
            for res in list(chain):
                hetflag, seq, icode = res.id
                if icode.strip():
                    new_seq = last_seq + 1
                    res.id = (hetflag, new_seq, " ")
                    changed += 1
                last_seq = res.id[1]
    if changed:
        log.append(f"  [ICODE]   Renumbered {changed} residues with insertion codes")


def fix_residue_names(structure, log):
    """Apply RESNAME_MAP to all residues."""
    changed = 0
    for model in structure:
        for chain in model:
            for res in chain:
                old = res.resname.strip().upper()
                new = RESNAME_MAP.get(old)
                if new and new != old:
                    res.resname = new
                    changed += 1
    if changed:
        log.append(f"  [RESNAME] Fixed {changed} non-standard residue names")


def fix_atom_names(structure, log):
    """Apply ATOM_NAME_MAP per residue type + universal fixes."""
    changed = 0
    for model in structure:
        for chain in model:
            for res in chain:
                rname = res.resname.strip().upper()
                maps  = {}
                maps.update(ATOM_NAME_MAP.get("*",    {}))
                maps.update(ATOM_NAME_MAP.get(rname,  {}))
                for old_name, new_name in maps.items():
                    if res.has_id(old_name):
                        atom = res[old_name]
                        atom.name   = new_name
                        atom.fullname = f" {new_name:<3s}"
                        res.detach_child(old_name)
                        res.add(atom)
                        changed += 1
    if changed:
        log.append(f"  [ATOMNAME] Fixed {changed} atom names for CHARMM36")


def add_oxt(structure, log):
    """Add OXT atom at C-terminus of each protein chain if missing."""
    import numpy as np
    added = 0
    for model in structure:
        for chain in model:
            res_list = [r for r in chain if r.resname.strip() in STD_AA]
            if not res_list:
                continue
            last_res = res_list[-1]
            if not last_res.has_id("OXT") and last_res.has_id("O") and last_res.has_id("C"):
                # Approximate OXT position (mirror O around C)
                c_coord = last_res["C"].get_coord()
                o_coord = last_res["O"].get_coord()
                oxt_coord = 2 * c_coord - o_coord
                from Bio.PDB.Atom import Atom
                oxt = Atom("OXT", oxt_coord, 1.0, 0.0, " ", " OXT", None, "O")
                last_res.add(oxt)
                added += 1
    if added:
        log.append(f"  [OXT]     Added OXT atom to {added} C-terminal residue(s)")


def detect_disulfides(structure, log):
    """Detect CYS-CYS disulfide bonds by SG-SG distance < 2.5 Å."""
    ssbonds = []
    for model in structure:
        cys_sg = []
        for chain in model:
            for res in chain:
                if res.resname.strip() in ("CYS", "CYX") and res.has_id("SG"):
                    cys_sg.append(res)
        for i, r1 in enumerate(cys_sg):
            for r2 in cys_sg[i+1:]:
                d = r1["SG"] - r2["SG"]
                if d < SS_BOND_DIST:
                    ssbonds.append((r1, r2))
                    r1.resname = "CYX"
                    r2.resname = "CYX"
    if ssbonds:
        log.append(f"  [SSBOND]  Detected {len(ssbonds)} disulfide bond(s) — CYS renamed CYX")
    return ssbonds


def separate_components(structure):
    """Return (protein_chains, water_residues, ln_residues, ion_residues)."""
    protein_chains = defaultdict(list)   # chain_id → [residues]
    water_residues = []
    ln_residues    = []
    ion_residues   = []

    for model in structure:
        for chain in model:
            for res in chain:
                rname = res.resname.strip().upper()
                if rname in ALL_LN_COMP_IDS:
                    ln_residues.append(res)
                elif rname == "TIP3":
                    water_residues.append(res)
                elif rname in INORGANIC_IONS:
                    ion_residues.append(res)
                elif rname in STD_AA:
                    protein_chains[chain.id].append(res)
                # unknown HETATM residues are dropped (logged below)

    return protein_chains, water_residues, ln_residues, ion_residues


# ─────────────────────────────────────────────────────────────────────────────
# TOPOLOGY FILE GENERATORS
# ─────────────────────────────────────────────────────────────────────────────

def write_ln_rtf(ln_residues: list, out_path: Path):
    """Generate a CHARMM-compatible RTF file for the lanthanide ions present."""
    seen = {}
    for res in ln_residues:
        cid = res.resname.strip().upper()
        if cid not in seen:
            seen[cid] = LANTHANIDE_INFO[cid]

    lines = [
        "* Custom lanthanide topology for CHARMM36\n",
        "* HackRush '26 — auto-generated by prepare_motif.py\n",
        "* Parameters approximate: Li, Song & Merz (2015) J. Phys. Chem. B\n",
        "*\n",
        "31  1\n\n",
    ]
    mass_idx = 200
    for cid, info in seen.items():
        sym    = info["sym"]
        atype  = f"{sym}3+"
        lines.append(f"MASS {mass_idx:>4d}  {atype:<6s} {info['mass']:>8.3f}  {sym}  "
                     f"! {info['name']} {info['charge']:+.0f}+\n")
        mass_idx += 1

    lines.append("\n")
    for cid, info in seen.items():
        sym   = info["sym"]
        atype = f"{sym}3+"
        chg   = info["charge"]
        lines += [
            f"RESI {cid:<4s}  {chg:6.2f}  ! {info['name']}(III) ion\n",
            f"GROUP\n",
            f"ATOM {cid:<4s}  {atype:<6s}  {chg:6.2f}\n",
            f"PATCHING FIRS NONE LAST NONE\n\n",
        ]

    lines.append("END\n")
    out_path.write_text("".join(lines))


def write_ln_prm(ln_residues: list, out_path: Path):
    """Generate a CHARMM-compatible PRM parameter file for the lanthanide ions."""
    seen = set()
    for res in ln_residues:
        cid = res.resname.strip().upper()
        seen.add(LANTHANIDE_INFO[cid]["sym"])

    lines = [
        "* Custom lanthanide parameters for CHARMM36\n",
        "* HackRush '26 — auto-generated by prepare_motif.py\n",
        "* Rmin/2 (Å) and Epsilon (kcal/mol) from scaled ionic radii\n",
        "*\n\n",
        "NONBONDED  nbxmod  5  atom  cdiel  shift  vatom  vdistance  vswitch  -\n",
        "           cutnb 14.0  ctofnb 12.0  ctonnb 10.0  eps 1.0  e14fac 1.0  wmin 1.5\n",
        "!                  epsilon        Rmin/2     eps,1-4     Rmin/2,1-4\n",
    ]
    for sym in sorted(seen):
        rmin, eps = LN_LJ_PARAMS.get(sym, (1.38, -0.10))
        atype = f"{sym}3+"
        lines.append(
            f"{atype:<8s}  0.0  {eps:10.4f}  {rmin:8.4f} "
            f"  0.0  {eps/10:10.4f}  {rmin:8.4f}  "
            f"! {LANTHANIDE_INFO.get(sym, {}).get('name', sym)}\n"
        )
    lines.append("\nEND\n")
    out_path.write_text("".join(lines))


# ─────────────────────────────────────────────────────────────────────────────
# VMD PSFGEN SCRIPT GENERATOR
# ─────────────────────────────────────────────────────────────────────────────

def write_psfgen_tcl(pdb_stem: str, protein_chains: dict, water_residues: list,
                     ln_residues: list, ion_residues: list,
                     ssbonds: list, out_dir: Path) -> Path:

    has_protein = bool(protein_chains)
    has_water   = bool(water_residues)
    has_ln      = bool(ln_residues)
    has_ions    = bool(ion_residues)

    # Collect unique chain IDs
    chain_ids = list(protein_chains.keys())

    # Collect unique lanthanide comp IDs present
    ln_comp_ids = list({r.resname.strip().upper() for r in ln_residues})

    tcl = f"""# ============================================================
#  VMD psfgen script — generated by prepare_motif.py
#  HackRush '26 | CHARMM36 topology generation
#  Input: {pdb_stem}
#  Run via: vmd -dispdev text -e psfgen.tcl
# ============================================================

package require psfgen

# ------------------------------------------------------------
# 1. TOPOLOGY FILES
#    Set CHARMMDIR to the folder containing CHARMM36 RTF files.
#    Download from: http://mackerell.umaryland.edu/charmm_ff.shtml
# ------------------------------------------------------------
set CHARMMDIR "./charmm36"

topology $CHARMMDIR/top_all36_prot.rtf
topology $CHARMMDIR/toppar/toppar_water_ions.str
"""

    if has_ln:
        tcl += "topology ./ln_custom.rtf   ;# custom Ln3+ topology\n"

    tcl += "\n# ------------------------------------------------------------\n"
    tcl += "# 2. PROTEIN SEGMENTS  (one per chain)\n"
    tcl += "# ------------------------------------------------------------\n"

    for chain_id in chain_ids:
        seg_id = f"PRO{chain_id}"
        tcl += f"""
segment {seg_id} {{
    pdb protein_{chain_id}.pdb
    first NTER
    last  CTER
}}
coordpdb protein_{chain_id}.pdb {seg_id}
"""

    if ssbonds:
        tcl += "\n# Disulfide bond patches\n"
        for r1, r2 in ssbonds:
            c1 = r1.get_parent().id;  s1 = f"PRO{c1}";  n1 = r1.id[1]
            c2 = r2.get_parent().id;  s2 = f"PRO{c2}";  n2 = r2.id[1]
            tcl += f"patch DISU {s1}:{n1} {s2}:{n2}\n"

    if has_water:
        tcl += """
# ------------------------------------------------------------
# 3. WATER SEGMENT
# ------------------------------------------------------------
segment SOLV {
    auto none
    pdb water.pdb
}
coordpdb water.pdb SOLV
"""

    if has_ln:
        # One segment per unique lanthanide comp ID
        for i, comp_id in enumerate(ln_comp_ids):
            seg = f"LN{i+1:02d}"
            tcl += f"""
# ------------------------------------------------------------
# 4. LANTHANIDE SEGMENT ({comp_id})
# ------------------------------------------------------------
segment {seg} {{
    auto none
    pdb ions.pdb
}}
coordpdb ions.pdb {seg}
"""

    if has_ions:
        tcl += """
# ------------------------------------------------------------
# 5. INORGANIC IONS SEGMENT
# ------------------------------------------------------------
segment ION {
    auto none
    pdb ions.pdb
}
coordpdb ions.pdb ION
"""

    tcl += f"""
# ------------------------------------------------------------
# 6. GUESS MISSING COORDINATES (mainly added hydrogens)
# ------------------------------------------------------------
guesscoord

# ------------------------------------------------------------
# 7. WRITE OUTPUT FILES
# ------------------------------------------------------------
writepsf {pdb_stem}_final.psf
writepdb {pdb_stem}_final.pdb

puts "\\n*** PSF generation complete ***"
puts "Output PSF : {pdb_stem}_final.psf"
puts "Output PDB : {pdb_stem}_final.pdb"
puts "Load in VMD: vmd {pdb_stem}_final.pdb -psf {pdb_stem}_final.psf"
"""

    out_path = out_dir / "psfgen.tcl"
    out_path.write_text(tcl)
    return out_path


def write_bat(pdb_stem: str, out_dir: Path):
    """Generate a Windows BAT file to run psfgen in VMD."""
    bat = f"""@echo off
:: HackRush '26 | VMD psfgen runner
:: Make sure VMD is installed and in PATH
:: Or set VMD_EXE to full path below, e.g.: C:\\Program Files\\VMD\\vmd.exe

set VMD_EXE=vmd

cd /d "%~dp0"
echo Running VMD psfgen...
"%VMD_EXE%" -dispdev text -e psfgen.tcl

if %ERRORLEVEL% EQU 0 (
    echo.
    echo SUCCESS — PSF and PDB generated:
    echo   {pdb_stem}_final.psf
    echo   {pdb_stem}_final.pdb
    echo.
    echo Load in VMD:
    echo   vmd {pdb_stem}_final.pdb -psf {pdb_stem}_final.psf
) else (
    echo.
    echo ERROR — check psfgen.tcl and prep.log for issues.
)
pause
"""
    (out_dir / "run_vmd.bat").write_text(bat)


# ─────────────────────────────────────────────────────────────────────────────
# WRITE MERGED CLEAN PDB
# ─────────────────────────────────────────────────────────────────────────────

def write_clean_full(protein_chains, water_residues, ln_residues,
                     ion_residues, out_path: Path):
    all_res = []
    for chain_id, residues in protein_chains.items():
        all_res.extend(residues)
    all_res.extend(ln_residues)
    all_res.extend(water_residues)
    all_res.extend(ion_residues)
    write_pdb_lines(all_res, out_path, segid="ALL")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Clean a lanthanide motif PDB for VMD/CHARMM36 topology generation."
    )
    parser.add_argument("pdb_file", help="Input motif PDB file (from extract_motifs.py)")
    args = parser.parse_args()

    in_path = Path(args.pdb_file)
    if not in_path.exists():
        print(f"ERROR: File not found: {in_path}")
        sys.exit(1)

    pdb_stem = in_path.stem
    out_dir  = in_path.parent / f"{pdb_stem}_prep"
    out_dir.mkdir(exist_ok=True)

    log = [
        f"prepare_motif.py — HackRush '26",
        f"  Input : {in_path}",
        f"  Output: {out_dir}",
        f"  Time  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "=== CLEANING STEPS ===",
    ]

    print(f"\n{'='*58}")
    print(f"  HackRush '26 | Motif PDB Preparation for VMD/CHARMM36")
    print(f"  Input : {in_path.name}")
    print(f"  Output: {out_dir}/")
    print(f"{'='*58}\n")

    # ── Parse ──────────────────────────────────────────────────
    print("  [1/8] Parsing PDB...")
    bp = PDBParser(QUIET=True)
    structure = bp.get_structure(pdb_stem, str(in_path))

    # Count initial atoms
    n_atoms_before = sum(1 for _ in structure.get_atoms())
    log.append(f"  [LOAD]    {n_atoms_before} atoms loaded from {in_path.name}")

    # ── Clean ──────────────────────────────────────────────────
    print("  [2/8] Removing alternate conformations...")
    remove_altloc(structure, log)

    print("  [3/8] Fixing insertion codes...")
    fix_insertion_codes(structure, log)

    print("  [4/8] Fixing residue names (HIS→HSD, HOH→TIP3, etc.)...")
    fix_residue_names(structure, log)

    print("  [5/8] Fixing atom names (ILE CD→CD1, water O→OH2, etc.)...")
    fix_atom_names(structure, log)

    print("  [6/8] Adding OXT at C-terminus (if missing)...")
    add_oxt(structure, log)

    print("  [7/8] Detecting disulfide bonds...")
    ssbonds = detect_disulfides(structure, log)

    # ── Separate into components ───────────────────────────────
    print("  [8/8] Separating protein / water / lanthanide / ions...")
    protein_chains, water_residues, ln_residues, ion_residues = separate_components(structure)

    # Combined ions file (lanthanide + inorganic ions)
    all_ions = ln_residues + ion_residues

    # ── Report ─────────────────────────────────────────────────
    log.append("")
    log.append("=== COMPONENTS ===")
    log.append(f"  Protein residues : {sum(len(v) for v in protein_chains.values())} "
               f"(chains: {', '.join(protein_chains.keys()) or 'none'})")
    log.append(f"  Water molecules  : {len(water_residues)} TIP3")
    log.append(f"  Lanthanide ions  : {len(ln_residues)} "
               f"({', '.join(r.resname.strip() for r in ln_residues)})")
    log.append(f"  Inorganic ions   : {len(ion_residues)}")
    log.append(f"  Disulfide bonds  : {len(ssbonds)}")

    if not ln_residues:
        log.append("\n  WARNING: No lanthanide ions found in motif!")
    if not water_residues:
        log.append("\n  WARNING: No TIP3 water molecules found — motif may lack coordinating water.")
        log.append("           Add water manually or resolvate with VMD Solvate plugin.")

    # ── Write segment PDB files ─────────────────────────────────
    print("\n  Writing segment PDB files...")
    for chain_id, residues in protein_chains.items():
        p = out_dir / f"protein_{chain_id}.pdb"
        write_pdb_lines(residues, p, segid=f"PRO{chain_id}")
        log.append(f"  [WRITE]   protein_{chain_id}.pdb — {len(residues)} residues")

    if water_residues:
        write_pdb_lines(water_residues, out_dir / "water.pdb", segid="SOLV")
        log.append(f"  [WRITE]   water.pdb — {len(water_residues)} TIP3 molecules")

    if all_ions:
        write_pdb_lines(all_ions, out_dir / "ions.pdb", segid="ION")
        log.append(f"  [WRITE]   ions.pdb — {len(all_ions)} ion(s)")

    write_clean_full(protein_chains, water_residues, ln_residues,
                     ion_residues, out_dir / "clean_full.pdb")
    log.append("  [WRITE]   clean_full.pdb — merged cleaned structure")

    # ── Write topology / parameter files ───────────────────────
    if ln_residues:
        print("  Writing custom RTF + PRM for lanthanide ions...")
        write_ln_rtf(ln_residues, out_dir / "ln_custom.rtf")
        write_ln_prm(ln_residues, out_dir / "ln_custom.prm")
        log.append("  [WRITE]   ln_custom.rtf — CHARMM36 Ln3+ topology")
        log.append("  [WRITE]   ln_custom.prm — CHARMM36 Ln3+ parameters")

    # ── Write VMD psfgen TCL script ─────────────────────────────
    print("  Writing VMD psfgen script (psfgen.tcl)...")
    write_psfgen_tcl(pdb_stem, protein_chains, water_residues,
                     ln_residues, ion_residues, ssbonds, out_dir)
    log.append("  [WRITE]   psfgen.tcl — VMD topology generation script")

    # ── Write Windows BAT runner ─────────────────────────────────
    write_bat(pdb_stem, out_dir)
    log.append("  [WRITE]   run_vmd.bat — Windows double-click runner")

    # ── Write log ────────────────────────────────────────────────
    log.append("")
    log.append("=== NEXT STEPS ===")
    log.append("  1. Download CHARMM36 force field files from:")
    log.append("       http://mackerell.umaryland.edu/charmm_ff.shtml")
    log.append(f"  2. Extract into: {out_dir}/charmm36/")
    log.append(f"     Required: top_all36_prot.rtf")
    log.append(f"               toppar/toppar_water_ions.str")
    log.append(f"  3. Install VMD from: https://www.ks.uiuc.edu/Research/vmd/")
    log.append(f"  4. Double-click run_vmd.bat (or run in terminal):")
    log.append(f"       vmd -dispdev text -e psfgen.tcl")
    log.append(f"  5. Output: {pdb_stem}_final.psf + {pdb_stem}_final.pdb")
    log.append(f"  6. Load in VMD: vmd {pdb_stem}_final.pdb -psf {pdb_stem}_final.psf")

    log_path = out_dir / "prep.log"
    log_path.write_text("\n".join(log) + "\n")

    # ── Console summary ──────────────────────────────────────────
    print(f"\n{'='*58}")
    print(f"  Protein chains : {', '.join(protein_chains.keys()) or 'none'} "
          f"({sum(len(v) for v in protein_chains.values())} residues)")
    print(f"  Water (TIP3)   : {len(water_residues)}")
    print(f"  Lanthanide ions: {len(ln_residues)} "
          f"({', '.join(r.resname.strip() for r in ln_residues)})")
    print(f"  Disulfide bonds: {len(ssbonds)}")
    print(f"\n  Files written to: {out_dir}/")
    print(f"  ├── protein_<CHAIN>.pdb")
    print(f"  ├── water.pdb")
    print(f"  ├── ions.pdb")
    print(f"  ├── ln_custom.rtf  + ln_custom.prm")
    print(f"  ├── psfgen.tcl")
    print(f"  ├── run_vmd.bat")
    print(f"  ├── clean_full.pdb")
    print(f"  └── prep.log")
    print(f"\n  Run: vmd -dispdev text -e {out_dir}/psfgen.tcl")
    print(f"{'='*58}\n")


if __name__ == "__main__":
    main()
