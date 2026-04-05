#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
end_to_end_motif_pipeline.py
HackRush 26 - Problem 9B
Lanthanide Binding Motif Pipeline: VMD (CHARMM36) + GROMACS

PIPELINE STAGES
---------------
Stage 1  Read IDs -> Download PDB -> Extract motifs
  - Coordination shell : all residues within 4.0 A of Ln3+
  - Water requirement  : at least 1 HOH within 3.0 A (skip site if absent)
  - Output             : motifs/<PDB>_<CH>_<RNUM>_<CID>.pdb + motif_summary.csv

Stage 2  Clean each motif
  - Remove non-A altlocs
  - Fix insertion codes (renumber)
  - CHARMM36 residue names (HIS->HSD, HOH->TIP3, etc.)
  - Fix atom names (ILE CD->CD1, water O->OH2, OT1->O, OT2->OXT)
  - Add missing OXT at C-terminus
  - Detect disulfides (SG-SG < 2.5 A) -> label CYX
  - Split: protein chains / water / lanthanide / other ions

Stage 3a  VMD / CHARMM36 output
  - protein_<CH>_charmm.pdb    CHARMM36 named protein segment
  - water_charmm.pdb            TIP3 water segment
  - ions_charmm.pdb             Ln3+ + inorganic ions
  - ln_custom.rtf               MASS + RESI blocks for CHARMM36
  - ln_custom.prm               NONBONDED LJ params (Rmin/2, kcal/mol)
  - psfgen.tcl                  VMD psfgen script -> writes .psf + .pdb
  - run_vmd.bat                 Windows one-click VMD runner

Stage 3b  GROMACS / CHARMM36 output
  - protein_gmx.pdb             HETATM-free, HIS names for pdb2gmx
  - ffnonbonded_ln.itp          [ atomtypes ] sigma/nm + eps/kJ/mol
  - ln_<CID>.itp                [ moleculetype ] + [ atoms ] per Ln
  - posre.itp                   position restraints (protein heavy atoms)
  - topol.top                   master topology with all #include entries
  - em.mdp                      energy minimisation parameters
  - md.mdp                      production MD parameters (100 ps template)
  - run_gmx.bat                 Windows GROMACS runner
  - run_gmx.sh                  Linux / WSL GROMACS runner
  - clean_full.pdb              merged motif PDB for visual inspection
  - prep.log                    log of all changes made

Usage
-----
  pip install biopython requests

  python end_to_end_motif_pipeline.py
  python end_to_end_motif_pipeline.py --ids 1H87 2XU3
  python end_to_end_motif_pipeline.py --target vmd
  python end_to_end_motif_pipeline.py --target gromacs
  python end_to_end_motif_pipeline.py --skip-extract
"""

import os
import sys
import csv
import time
import argparse
import warnings
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import Optional, List, Dict, Tuple

import requests

try:
    from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio.PDB.Atom import Atom as BioAtom
except ImportError:
    sys.exit("ERROR: Install biopython first:  pip install biopython requests")

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# ---------------------------------------------------------------------------
# SECTION 1 - CONSTANTS
# ---------------------------------------------------------------------------

SCRIPT_DIR    = Path(os.path.dirname(os.path.abspath(__file__)))
IDS_FILE      = SCRIPT_DIR / "lanthanide_pdb_ids.txt"
PDB_DIR       = SCRIPT_DIR / "pdb_files"
MOTIF_DIR     = SCRIPT_DIR / "motifs"
MOTIF_SUMMARY = SCRIPT_DIR / "motif_summary.csv"
RCSB_URL      = "https://files.rcsb.org/download/{}.pdb"

COORD_CUTOFF = 4.0   # Angstrom - coordination shell radius
WATER_CUTOFF = 3.0   # Angstrom - water proximity requirement
SS_CUTOFF    = 2.5   # Angstrom - SG-SG disulfide threshold

# Lanthanide data: comp_id -> {sym, name, mass, charge, atnum}
_LN_RAW = [
    ("LA",  "LA", "Lanthanum",    138.905, 3.0, 57),
    ("LA3", "LA", "Lanthanum",    138.905, 3.0, 57),
    ("CE",  "CE", "Cerium",       140.116, 3.0, 58),
    ("CE3", "CE", "Cerium",       140.116, 3.0, 58),
    ("CE4", "CE", "Cerium",       140.116, 4.0, 58),
    ("PR",  "PR", "Praseodymium", 140.908, 3.0, 59),
    ("PR3", "PR", "Praseodymium", 140.908, 3.0, 59),
    ("ND",  "ND", "Neodymium",    144.242, 3.0, 60),
    ("ND3", "ND", "Neodymium",    144.242, 3.0, 60),
    ("PM",  "PM", "Promethium",   145.000, 3.0, 61),
    ("SM",  "SM", "Samarium",     150.360, 3.0, 62),
    ("SM3", "SM", "Samarium",     150.360, 3.0, 62),
    ("EU",  "EU", "Europium",     151.964, 3.0, 63),
    ("EU3", "EU", "Europium",     151.964, 3.0, 63),
    ("GD",  "GD", "Gadolinium",   157.250, 3.0, 64),
    ("GD3", "GD", "Gadolinium",   157.250, 3.0, 64),
    ("TB",  "TB", "Terbium",      158.925, 3.0, 65),
    ("TB3", "TB", "Terbium",      158.925, 3.0, 65),
    ("DY",  "DY", "Dysprosium",   162.500, 3.0, 66),
    ("DY3", "DY", "Dysprosium",   162.500, 3.0, 66),
    ("HO",  "HO", "Holmium",      164.930, 3.0, 67),
    ("HO3", "HO", "Holmium",      164.930, 3.0, 67),
    ("ER",  "ER", "Erbium",       167.259, 3.0, 68),
    ("ER3", "ER", "Erbium",       167.259, 3.0, 68),
    ("TM",  "TM", "Thulium",      168.934, 3.0, 69),
    ("TM3", "TM", "Thulium",      168.934, 3.0, 69),
    ("YB",  "YB", "Ytterbium",    173.045, 3.0, 70),
    ("YB3", "YB", "Ytterbium",    173.045, 3.0, 70),
    ("LU",  "LU", "Lutetium",     174.967, 3.0, 71),
    ("LU3", "LU", "Lutetium",     174.967, 3.0, 71),
]
LN_INFO: Dict[str, dict] = {
    cid: {"sym": sym, "name": name, "mass": mass, "charge": chg, "atnum": atnum}
    for cid, sym, name, mass, chg, atnum in _LN_RAW
}
ALL_LN = set(LN_INFO.keys())

# LJ params: sym -> (Rmin/2 in Angstrom, |eps| in kcal/mol) - CHARMM convention
LN_LJ: Dict[str, Tuple[float, float]] = {
    "LA": (1.460, 0.1220), "CE": (1.447, 0.1200), "PR": (1.435, 0.1180),
    "ND": (1.424, 0.1160), "PM": (1.413, 0.1140), "SM": (1.403, 0.1120),
    "EU": (1.393, 0.1100), "GD": (1.384, 0.1085), "TB": (1.376, 0.1070),
    "DY": (1.368, 0.1055), "HO": (1.360, 0.1040), "ER": (1.353, 0.1025),
    "TM": (1.346, 0.1010), "YB": (1.340, 0.0995), "LU": (1.334, 0.0980),
}

def gmx_lj(sym: str) -> Tuple[float, float]:
    """Convert CHARMM Rmin/2 + eps(kcal/mol) to GROMACS sigma(nm) + eps(kJ/mol)."""
    r2, e_kcal = LN_LJ.get(sym, (1.38, 0.10))
    sigma_nm = 2.0 * r2 * (2.0 ** (-1.0 / 6.0)) / 10.0
    eps_kj   = e_kcal * 4.184
    return sigma_nm, eps_kj

STD_AA = {
    "ALA","ARG","ASN","ASP","CYS","CYX","GLN","GLU","GLY",
    "HSD","HSE","HSP","HIS","ILE","LEU","LYS","MET","PHE",
    "PRO","SER","THR","TRP","TYR","VAL","ACE","CT3","NME",
}
WATER_NAMES = {"HOH","WAT","H2O","SOL","TIP3","DOD","TIP"}
INORGANIC   = {"SOD","CLA","POT","CAL","ZN2","FE2","FE3","MG","MN","NI","CU"}

CHARMM_RESMAP: Dict[str, str] = {
    "HIS":"HSD","HID":"HSD","HIE":"HSE","HIP":"HSP",
    "HOH":"TIP3","WAT":"TIP3","H2O":"TIP3","DOD":"TIP3",
    "MSE":"MET","CSE":"CYS","HYP":"PRO",
    "NA":"SOD","K":"POT","CL":"CLA","CA":"CAL","ZN":"ZN2","FE":"FE2",
}
GROMACS_RESMAP: Dict[str, str] = {
    "HIS":"HIS","HID":"HIS","HIE":"HIS","HIP":"HIS",
    "HSD":"HIS","HSE":"HIS","HSP":"HIS",
    "HOH":"SOL","WAT":"SOL","H2O":"SOL","TIP3":"SOL","DOD":"SOL",
    "MSE":"MET","CSE":"CYS","HYP":"PRO",
    "NA":"NA+","K":"K+","CL":"CL-","CA":"CA2+","ZN":"ZN2+","FE":"FE2+",
}
ATOM_FIXES: Dict[str, Dict[str, str]] = {
    "*":   {"OT1":"O","OT2":"OXT"},
    "ILE": {"CD":"CD1"},
    "TIP3":{"O":"OH2","OW":"OH2"},
}

# ---------------------------------------------------------------------------
# SECTION 2 - HELPERS
# ---------------------------------------------------------------------------

class MotifSelect(Select):
    """Biopython PDBIO selector: keeps only residues in the motif set."""
    def __init__(self, model_id: int, motif_set: set):
        self.model_id  = model_id
        self.motif_set = motif_set
    def accept_model(self, m) -> int:
        return 1 if m.id == self.model_id else 0
    def accept_residue(self, r) -> int:
        return 1 if (r.get_parent().id, r.id) in self.motif_set else 0
    def accept_atom(self, a) -> int:
        return 1


def pdb_line(serial: int, name: str, resname: str, chain: str,
             resseq: int, x: float, y: float, z: float,
             occ: float = 1.0, bfac: float = 0.0,
             segid: str = "", element: str = "",
             hetatm: bool = False) -> str:
    rec   = "HETATM" if hetatm else "ATOM  "
    aname = (name + "   ")[:4] if len(name) >= 4 else (" " + name + "   ")[:4]
    return (
        "{rec}{serial:5d} {aname:4s} {resname:>3s} {chain:1s}{resseq:4d}    "
        "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}      "
        "{segid:<4s}  {element:>2s}\n"
    ).format(rec=rec, serial=serial, aname=aname, resname=resname,
             chain=chain, resseq=resseq, x=x, y=y, z=z,
             occ=occ, bfac=bfac, segid=segid, element=element)


def write_segment_pdb(residues: list, path: Path,
                      segid: str = "", force_hetatm: bool = False) -> None:
    lines  = []
    serial = 1
    for res in residues:
        rn   = res.resname.strip()
        ch   = res.get_parent().id
        rseq = res.id[1]
        het  = force_hetatm or (rn not in STD_AA and rn not in {"TIP3"})
        for atom in res.get_atoms():
            c  = atom.get_coord()
            el = (atom.element or atom.name[0]).strip()
            lines.append(pdb_line(
                serial, atom.name, rn, ch, rseq,
                float(c[0]), float(c[1]), float(c[2]),
                float(getattr(atom, "occupancy", 1.0)),
                float(getattr(atom, "bfactor",   0.0)),
                segid, el, het))
            serial += 1
    lines += ["TER\n", "END\n"]
    path.write_text("".join(lines), encoding="utf-8")

# ---------------------------------------------------------------------------
# SECTION 3 - DOWNLOAD
# ---------------------------------------------------------------------------

def download_pdb(pdb_id: str) -> Optional[Path]:
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    p = PDB_DIR / (pdb_id + ".pdb")
    if p.exists() and p.stat().st_size > 500:
        return p
    for attempt in range(3):
        try:
            r = requests.get(RCSB_URL.format(pdb_id), timeout=60)
            if r.status_code == 200:
                p.write_bytes(r.content)
                return p
            break
        except Exception:
            time.sleep(2 ** attempt)
    return None

# ---------------------------------------------------------------------------
# SECTION 4 - MOTIF EXTRACTION
# ---------------------------------------------------------------------------

def extract_motifs(pdb_id: str, pdb_path: Path) -> List[dict]:
    parser = PDBParser(QUIET=True)
    try:
        struct = parser.get_structure(pdb_id, str(pdb_path))
    except Exception as ex:
        print("      parse error: " + str(ex))
        return []

    records = []
    MOTIF_DIR.mkdir(parents=True, exist_ok=True)

    for model in struct:
        ns = NeighborSearch(list(model.get_atoms()))
        for chain in model:
            for res in chain:
                cid = res.resname.strip().upper()
                if cid not in ALL_LN:
                    continue
                ln_atoms = list(res.get_atoms())
                if not ln_atoms:
                    continue

                ln_coord = ln_atoms[0].coord

                shell_4 = ns.search(ln_coord, COORD_CUTOFF, level="R")
                water_3 = [r for r in ns.search(ln_coord, WATER_CUTOFF, level="R")
                           if r.resname.strip().upper() in WATER_NAMES]

                if not water_3:
                    continue  # no water within 3.0 A -> skip

                motif_set = {(chain.id, res.id)}
                for r in shell_4:
                    motif_set.add((r.get_parent().id, r.id))
                for r in water_3:
                    motif_set.add((r.get_parent().id, r.id))

                prot_shell  = [r for r in shell_4
                               if r.resname.strip().upper() not in ALL_LN
                               and r.resname.strip().upper() not in WATER_NAMES]
                water_shell = [r for r in shell_4
                               if r.resname.strip().upper() in WATER_NAMES]

                fname = "{}_{}_{:d}_{}.pdb".format(
                    pdb_id, chain.id, res.id[1], cid)
                io = PDBIO()
                io.set_structure(struct)
                io.save(str(MOTIF_DIR / fname),
                        MotifSelect(model.id, motif_set))

                records.append({
                    "pdb_id":         pdb_id,
                    "motif_file":     fname,
                    "comp_id":        cid,
                    "element":        LN_INFO[cid]["name"],
                    "chain":          chain.id,
                    "resnum":         res.id[1],
                    "coord_residues": len(prot_shell),
                    "waters_3A":      len(water_3),
                    "waters_4A":      len(water_shell),
                })
    return records

# ---------------------------------------------------------------------------
# SECTION 5 - STRUCTURE CLEANING
# ---------------------------------------------------------------------------

def clean_structure(struct, log: List[str]) -> list:
    """
    Apply all cleaning steps in-place.
    Returns list of (res1, res2) disulfide pairs detected.
    """

    # A: Remove non-A altlocs
    removed = 0
    for model in struct:
        for chain in model:
            for res in chain:
                drop = [a.id for a in res if a.altloc not in (" ", "A", "")]
                for aid in drop:
                    res.detach_child(aid)
                    removed += 1
                for atom in res:
                    atom.altloc = " "
    if removed:
        log.append("  [ALTLOC]   removed {:d} non-A altloc atoms".format(removed))

    # B: Fix insertion codes - renumber sequentially
    ic_count = 0
    for model in struct:
        for chain in model:
            cursor = 0
            for res in list(chain):
                hf, seq, ic = res.id
                if ic.strip():
                    res.id    = (hf, cursor + 1, " ")
                    ic_count += 1
                cursor = res.id[1]
    if ic_count:
        log.append("  [ICODE]    renumbered {:d} insertion-code residues".format(ic_count))

    # C: Fix residue names to CHARMM36 convention
    rn_count = 0
    for model in struct:
        for chain in model:
            for res in chain:
                old = res.resname.strip().upper()
                new = CHARMM_RESMAP.get(old)
                if new and new != old:
                    res.resname = new
                    rn_count  += 1
    if rn_count:
        log.append("  [RESNAME]  fixed {:d} residue names (CHARMM36)".format(rn_count))

    # D: Fix atom names
    an_count = 0
    for model in struct:
        for chain in model:
            for res in chain:
                rn  = res.resname.strip().upper()
                mp  = {}
                mp.update(ATOM_FIXES.get("*",  {}))
                mp.update(ATOM_FIXES.get(rn,   {}))
                for old_name, new_name in list(mp.items()):
                    if res.has_id(old_name):
                        atom          = res[old_name]
                        atom.name     = new_name
                        atom.fullname = (" " + new_name + "   ")[:4]
                        res.detach_child(old_name)
                        res.add(atom)
                        an_count += 1
    if an_count:
        log.append("  [ATOMNAME] fixed {:d} atom names".format(an_count))

    # E: Add OXT at C-terminus if missing
    oxt_count = 0
    for model in struct:
        for chain in model:
            prot = [r for r in chain if r.resname.strip() in STD_AA]
            if not prot:
                continue
            last = prot[-1]
            if (not last.has_id("OXT")
                    and last.has_id("O")
                    and last.has_id("C")):
                c_coord   = last["C"].coord
                o_coord   = last["O"].coord
                oxt_coord = 2.0 * c_coord - o_coord
                oxt = BioAtom("OXT", oxt_coord, 1.0, 0.0, " ", " OXT", None, "O")
                last.add(oxt)
                oxt_count += 1
    if oxt_count:
        log.append("  [OXT]      added OXT to {:d} C-termini".format(oxt_count))

    # F: Detect disulfides - label as CYX
    ssbonds = []
    for model in struct:
        cys_list = [r for ch in model for r in ch
                    if r.resname.strip() in ("CYS", "CYX") and r.has_id("SG")]
        for i, r1 in enumerate(cys_list):
            for r2 in cys_list[i + 1:]:
                if r1["SG"] - r2["SG"] < SS_CUTOFF:
                    r1.resname = "CYX"
                    r2.resname = "CYX"
                    ssbonds.append((r1, r2))
    if ssbonds:
        log.append("  [SS]       {:d} disulfide(s) detected -> CYX".format(len(ssbonds)))

    return ssbonds


def separate(struct) -> Tuple[dict, list, list, list]:
    """
    Split cleaned structure into groups:
      protein_chains : {chain_id: [residue, ...]}
      waters         : [residue, ...]
      ln_ions        : [residue, ...]
      other_ions     : [residue, ...]
    """
    protein_chains: Dict[str, list] = defaultdict(list)
    waters:     list = []
    ln_ions:    list = []
    other_ions: list = []

    for model in struct:
        for chain in model:
            for res in chain:
                rn = res.resname.strip().upper()
                if rn in ALL_LN:
                    ln_ions.append(res)
                elif rn == "TIP3":
                    waters.append(res)
                elif rn in INORGANIC:
                    other_ions.append(res)
                elif rn in STD_AA:
                    protein_chains[chain.id].append(res)

    return dict(protein_chains), waters, ln_ions, other_ions

# ---------------------------------------------------------------------------
# SECTION 6 - VMD / CHARMM36 OUTPUT
# ---------------------------------------------------------------------------

def write_rtf(ln_residues: list, path: Path) -> None:
    """CHARMM36 RTF: MASS + RESI blocks for every distinct Ln3+ present."""
    seen: Dict[str, dict] = {}
    for res in ln_residues:
        cid = res.resname.strip().upper()
        if cid not in seen:
            seen[cid] = LN_INFO[cid]

    lines = [
        "* Custom Ln3+ residue topology for CHARMM36\n",
        "* Generated by end_to_end_motif_pipeline.py\n",
        "*\n",
        "31 1\n\n",
    ]
    for idx, (cid, info) in enumerate(seen.items()):
        atype = info["sym"] + "3P"
        lines.append(
            "MASS {:4d}  {:<8s} {:9.4f}  {}  ! {}\n".format(
                200 + idx, atype, info["mass"], info["sym"], info["name"])
        )
    lines.append("\n")
    for cid, info in seen.items():
        atype = info["sym"] + "3P"
        lines += [
            "RESI {:<4s}  {:6.3f}   ! {} ion\n".format(cid, info["charge"], info["name"]),
            "GROUP\n",
            "ATOM {:<4s}  {:<8s}  {:6.3f}\n".format(cid, atype, info["charge"]),
            "PATCHING FIRS NONE LAST NONE\n\n",
        ]
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def write_prm(ln_residues: list, path: Path) -> None:
    """CHARMM36 PRM: NONBONDED LJ parameters for Ln3+ atom types."""
    syms = sorted({LN_INFO[r.resname.strip().upper()]["sym"] for r in ln_residues})
    lines = [
        "* Custom Ln3+ LJ parameters for CHARMM36\n",
        "* Rmin/2 (Angstrom)  |eps| (kcal/mol)\n",
        "* Generated by end_to_end_motif_pipeline.py\n",
        "*\n\n",
        "NONBONDED nbxmod 5 atom cdiel shift vatom vdistance vswitch -\n",
        "          cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0\n\n",
    ]
    for sym in syms:
        r2, eps = LN_LJ.get(sym, (1.38, 0.10))
        atype   = sym + "3P"
        lines.append(
            "{:<8s}  0.0  {:9.4f}  {:8.4f}   0.0  {:9.4f}  {:8.4f}  ! {}\n".format(
                atype, -eps, r2, -eps / 10.0, r2, sym)
        )
    lines.append("\nEND\n")
    path.write_text("".join(lines), encoding="utf-8")


def write_psfgen_tcl(stem: str, pchains: dict, has_water: bool,
                     ln_residues: list, has_other_ions: bool,
                     ssbonds: list, out_dir: Path) -> None:
    ln_comp_ids = sorted({r.resname.strip().upper() for r in ln_residues})
    lines = [
        "# psfgen.tcl - {}\n".format(stem),
        "# VMD CHARMM36 PSF/PDB generation\n",
        "# Run: vmd -dispdev text -e psfgen.tcl\n",
        "#\n",
        "# Place CHARMM36 topology files in ./charmm36/ before running.\n",
        "# Required: top_all36_prot.rtf, toppar/toppar_water_ions.str\n",
        "\npackage require psfgen\n\n",
        "set CHARMMDIR \"./charmm36\"\n",
        "topology $CHARMMDIR/top_all36_prot.rtf\n",
        "topology $CHARMMDIR/toppar/toppar_water_ions.str\n",
    ]
    if ln_residues:
        lines.append("topology ./ln_custom.rtf\n")

    lines.append("\n# --- Protein segments ---\n")
    for ch in pchains:
        seg = "PRO" + ch
        lines += [
            "segment {} {{\n".format(seg),
            "    pdb protein_{}_charmm.pdb\n".format(ch),
            "    first NTER\n",
            "    last  CTER\n",
            "}\n",
            "coordpdb protein_{}_charmm.pdb {}\n\n".format(ch, seg),
        ]

    for r1, r2 in ssbonds:
        c1 = r1.get_parent().id
        c2 = r2.get_parent().id
        lines.append("patch DISU PRO{}:{} PRO{}:{}\n".format(
            c1, r1.id[1], c2, r2.id[1]))

    if has_water:
        lines += [
            "\n# --- Water segment ---\n",
            "segment SOLV {\n",
            "    auto none\n",
            "    pdb water_charmm.pdb\n",
            "}\n",
            "coordpdb water_charmm.pdb SOLV\n",
        ]

    for idx, cid in enumerate(ln_comp_ids):
        seg = "LN{:02d}".format(idx + 1)
        lines += [
            "\n# --- Lanthanide {} ---\n".format(cid),
            "segment {} {{\n".format(seg),
            "    auto none\n",
            "    pdb ions_charmm.pdb\n",
            "}\n",
            "coordpdb ions_charmm.pdb {}\n".format(seg),
        ]

    if has_other_ions:
        lines += [
            "\n# --- Inorganic ions ---\n",
            "segment ION {\n",
            "    auto none\n",
            "    pdb ions_charmm.pdb\n",
            "}\n",
            "coordpdb ions_charmm.pdb ION\n",
        ]

    lines += [
        "\nguesscoord\n",
        "writepsf {}_vmd.psf\n".format(stem),
        "writepdb {}_vmd.pdb\n".format(stem),
        'puts "\\nDone: {}_vmd.psf and {}_vmd.pdb written."\n'.format(stem, stem),
    ]
    (out_dir / "psfgen.tcl").write_text("".join(lines), encoding="utf-8")


def write_vmd_bat(stem: str, out_dir: Path) -> None:
    content = (
        "@echo off\n"
        "cd /d \"%~dp0\"\n"
        "echo Running VMD psfgen for {}...\n"
        "vmd -dispdev text -e psfgen.tcl\n"
        "pause\n"
    ).format(stem)
    (out_dir / "run_vmd.bat").write_text(content, encoding="utf-8")

# ---------------------------------------------------------------------------
# SECTION 7 - GROMACS OUTPUT
# ---------------------------------------------------------------------------

def write_protein_gmx(pchains: dict, path: Path) -> None:
    """Protein-only PDB with GROMACS residue names (for gmx pdb2gmx)."""
    lines  = []
    serial = 1
    for ch_id, residues in pchains.items():
        for res in residues:
            rn_orig = res.resname.strip().upper()
            rn_gmx  = GROMACS_RESMAP.get(rn_orig, res.resname.strip())
            rseq    = res.id[1]
            for atom in res.get_atoms():
                c  = atom.get_coord()
                el = (atom.element or atom.name[0]).strip()
                lines.append(pdb_line(
                    serial, atom.name, rn_gmx, res.get_parent().id, rseq,
                    float(c[0]), float(c[1]), float(c[2]),
                    float(getattr(atom, "occupancy", 1.0)),
                    float(getattr(atom, "bfactor",   0.0)),
                    "", el, False))
                serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def write_ffnonbonded(ln_residues: list, path: Path) -> None:
    """GROMACS atomtypes section for Ln3+ (sigma nm, eps kJ/mol)."""
    syms = sorted({LN_INFO[r.resname.strip().upper()]["sym"] for r in ln_residues})
    lines = [
        "; Custom Ln3+ atom types for GROMACS CHARMM36\n",
        "; sigma (nm)  epsilon (kJ/mol)\n",
        "; Generated by end_to_end_motif_pipeline.py\n\n",
        "[ atomtypes ]\n",
        ";name       at.num   mass         charge  ptype  sigma         epsilon\n",
    ]
    for sym in syms:
        info     = next(v for k, v in LN_INFO.items() if v["sym"] == sym)
        sigma, e = gmx_lj(sym)
        atype    = sym + "3P"
        lines.append(
            "{:<10s}  {:4d}  {:10.4f}   {:6.3f}   A  {:12.8f}  {:12.6f}  ; {}\n".format(
                atype, info["atnum"], info["mass"],
                info["charge"], sigma, e, info["name"])
        )
    path.write_text("".join(lines), encoding="utf-8")


def write_ln_itp(ln_residues: list, out_dir: Path) -> List[str]:
    """One .itp per distinct Ln comp_id. Returns list of written filenames."""
    seen: Dict[str, dict] = {}
    for res in ln_residues:
        cid = res.resname.strip().upper()
        if cid not in seen:
            seen[cid] = LN_INFO[cid]
    written = []
    for cid, info in seen.items():
        atype = info["sym"] + "3P"
        fname = "ln_{}.itp".format(cid)
        lines = [
            "; Lanthanide topology: {} ({})\n".format(info["name"], cid),
            "; Generated by end_to_end_motif_pipeline.py\n\n",
            "[ moleculetype ]\n",
            "; Name            nrexcl\n",
            "{:<16s}  1\n\n".format(cid),
            "[ atoms ]\n",
            "; nr  type        resnr  residue  atom    cgnr  charge      mass\n",
            "  1   {:<10s}  1      {:<6s}  {:<6s}  1  {:8.4f}  {:10.4f}\n\n".format(
                atype, cid, cid, info["charge"], info["mass"]),
        ]
        (out_dir / fname).write_text("".join(lines), encoding="utf-8")
        written.append(fname)
    return written


def write_posre(pchains: dict, path: Path) -> None:
    """Position restraints ITP for all protein heavy atoms."""
    lines = [
        "; Position restraints on protein heavy atoms\n",
        "; Generated by end_to_end_motif_pipeline.py\n\n",
        "[ position_restraints ]\n",
        ";  i  funct   fcx    fcy    fcz\n",
    ]
    idx = 1
    for ch_id, residues in pchains.items():
        for res in residues:
            for atom in res.get_atoms():
                el = (atom.element or atom.name[0]).strip().upper()
                if el not in ("H", "D", ""):
                    lines.append("  {:5d}   1   1000   1000   1000\n".format(idx))
                idx += 1
    path.write_text("".join(lines), encoding="utf-8")


def write_topol(stem: str, pchains: dict, waters: list,
                ln_residues: list, other_ions: list,
                ln_itp_files: List[str], path: Path) -> None:
    """Master GROMACS topology file."""
    chains  = list(pchains.keys())
    ln_mols: Dict[str, int] = {}
    for res in ln_residues:
        cid = res.resname.strip().upper()
        ln_mols[cid] = ln_mols.get(cid, 0) + 1

    lines = [
        "; GROMACS topology: {}\n".format(stem),
        "; Force field: CHARMM36\n",
        "; Generated by end_to_end_motif_pipeline.py\n\n",
        '#include "charmm36.ff/forcefield.itp"\n',
        '#include "./ffnonbonded_ln.itp"\n\n',
        "; Protein topology (from gmx pdb2gmx - see run_gmx.bat)\n",
        '#include "./{}_protein.top"\n\n'.format(stem),
        "; Lanthanide molecule topologies\n",
    ]
    for fname in ln_itp_files:
        lines.append('#include "./{}"\n'.format(fname))

    lines += [
        "\n; Position restraints (active when -DPOSRES is passed to grompp)\n",
        "#ifdef POSRES\n",
        '#include "./posre.itp"\n',
        "#endif\n\n",
        "[ system ]\n",
        "{} lanthanide motif\n\n".format(stem),
        "[ molecules ]\n",
        "; Compound              #mols\n",
    ]
    for ch in chains:
        lines.append("Protein_{:<13s}  1\n".format(ch))
    for cid, n in ln_mols.items():
        lines.append("{:<22s}  {:d}\n".format(cid, n))
    if waters:
        lines.append("SOL                    {:d}\n".format(len(waters)))

    path.write_text("".join(lines), encoding="utf-8")


def write_mdp(path: Path, mdp_type: str = "em") -> None:
    if mdp_type == "em":
        text = (
            "; Energy minimisation parameters\n"
            "; Generated by end_to_end_motif_pipeline.py\n\n"
            "integrator    = steep\n"
            "emtol         = 1000.0\n"
            "emstep        = 0.01\n"
            "nsteps        = 50000\n\n"
            "nstlist       = 1\n"
            "cutoff-scheme = Verlet\n"
            "ns_type       = grid\n"
            "coulombtype   = PME\n"
            "rcoulomb      = 1.0\n"
            "rvdw          = 1.0\n"
            "pbc           = xyz\n"
        )
    else:
        text = (
            "; Production MD - 100 ps template\n"
            "; Generated by end_to_end_motif_pipeline.py\n\n"
            "integrator    = md\n"
            "dt            = 0.002\n"
            "nsteps        = 50000\n\n"
            "nstxout       = 500\n"
            "nstvout       = 500\n"
            "nstenergy     = 500\n"
            "nstlog        = 500\n\n"
            "cutoff-scheme = Verlet\n"
            "ns_type       = grid\n"
            "coulombtype   = PME\n"
            "rcoulomb      = 1.0\n"
            "rvdw          = 1.0\n\n"
            "tcoupl        = V-rescale\n"
            "tc-grps       = Protein  Non-Protein\n"
            "tau_t         = 0.1      0.1\n"
            "ref_t         = 300      300\n\n"
            "pcoupl        = Parrinello-Rahman\n"
            "pcoupltype    = isotropic\n"
            "tau_p         = 2.0\n"
            "ref_p         = 1.0\n"
            "compressibility = 4.5e-5\n\n"
            "pbc           = xyz\n"
            "gen_vel       = yes\n"
            "gen_temp      = 300\n"
        )
    path.write_text(text, encoding="utf-8")


def write_gmx_runners(stem: str, out_dir: Path) -> None:
    bat = (
        "@echo off\n"
        "rem GROMACS runner for {}\n"
        "cd /d \"%~dp0\"\n\n"
        "echo [1/4] Running pdb2gmx...\n"
        "gmx pdb2gmx -f protein_gmx.pdb -o {}_protein.gro -p {}_protein.top "
        "-i posre.itp -ff charmm36 -water tip3p -ter -ignh\n\n"
        "echo [2/4] Setting box...\n"
        "gmx editconf -f {}_protein.gro -o {}_box.gro -c -d 1.0 -bt cubic\n\n"
        "echo [3/4] Energy minimisation...\n"
        "gmx grompp -f em.mdp -c {}_box.gro -p topol.top -o em.tpr\n"
        "gmx mdrun  -v -deffnm em\n\n"
        "echo [4/4] Done. Check em.log for convergence.\n"
        "pause\n"
    ).format(stem, stem, stem, stem, stem, stem)

    sh = (
        "#!/bin/bash\n"
        "# GROMACS runner for {}\n"
        "# Linux/WSL: chmod +x run_gmx.sh && ./run_gmx.sh\n"
        "set -e\n"
        'cd "$(dirname \"$0\")"\n\n'
        "echo '[1/4] Running pdb2gmx...'\n"
        "gmx pdb2gmx -f protein_gmx.pdb -o {}_protein.gro -p {}_protein.top "
        "-i posre.itp -ff charmm36 -water tip3p -ter -ignh\n\n"
        "echo '[2/4] Setting box...'\n"
        "gmx editconf -f {}_protein.gro -o {}_box.gro -c -d 1.0 -bt cubic\n\n"
        "echo '[3/4] Energy minimisation...'\n"
        "gmx grompp -f em.mdp -c {}_box.gro -p topol.top -o em.tpr\n"
        "gmx mdrun -v -deffnm em\n\n"
        "echo '[4/4] Done.'\n"
    ).format(stem, stem, stem, stem, stem, stem)

    (out_dir / "run_gmx.bat").write_text(bat, encoding="utf-8")
    (out_dir / "run_gmx.sh").write_text(sh,  encoding="utf-8")

# ---------------------------------------------------------------------------
# SECTION 8 - PER-MOTIF ORCHESTRATION
# ---------------------------------------------------------------------------

def prep_motif(motif_path: Path, target: str) -> Path:
    stem    = motif_path.stem
    out_dir = motif_path.parent / (stem + "_prep")
    out_dir.mkdir(parents=True, exist_ok=True)

    log = [
        "=" * 60,
        "Motif  : " + motif_path.name,
        "Time   : " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Target : " + target,
        "=" * 60,
        "",
    ]

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure(stem, str(motif_path))
    log.append("Atoms loaded : {:d}".format(sum(1 for _ in struct.get_atoms())))

    ssbonds = clean_structure(struct, log)
    pchains, waters, ln_ions, other_ions = separate(struct)
    all_ions = ln_ions + other_ions

    log += [
        "",
        "Protein residues : {:d}  (chains: {})".format(
            sum(len(v) for v in pchains.values()), list(pchains.keys())),
        "Water (TIP3)     : {:d}".format(len(waters)),
        "Lanthanide ions  : {:d}  {}".format(
            len(ln_ions), [r.resname.strip() for r in ln_ions]),
        "Inorganic ions   : {:d}".format(len(other_ions)),
        "Disulfides       : {:d}".format(len(ssbonds)),
        "",
    ]

    if not waters:
        log.append("WARNING: No water molecules in motif.")

    # Shared segment PDBs (used by both VMD and GROMACS paths)
    for ch, residues in pchains.items():
        write_segment_pdb(residues,
                          out_dir / ("protein_" + ch + "_charmm.pdb"),
                          segid="PRO" + ch)
    if waters:
        write_segment_pdb(waters, out_dir / "water_charmm.pdb", segid="SOLV")
    if all_ions:
        write_segment_pdb(all_ions, out_dir / "ions_charmm.pdb",
                          segid="ION", force_hetatm=True)

    # Full merged PDB for inspection
    all_res = [r for rs in pchains.values() for r in rs]
    all_res += ln_ions + waters + other_ions
    write_segment_pdb(all_res, out_dir / "clean_full.pdb", segid="ALL")

    # Stage 3a: VMD / CHARMM36
    if target in ("vmd", "both"):
        log.append("-- VMD / CHARMM36 files --")
        if ln_ions:
            write_rtf(ln_ions, out_dir / "ln_custom.rtf")
            write_prm(ln_ions, out_dir / "ln_custom.prm")
            log.append("  ln_custom.rtf  ln_custom.prm")
        write_psfgen_tcl(stem, pchains, bool(waters),
                         ln_ions, bool(other_ions), ssbonds, out_dir)
        write_vmd_bat(stem, out_dir)
        log.append("  psfgen.tcl  run_vmd.bat")
        for ch in pchains:
            log.append("  protein_{}_charmm.pdb".format(ch))
        if waters:
            log.append("  water_charmm.pdb")
        if all_ions:
            log.append("  ions_charmm.pdb")

    # Stage 3b: GROMACS
    if target in ("gromacs", "both"):
        log.append("-- GROMACS / CHARMM36 files --")
        write_protein_gmx(pchains, out_dir / "protein_gmx.pdb")
        log.append("  protein_gmx.pdb")
        if ln_ions:
            write_ffnonbonded(ln_ions, out_dir / "ffnonbonded_ln.itp")
            ln_itps = write_ln_itp(ln_ions, out_dir)
            log.append("  ffnonbonded_ln.itp  " + "  ".join(ln_itps))
        else:
            ln_itps = []
        write_posre(pchains, out_dir / "posre.itp")
        write_topol(stem, pchains, waters, ln_ions, other_ions,
                    ln_itps, out_dir / "topol.top")
        write_mdp(out_dir / "em.mdp", "em")
        write_mdp(out_dir / "md.mdp", "md")
        write_gmx_runners(stem, out_dir)
        log.append("  posre.itp  topol.top  em.mdp  md.mdp")
        log.append("  run_gmx.bat  run_gmx.sh")

    log += [
        "",
        "-- NEXT STEPS --",
        "VMD:     cd motifs/{}_prep && run_vmd.bat".format(stem),
        "         (place CHARMM36 RTF files in ./charmm36/ first)",
        "GROMACS: cd motifs/{}_prep && run_gmx.bat  (Windows)".format(stem),
        "         OR  ./run_gmx.sh  (Linux/WSL)",
        "         charmm36.ff must be in GROMACS data path or current folder",
    ]

    (out_dir / "prep.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    return out_dir

# ---------------------------------------------------------------------------
# SECTION 9 - MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Lanthanide motif pipeline: VMD + GROMACS output"
    )
    ap.add_argument("--ids", nargs="*",
                    help="PDB IDs to process (default: read lanthanide_pdb_ids.txt)")
    ap.add_argument("--target", choices=["vmd", "gromacs", "both"],
                    default="both",
                    help="Output target  (default: both)")
    ap.add_argument("--skip-extract", action="store_true",
                    help="Skip download+extraction; prep existing motifs/ only")
    args = ap.parse_args()

    print("")
    print("=" * 65)
    print("  HackRush 26 | Lanthanide Motif Pipeline")
    print("  Target : " + args.target.upper())
    print("=" * 65)

    PDB_DIR.mkdir(parents=True, exist_ok=True)
    MOTIF_DIR.mkdir(parents=True, exist_ok=True)

    all_records: List[dict] = []

    # Stage 1: Download and extract
    if not args.skip_extract:
        ids: List[str] = list(args.ids or [])
        if not ids and IDS_FILE.exists():
            raw = IDS_FILE.read_text(encoding="utf-8").splitlines()
            ids = [line.strip().upper() for line in raw
                   if line.strip() and not line.startswith("#")]
        if not ids:
            print("  No IDs found. Use --ids or populate lanthanide_pdb_ids.txt")
        else:
            print("")
            print("[STAGE 1] Extracting motifs from {:d} structures".format(len(ids)))
            print("-" * 65)
            for i, pdb_id in enumerate(ids, 1):
                print("  [{:4d}/{}] {}".format(i, len(ids), pdb_id),
                      end="  ", flush=True)
                pdb_path = download_pdb(pdb_id)
                if pdb_path is None:
                    print("[ERR] download failed")
                    continue
                motifs = extract_motifs(pdb_id, pdb_path)
                if motifs:
                    print("[OK] {:d} motif(s)".format(len(motifs)))
                else:
                    print("-- no valid motifs")
                all_records.extend(motifs)
                time.sleep(0.05)

            if all_records:
                with open(MOTIF_SUMMARY, "w", newline="", encoding="utf-8") as fh:
                    w = csv.DictWriter(fh, fieldnames=list(all_records[0].keys()))
                    w.writeheader()
                    w.writerows(all_records)
                print("")
                print("  motif_summary.csv -> {:d} motifs".format(len(all_records)))
    else:
        print("[STAGE 1] Skipped (--skip-extract)")

    # Stage 2+3: Clean and generate output
    motif_files = sorted(MOTIF_DIR.glob("*.pdb"))
    print("")
    print("[STAGE 2+3] Cleaning + generating {} output".format(
        args.target.upper()))
    print("  {:d} motif PDB(s) found in motifs/".format(len(motif_files)))
    print("-" * 65)

    for j, mp in enumerate(motif_files, 1):
        print("  [{:3d}/{}] {}".format(j, len(motif_files), mp.name))
        od = prep_motif(mp, target=args.target)
        print("         -> " + od.name + "/")

    print("")
    print("=" * 65)
    print("  Done: {:d} motif(s) prepared".format(len(motif_files)))
    print("-" * 65)
    if args.target in ("vmd", "both"):
        print("  VMD:     cd motifs/<motif>_prep && run_vmd.bat")
        print("           (place CHARMM36 RTF in ./charmm36/ first)")
    if args.target in ("gromacs", "both"):
        print("  GROMACS: cd motifs/<motif>_prep && run_gmx.bat  (Windows)")
        print("           OR  ./run_gmx.sh  (Linux/WSL)")
    print("=" * 65)
    print("")


if __name__ == "__main__":
    main()
