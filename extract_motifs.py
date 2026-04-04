#!/usr/bin/env python3
"""
HackRush '26 | Motif Coordinate Extractor
==========================================
Reads: lanthanide_pdb_ids.txt  (from Checkpoint 1)
Does:
  1. Downloads each PDB file from RCSB
  2. Parses with Biopython PDBParser
  3. For each Ln³⁺ ion in the structure:
       a. NeighborSearch → residues within 4.0 Å  (coordinating shell)
       b. NeighborSearch → HOH within 3.0 Å       (water filter — REQUIRED)
       c. If HOH found → extract & save motif PDB
Writes:
  pdb_files/          raw downloaded PDB files (cached)
  motifs/             one PDB file per valid lanthanide binding motif
  motif_summary.csv   summary table of all motifs

Naming convention:
  {PDB_ID}_{CHAIN}_{RES_NUM}_{COMP_ID}.pdb
  e.g.  1H87_A_501_GD3.pdb

Install: pip install biopython requests
Run:     python extract_motifs.py
"""

import os, csv, time, warnings
import requests
from pathlib import Path
from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# ─────────────────────────────────────────────────────────────
# CONFIG — thresholds from problem statement
# ─────────────────────────────────────────────────────────────
COORD_CUTOFF = 4.0   # Å  coordinating residues
WATER_CUTOFF = 3.0   # Å  water molecule (HOH) — REQUIRED for valid motif

SCRIPT_DIR  = Path(os.path.dirname(os.path.abspath(__file__)))
IDS_FILE    = SCRIPT_DIR / "lanthanide_pdb_ids.txt"
PDB_DIR     = SCRIPT_DIR / "pdb_files"
MOTIF_DIR   = SCRIPT_DIR / "motifs"
SUMMARY_CSV = SCRIPT_DIR / "motif_summary.csv"

RCSB_PDB_URL = "https://files.rcsb.org/download/{}.pdb"

ALL_LN_COMP_IDS = {
    "LA","LA3","CE","CE3","CE4","PR","PR3","ND","ND3","PM",
    "SM","SM3","EU","EU3","GD","GD3","TB","TB3","DY","DY3",
    "HO","HO3","ER","ER3","TM","TM3","YB","YB3","LU","LU3",
}
LN_ELEMENT_NAME = {
    "LA":"Lanthanum","LA3":"Lanthanum","CE":"Cerium","CE3":"Cerium","CE4":"Cerium",
    "PR":"Praseodymium","PR3":"Praseodymium","ND":"Neodymium","ND3":"Neodymium",
    "PM":"Promethium","SM":"Samarium","SM3":"Samarium","EU":"Europium","EU3":"Europium",
    "GD":"Gadolinium","GD3":"Gadolinium","TB":"Terbium","TB3":"Terbium",
    "DY":"Dysprosium","DY3":"Dysprosium","HO":"Holmium","HO3":"Holmium",
    "ER":"Erbium","ER3":"Erbium","TM":"Thulium","TM3":"Thulium",
    "YB":"Ytterbium","YB3":"Ytterbium","LU":"Lutetium","LU3":"Lutetium",
}


# ─────────────────────────────────────────────────────────────
# SELECT: which residues go into the motif PDB
# ─────────────────────────────────────────────────────────────

class MotifSelect(Select):
    """
    Biopython PDBIO Select subclass.
    Accepts only residues in the motif set and the first model.
    """
    def __init__(self, model_id: int, motif_set: set):
        self.model_id  = model_id
        self.motif_set = motif_set  # set of (chain_id, residue.id_tuple)

    def accept_model(self, model):
        return 1 if model.id == self.model_id else 0

    def accept_chain(self, chain):
        return 1  # chain filtering done at residue level

    def accept_residue(self, residue):
        key = (residue.get_parent().id, residue.id)
        return 1 if key in self.motif_set else 0

    def accept_atom(self, atom):
        return 1  # keep all atoms of accepted residues


# ─────────────────────────────────────────────────────────────
# DOWNLOAD
# ─────────────────────────────────────────────────────────────

def download_pdb(pdb_id: str) -> Path | None:
    path = PDB_DIR / f"{pdb_id}.pdb"
    if path.exists() and path.stat().st_size > 100:
        return path  # use cached file
    try:
        r = requests.get(RCSB_PDB_URL.format(pdb_id), timeout=30)
        if r.status_code == 200:
            path.write_bytes(r.content)
            return path
        print(f" HTTP {r.status_code}", end="")
    except Exception as e:
        print(f" download error: {e}", end="")
    return None


# ─────────────────────────────────────────────────────────────
# EXTRACT MOTIFS  (core logic)
# ─────────────────────────────────────────────────────────────

def extract_motifs(pdb_id: str, pdb_path: Path) -> list:
    """
    Parses PDB with Bio.PDB.PDBParser.
    Uses NeighborSearch to find coordinating shell around each Ln³⁺ ion.

    Motif is valid when:
      ✓ At least 1 HOH within WATER_CUTOFF (3.0 Å)

    Motif contains:
      • The Ln³⁺ residue itself
      • All residues (protein + HOH) within COORD_CUTOFF (4.0 Å)

    Each valid motif is saved as an individual .pdb file.
    Returns list of summary dicts.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, str(pdb_path))
    except Exception as e:
        print(f" parse error: {e}", end="")
        return []

    records = []

    for model in structure:
        # Build NeighborSearch over all atoms in this model
        all_atoms = list(model.get_atoms())
        if not all_atoms:
            continue
        ns = NeighborSearch(all_atoms)

        for chain in model:
            for residue in chain:
                comp_id = residue.resname.strip().upper()
                if comp_id not in ALL_LN_COMP_IDS:
                    continue

                # ── Lanthanide residue found ─────────────────────
                ln_atoms = list(residue.get_atoms())
                if not ln_atoms:
                    continue
                ln_coord = ln_atoms[0].coord   # Ln³⁺ is a single atom

                # ── Step A: residues within 4.0 Å ───────────────
                coord_residues = ns.search(ln_coord, COORD_CUTOFF, level="R")

                # ── Step B: HOH within 3.0 Å (required filter) ──
                water_nearby = [
                    r for r in ns.search(ln_coord, WATER_CUTOFF, level="R")
                    if r.resname.strip().upper() == "HOH"
                ]

                if not water_nearby:
                    continue  # no coordinating water → skip

                # ── Build motif residue set ──────────────────────
                motif_set = set()
                motif_set.add((chain.id, residue.id))      # the Ln³⁺ itself
                for r in coord_residues:
                    motif_set.add((r.get_parent().id, r.id))
                for r in water_nearby:
                    motif_set.add((r.get_parent().id, r.id))  # ensure HOH included

                # ── Categorise coordinating residues ────────────
                protein_shell = [
                    r for r in coord_residues
                    if r.resname.strip().upper() not in ALL_LN_COMP_IDS
                    and r.resname.strip().upper() != "HOH"
                ]
                water_shell = [
                    r for r in coord_residues
                    if r.resname.strip().upper() == "HOH"
                ]

                # ── Save motif PDB ───────────────────────────────
                fname = f"{pdb_id}_{chain.id}_{residue.id[1]}_{comp_id}.pdb"
                out   = MOTIF_DIR / fname
                io    = PDBIO()
                io.set_structure(structure)
                io.save(str(out), MotifSelect(model.id, motif_set))

                records.append({
                    "pdb_id":              pdb_id,
                    "motif_file":          fname,
                    "ln_element":          LN_ELEMENT_NAME.get(comp_id, comp_id),
                    "comp_id":             comp_id,
                    "chain":               chain.id,
                    "ln_residue_num":      residue.id[1],
                    "ln_x":                round(float(ln_coord[0]), 3),
                    "ln_y":                round(float(ln_coord[1]), 3),
                    "ln_z":                round(float(ln_coord[2]), 3),
                    "coord_residues_4A":   len(protein_shell),
                    "water_mols_3A":       len(water_nearby),
                    "total_waters_4A":     len(water_shell),
                    "coordinating_aa":     "|".join(
                        f"{r.resname}{r.id[1]}" for r in protein_shell
                    ),
                })

    return records


# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────

def load_ids() -> list:
    if not IDS_FILE.exists():
        print(f"  ! {IDS_FILE} not found. Run checkpoint1 first.")
        return []
    ids = []
    with open(IDS_FILE) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line.upper())
    print(f"  Loaded {len(ids)} PDB IDs from {IDS_FILE.name}")
    return ids


def main():
    print("\n" + "="*62)
    print("  HackRush '26 | Lanthanide Binding Motif Extractor")
    print(f"  Coord cutoff : {COORD_CUTOFF} Å  |  Water cutoff : {WATER_CUTOFF} Å")
    print("="*62)

    PDB_DIR.mkdir(exist_ok=True)
    MOTIF_DIR.mkdir(exist_ok=True)

    pdb_ids = load_ids()
    if not pdb_ids:
        return

    all_motifs = []
    no_motif   = []
    failed     = []

    for idx, pdb_id in enumerate(pdb_ids, 1):
        print(f"  [{idx:>4}/{len(pdb_ids)}]  {pdb_id:<6}", end="  ")

        # Step 1: download PDB
        pdb_path = download_pdb(pdb_id)
        if not pdb_path:
            print("✗ download failed")
            failed.append(pdb_id)
            continue

        # Step 2: extract motifs
        motifs = extract_motifs(pdb_id, pdb_path)

        if motifs:
            print(f"✓  {len(motifs)} motif(s)  "
                  f"[{', '.join(m['comp_id'] for m in motifs[:3])}]")
            all_motifs.extend(motifs)
        else:
            print("–  no valid motif (no HOH ≤3Å)")
            no_motif.append(pdb_id)

        time.sleep(0.05)   # be polite to RCSB server

    # ── Save summary CSV ────────────────────────────────────
    fieldnames = [
        "pdb_id","motif_file","ln_element","comp_id","chain",
        "ln_residue_num","ln_x","ln_y","ln_z",
        "coord_residues_4A","water_mols_3A","total_waters_4A","coordinating_aa"
    ]
    with open(SUMMARY_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_motifs)

    # ── Final report ────────────────────────────────────────
    print("\n" + "="*62)
    print(f"  PDB structures processed   : {len(pdb_ids)}")
    print(f"  Structures with motif(s)   : {len(pdb_ids)-len(no_motif)-len(failed)}")
    print(f"  Valid motifs extracted     : {len(all_motifs)}")
    print(f"  No valid motif (no HOH)    : {len(no_motif)}")
    print(f"  Failed downloads           : {len(failed)}")
    print(f"\n  Motif PDB files → {MOTIF_DIR}/")
    print(f"  Summary CSV     → {SUMMARY_CSV}")
    print("="*62)

    if all_motifs:
        print(f"\n  Sample motifs:")
        print(f"  {'PDB':<6} {'Comp':<6} {'Ch':<3} {'Res#':<6} "
              f"{'Coord':<6} {'HOH':<5} Coordinating AAs")
        print("  " + "-"*65)
        for m in all_motifs[:12]:
            print(f"  {m['pdb_id']:<6} {m['comp_id']:<6} {m['chain']:<3} "
                  f"{m['ln_residue_num']:<6} {m['coord_residues_4A']:<6} "
                  f"{m['water_mols_3A']:<5} {m['coordinating_aa'][:30]}")


if __name__ == "__main__":
    main()
