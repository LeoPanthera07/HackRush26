# 🧬 HackRush '26 — Problem 9B
## Extraction of Lanthanide-Binding Motifs from Protein Structures

> **Team:** Rick and Morty
> 
> **Members:** Avishka Jindal · Mihir Gor
> 
> **Event:** HackRush 2026, IIT Gandhinagar
> 
> **Stakeholders:** Aditya Mehta · Adwaith P

---

<img width="706" height="698" alt="image" src="https://github.com/user-attachments/assets/1bdaa31f-252e-4042-bf5d-fcae3ceb6def" />

Lanthanide Binding Motif — 1GGY YB site

*YB³⁺ binding motif extracted from PDB: 1GGY, Chain A, Residue 732 — visualized in VMD 2.0.0a6*

---

## 📌 Problem Overview

Lanthanide ions (La³⁺ → Lu³⁺) are widely used in structural biology, bioimaging, and catalysis due to their unique coordination chemistry. This project builds a **fully automated pipeline** that:

1. Retrieves and verifies PDB structures containing confirmed Ln³⁺ ions
2. Extracts lanthanide-binding motifs (coordination shell + water)
3. Cleans and preprocesses structures for molecular dynamics
4. Generates simulation-ready files for **VMD/CHARMM36** and **GROMACS**
5. Handles errors in `pdb2gmx` with minimal manual intervention

---

## 📁 Repository Structure

```
HackRush26/
│
├── checkpoint1_retrieve_pdb_ids.py   # Stage 1: RCSB search + verification
├── lanthanide_pdb_ids.txt            # Verified PDB IDs with confirmed Ln³⁺
├── motif_locations.csv               # Chain / residue / comp_id for every Ln³⁺ site
│
├── output/
│   ├── 1GGY_A_732_YB_prep/          # Prepared simulation files for 1GGY YB site
│   │   ├── protein_gmx.pdb          # GROMACS-ready protein
│   │   ├── ffnonbonded_ln.itp       # Custom Ln³⁺ LJ parameters
│   │   ├── topol.top                # Master GROMACS topology
│   │   ├── em.mdp / md.mdp          # Energy minimisation + MD run params
│   │   ├── psfgen.tcl               # VMD psfgen script
│   │   ├── ln_custom.rtf/.prm       # CHARMM36 RTF/PRM for YB³⁺
│   │   └── prep.log                 # Audit log of all automated changes
│   └── ...
│
├── 1GGY_A_732_YB_prep.zip            # Packaged checkpoint files
├── RandM_HackRush26_9B_Report.docx   # Full technical report
└── README.md
```

---

## 🔬 Pipeline Overview

### Checkpoint 1 — PDB ID Retrieval (`checkpoint1_retrieve_pdb_ids.py`)

A **zero-false-positive** three-stage extractor:

| Stage | Method | Purpose |
|-------|--------|---------|
| **Search** | RCSB full-text API | Collect raw candidates for all 15 Ln element names + comp IDs |
| **Verify** | RCSB DataQuery on `pdbx_entity_nonpoly.comp_id` | Keep only entries with actual Ln³⁺ present |
| **Locate** | `DataSchema.find_field_names()` (runtime discovery) | Find chain ID + residue number fields without hardcoding |

Covers all 15 lanthanides: `La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu`

---

### Full Pipeline — Motif Extraction & Prep (`end_to_end_motif_pipeline.py`)

```
PDB Download → Motif Extraction → Structure Cleaning → VMD + GROMACS Output
```

**Motif criteria:**
- All residues within **4.0 Å** of the Ln³⁺ ion (coordination shell)
- At least **1 water molecule within 3.0 Å** (required for valid motif)

**Cleaning steps:**
- Remove non-A alternate conformations
- Fix insertion codes (renumber sequentially)
- CHARMM36 residue renaming: `HIS→HSD`, `HOH→TIP3`, `MSE→MET` …
- Atom name fixes: `ILE CD→CD1`, `TIP3 O→OH2`, `OT1→O`, `OT2→OXT`
- Add missing OXT at C-terminus
- Detect disulfides (SG–SG < 2.5 Å) → label CYX

---

## ⚠️ Error-Handling Framework

A key contribution of this project is an automated error-handling system for `gmx pdb2gmx`:

| Error | Detection | Automated Fix |
|-------|-----------|---------------|
| Unknown HETATM / residue | Residue not in CHARMM36 topology | Rename via alias map; Ln³⁺ written as standalone ITP |
| Alternate conformations | altLoc ≠ `' '` or `'A'` | Keep only altLoc A |
| Insertion codes | `res.id[2] ≠ ' '` | Renumber sequentially |
| Missing OXT terminus | Check last residue per chain | Compute & add OXT from C/CA positions |
| Non-standard HIS | HID / HIE / HIP / HSD | Map all → HIS for GROMACS |
| Disulfide bonds | SG–SG < 2.5 Å | Rename CYS→CYX; add SSBOND record |
| Water name mismatch | HOH / WAT / TIP3 / DOD | Rename all → SOL (GROMACS) |
| Missing Ln³⁺ FF params | Not in any standard FF | Generate `ffnonbonded_ln.itp` from tabulated Rmin/2 + ε |
| Chain fragments | < 3 residues or missing N/C | Flag in `prep.log`; skip fragment; continue |

### LJ Parameter Conversion (CHARMM → GROMACS)

```
σ (nm)      = 2 × Rmin/2 × 2^(−1/6) / 10
ε (kJ/mol)  = ε (kcal/mol) × 4.184
```

All 15 lanthanides have tabulated parameters built into the pipeline.

---

## 🚀 Usage

### Install dependencies
```bash
pip install requests rcsb-api biopython
```

### Checkpoint 1 — Retrieve PDB IDs
```bash
python checkpoint1_retrieve_pdb_ids.py
# → lanthanide_pdb_ids.txt
# → motif_locations.csv
```

### Full Pipeline
```bash
# All IDs from lanthanide_pdb_ids.txt
python end_to_end_motif_pipeline.py

# Specific PDB IDs
python end_to_end_motif_pipeline.py --ids 1GGY 2XU3

# GROMACS output only
python end_to_end_motif_pipeline.py --target gromacs

# VMD output only
python end_to_end_motif_pipeline.py --target vmd
```

### GROMACS Simulation
```bash
cd output/1GGY_A_732_YB_prep

gmx pdb2gmx -f protein_gmx.pdb -o protein_proc.gro \
            -water tip3p -ff charmm36 -ter -ignh
```
> Force field: **CHARMM36** | Water model: **TIP3P**
> Terminals: **NH₃⁺** (N-term) · **COO⁻** (C-term) via `-ter` flag

---

## 📊 Results — Checkpoint 1

Output files from the verified RCSB search are included in this repo:

- **`lanthanide_pdb_ids.txt`** — Verified PDB IDs (zero false positives)
- **`motif_locations.csv`** — Per-site: PDB ID, comp_id, chain, residue number, chemical name

---

## 📄 Report

Full technical report with methodology, error-handling tables, and design decisions:
**`RandM_HackRush26_9B_Report.docx`**

---

## 🛠 Dependencies

| Package | Purpose |
|---------|---------|
| `biopython` | PDB parsing, NeighborSearch, structure I/O |
| `requests` | RCSB PDB file download |
| `rcsb-api` | RCSB Search + Data API (DataQuery, DataSchema) |
| VMD 2.0+ | PSF generation, visualization |
| GROMACS | Topology + simulation |

---

*HackRush 2026 · IIT Gandhinagar · Team Rick and Morty*
