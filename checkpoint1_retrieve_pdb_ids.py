#!/usr/bin/env python3
"""
HackRush '26 | Checkpoint 1
Zero-false-positive Lanthanide PDB ID Extractor

STRATEGY:
  Stage 1 — SEARCH  : full_text on comp IDs → raw candidates
  Stage 2 — VERIFY  : DataQuery on pdbx_entity_nonpoly.comp_id ONLY
                       (minimal query, guaranteed valid field)
                       → keep ONLY entries with actual Ln3+ comp_id
  Stage 3 — LOCATE  : DataSchema.find_field_names() discovers correct
                       chain/residue fields at runtime (no hardcoding)
                       → motif_locations.csv

Install: pip install requests rcsb-api
Run:     python checkpoint1_retrieve_pdb_ids.py
"""

import os, csv, time
import requests
from pathlib import Path

SCRIPT_DIR  = Path(os.path.dirname(os.path.abspath(__file__)))
IDS_FILE    = SCRIPT_DIR / "lanthanide_pdb_ids.txt"
CSV_FILE    = SCRIPT_DIR / "motif_locations.csv"

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
HEADERS     = {"Content-Type": "application/json"}
BATCH_SIZE  = 50

LANTHANIDES = {
    "La": ("lanthanum",    ["LA",  "LA3"        ]),
    "Ce": ("cerium",       ["CE",  "CE3",  "CE4"]),
    "Pr": ("praseodymium", ["PR",  "PR3"        ]),
    "Nd": ("neodymium",    ["ND",  "ND3"        ]),
    "Pm": ("promethium",   ["PM"               ]),
    "Sm": ("samarium",     ["SM",  "SM3"        ]),
    "Eu": ("europium",     ["EU",  "EU3"        ]),
    "Gd": ("gadolinium",   ["GD",  "GD3"        ]),
    "Tb": ("terbium",      ["TB",  "TB3"        ]),
    "Dy": ("dysprosium",   ["DY",  "DY3"        ]),
    "Ho": ("holmium",      ["HO",  "HO3"        ]),
    "Er": ("erbium",       ["ER",  "ER3"        ]),
    "Tm": ("thulium",      ["TM",  "TM3"        ]),
    "Yb": ("ytterbium",    ["YB",  "YB3"        ]),
    "Lu": ("lutetium",     ["LU",  "LU3"        ]),
}
ALL_LN_COMP_IDS = {c for _, (_, cids) in LANTHANIDES.items() for c in cids}


# ─────────────────────────────────────────────────────────────
# STAGE 1 — SEARCH (raw HTTP)
# ─────────────────────────────────────────────────────────────

def _parse(rs: list) -> list:
    out = []
    for x in rs:
        v = x if isinstance(x, str) else (x.get("identifier") or x.get("rcsb_id") or "")
        if v:
            out.append(v.strip().upper())
    return out

def fulltext(term: str) -> list:
    try:
        r = requests.post(RCSB_SEARCH, headers=HEADERS, timeout=60, json={
            "query": {"type": "terminal", "service": "full_text",
                      "parameters": {"value": f'"{term}"'}},
            "return_type": "entry",
            "request_options": {"return_all_hits": True}
        })
        return _parse(r.json().get("result_set", [])) if r.status_code == 200 else []
    except Exception:
        return []

def get_candidates() -> list:
    all_ids = set()
    print(f"  {'Sym':<4} {'Element':<20} Candidates")
    print("  " + "-"*38)
    for sym, (name, comp_ids) in LANTHANIDES.items():
        ids = set()
        for term in [name] + comp_ids:
            ids.update(fulltext(term))
            time.sleep(0.15)
        all_ids.update(ids)
        print(f"  [{sym:<2}] {name:<20} {len(ids):>6}")
    return sorted(all_ids)


# ─────────────────────────────────────────────────────────────
# STAGE 2 — VERIFY  (DataQuery with ONLY comp_id — no bad fields)
# ─────────────────────────────────────────────────────────────

def verify(candidates: list) -> tuple:
    """
    Uses the MINIMAL valid DataQuery (rcbs.pdf):
      DataQuery(input_type="entry",
                input_ids=[...],
                return_data_list=["rcsb_id",
                  "nonpolymer_entities.pdbx_entity_nonpoly.comp_id"])
    Keeps only entries where comp_id ∈ ALL_LN_COMP_IDS.
    Returns (verified_ids, {pdb_id: [comp_ids found]})
    """
    from rcsbapi.data import DataQuery

    verified     = []
    comp_id_map  = {}   # pdb_id → list of Ln comp_ids found
    total        = len(candidates)

    for i in range(0, total, BATCH_SIZE):
        batch = candidates[i : i + BATCH_SIZE]
        end   = min(i + BATCH_SIZE, total)
        print(f"  [{i+1:>5}–{end:>5} / {total}]", end="  ", flush=True)

        try:
            result = DataQuery(
                input_type="entry",
                input_ids=batch,
                # ONLY comp_id — the one field guaranteed to be valid
                return_data_list=[
                    "rcsb_id",
                    "nonpolymer_entities.pdbx_entity_nonpoly.comp_id",
                ],
                suppress_autocomplete_warning=True
            ).exec()
        except Exception as e:
            print(f"DataQuery error: {e}")
            time.sleep(1)
            continue

        entries = (result.get("data") or {}).get("entries") or []
        passed  = 0

        for entry in entries:
            pdb_id = (entry.get("rcsb_id") or "").upper()
            npe    = entry.get("nonpolymer_entities") or []

            found_ln = []
            for ent in npe:
                cid = ((ent.get("pdbx_entity_nonpoly") or {}).get("comp_id") or "").upper()
                if cid in ALL_LN_COMP_IDS:
                    found_ln.append(cid)

            if found_ln:
                verified.append(pdb_id)
                comp_id_map[pdb_id] = found_ln
                passed += 1

        print(f"passed  {passed} / {len(batch)}")
        time.sleep(0.2)

    return sorted(set(verified)), comp_id_map


# ─────────────────────────────────────────────────────────────
# STAGE 3 — LOCATE  (DataSchema discovers correct field names)
# ─────────────────────────────────────────────────────────────

def discover_location_fields() -> tuple:
    """
    Uses DataSchema.find_field_names() from rcbs.pdf to discover
    the correct chain_id and seq_num field paths at runtime.
    Never hardcodes field names that might be wrong.
    """
    from rcsbapi.data import DataSchema
    schema = DataSchema()

    # Find chain ID field under nonpolymer entities
    chain_field = None
    for search_term in ["auth_asym_id", "asym_id"]:
        try:
            matches = schema.find_field_names(search_term)
            for m in matches:
                if "nonpolymer" in m.lower() and "container_identifiers" in m.lower():
                    chain_field = m
                    break
            if chain_field:
                break
        except ValueError:
            pass

    # Find sequence/residue number field under nonpolymer entities
    seq_field = None
    for search_term in ["auth_seq_id", "seq_num", "seq_id"]:
        try:
            matches = schema.find_field_names(search_term)
            for m in matches:
                if "nonpolymer" in m.lower() and "container_identifiers" in m.lower():
                    seq_field = m
                    break
            if seq_field:
                break
        except ValueError:
            pass

    print(f"  chain_id field : {chain_field or 'not found'}")
    print(f"  seq_num field  : {seq_field  or 'not found'}")
    return chain_field, seq_field


def get_locations(verified: list, comp_id_map: dict) -> list:
    """
    Queries chain_id + seq_num for verified entries using
    fields discovered by DataSchema.find_field_names().
    """
    from rcsbapi.data import DataQuery

    chain_f, seq_f = discover_location_fields()

    # Build return_data_list dynamically — only include fields that were found
    return_fields = [
        "rcsb_id",
        "nonpolymer_entities.pdbx_entity_nonpoly.comp_id",
        "nonpolymer_entities.pdbx_entity_nonpoly.name",
    ]
    if chain_f:
        return_fields.append(chain_f)
    if seq_f:
        return_fields.append(seq_f)

    rows  = []
    total = len(verified)

    for i in range(0, total, BATCH_SIZE):
        batch = verified[i : i + BATCH_SIZE]
        end   = min(i + BATCH_SIZE, total)
        print(f"  [{i+1:>4}–{end:>4} / {total}]", end="  ", flush=True)

        try:
            result = DataQuery(
                input_type="entry",
                input_ids=batch,
                return_data_list=return_fields,
                suppress_autocomplete_warning=True
            ).exec()
        except Exception as e:
            print(f"location query error: {e}")
            # Still add basic rows from comp_id_map
            for pid in batch:
                for cid in comp_id_map.get(pid, []):
                    rows.append({"pdb_id": pid, "comp_id": cid,
                                 "chem_name": "", "chain_id": "?", "seq_num": "?"})
            continue

        entries = (result.get("data") or {}).get("entries") or []
        for entry in entries:
            pdb_id = (entry.get("rcsb_id") or "").upper()
            for ent in (entry.get("nonpolymer_entities") or []):
                comp  = ent.get("pdbx_entity_nonpoly") or {}
                cid   = (comp.get("comp_id") or "").upper()
                if cid not in ALL_LN_COMP_IDS:
                    continue
                # Extract location using discovered field paths
                cid_info = {}
                if chain_f:
                    # chain_f is like "nonpolymer_entities.rcsb_nonpolymer_entity_container_identifiers.auth_asym_id"
                    # The key in ent is the last component after the entity-level prefix
                    parts = chain_f.split(".")
                    if len(parts) >= 2:
                        container = ent.get(parts[-2]) or {}
                        cid_info["chain_id"] = container.get(parts[-1], "?")
                if seq_f:
                    parts = seq_f.split(".")
                    if len(parts) >= 2:
                        container = ent.get(parts[-2]) or {}
                        cid_info["seq_num"] = container.get(parts[-1], "?")

                rows.append({
                    "pdb_id":    pdb_id,
                    "comp_id":   cid,
                    "chem_name": comp.get("name", ""),
                    "chain_id":  cid_info.get("chain_id", "?"),
                    "seq_num":   cid_info.get("seq_num",  "?"),
                })

        print(f"done ({len(batch)})")
        time.sleep(0.2)

    return rows


# ─────────────────────────────────────────────────────────────
# SAVE + MAIN
# ─────────────────────────────────────────────────────────────

def save(verified: list, rows: list, n_cand: int):
    with open(IDS_FILE, "w", encoding="utf-8") as f:
        f.write("# HackRush 2026 | Problem 9B | Checkpoint 1\n")
        f.write("# VERIFIED: only structures confirmed to contain Ln3+ ion\n")
        f.write(f"# Candidates : {n_cand} | Verified : {len(verified)}\n")
        f.write("# " + "-"*45 + "\n")
        for pid in verified:
            f.write(pid + "\n")
    print(f"  IDs file   → {IDS_FILE}")

    with open(CSV_FILE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id","comp_id","chem_name","chain_id","seq_num"])
        w.writeheader()
        w.writerows(rows)
    print(f"  Motif CSV  → {CSV_FILE}  ({len(rows)} Ln3+ sites)")


def main():
    print("\n" + "="*60)
    print("  HackRush '26 | Checkpoint 1 | Zero-FP Lanthanide PDBs")
    print("="*60)

    print("\n[1/3] Searching candidates...")
    print("-"*60)
    candidates = get_candidates()
    print(f"\n  Raw candidates : {len(candidates)}")

    print("\n[2/3] Verifying (DataQuery: comp_id check only)...")
    print("-"*60)
    verified, comp_id_map = verify(candidates)
    print(f"\n  ✓ Verified structures with Ln3+ confirmed : {len(verified)}")

    print("\n[3/3] Locating motifs (DataSchema field discovery)...")
    print("-"*60)
    rows = get_locations(verified, comp_id_map)

    save(verified, rows, len(candidates))

    print("\n" + "="*60)
    print(f"  Confirmed Ln3+ structures : {len(verified)}")
    print(f"  Lanthanide binding sites  : {len(rows)}")
    if rows:
        print()
        print(f"  {'PDB':<6} {'CompID':<8} {'Chain':<7} {'Res':<7} Name")
        print("  " + "-"*52)
        for r in rows[:10]:
            print(f"  {r['pdb_id']:<6} {r['comp_id']:<8} "
                  f"{r['chain_id']:<7} {r['seq_num']:<7} {r['chem_name'][:25]}")
    print("="*60)


if __name__ == "__main__":
    main()
