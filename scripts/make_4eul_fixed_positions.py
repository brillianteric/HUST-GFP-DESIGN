#!/usr/bin/env python3
"""Build 4EUL fixed positions for ProteinMPNN.

The fixed set combines literature-supported residues and residues whose
heavy atoms are within the chromophore shell cutoff.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


CHROMO_RESNAMES = {"CRO", "CR2", "NRQ"}
EXCLUDE_RESNAMES = {
    "HOH", "WAT", "DOD", "NA", "K", "CL", "MG", "MN", "CA", "ZN",
    "FE", "CU", "CO", "NI", "SO4", "PO4", "GOL", "EDO",
}

TARGET = "4EUL"
CHAIN = "A"

# ProteinMPNN parsed indices, 1-based.
LITERATURE_FIXED = [62, 63, 64, 65, 67, 94, 146, 201, 203, 220]


def infer_element(atom_name: str, element_field: str) -> str:
    element = element_field.strip()
    if element:
        return element.upper()
    atom_name = atom_name.strip()
    if not atom_name:
        return ""
    if atom_name[0].isdigit() and len(atom_name) > 1:
        return atom_name[1].upper()
    return atom_name[0].upper()


def dist2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b))


def parse_pdb_atoms(pdb_path: Path) -> dict:
    residues = {}
    with pdb_path.open() as handle:
        for line in handle:
            record = line[0:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip() or "_"
            if chain != CHAIN:
                continue
            try:
                resseq = int(line[22:26])
                icode = line[26].strip()
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue

            key = (chain, resseq, icode, resname)
            residues.setdefault(
                key,
                {
                    "chain": chain,
                    "resseq": resseq,
                    "icode": icode,
                    "resname": resname,
                    "atoms": {},
                    "heavy_atoms": [],
                },
            )
            residues[key]["atoms"][atom_name] = coord
            element = infer_element(atom_name, line[76:78] if len(line) >= 78 else "")
            if element != "H":
                residues[key]["heavy_atoms"].append((atom_name, coord))
    return residues


def load_parsed_4eul(parsed_jsonl: Path) -> dict:
    with parsed_jsonl.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            obj = json.loads(line)
            if obj.get("name") == TARGET:
                return obj
    raise ValueError(f"{TARGET} not found in {parsed_jsonl}")


def is_nan_coord(xyz: list[float]) -> bool:
    return any(isinstance(value, float) and math.isnan(value) for value in xyz)


def build_ca_index_map(parsed_obj: dict, residues: dict, tol: float = 0.15) -> dict:
    coords = parsed_obj[f"coords_chain_{CHAIN}"][f"CA_chain_{CHAIN}"]
    parsed_ca = [
        (idx, tuple(float(v) for v in xyz))
        for idx, xyz in enumerate(coords, start=1)
        if xyz is not None and not is_nan_coord(xyz)
    ]

    mapping = {}
    tol2 = tol * tol
    for key, residue in residues.items():
        if "CA" not in residue["atoms"]:
            continue
        best = min(parsed_ca, key=lambda item: dist2(residue["atoms"]["CA"], item[1]))
        if dist2(residue["atoms"]["CA"], best[1]) <= tol2:
            mapping[key] = best[0]
    return mapping


def is_shell_candidate(residue: dict) -> bool:
    resname = residue["resname"]
    if resname in CHROMO_RESNAMES or resname in EXCLUDE_RESNAMES:
        return False
    return "CA" in residue["atoms"] and ("N" in residue["atoms"] or "C" in residue["atoms"])


def residue_label(residue: dict) -> str:
    icode = residue["icode"] if residue["icode"] else ""
    return f"{residue['chain']}:{residue['resname']}{residue['resseq']}{icode}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, required=True, help="Prepared single-chain 4EUL PDB")
    parser.add_argument("--parsed-jsonl", type=Path, required=True)
    parser.add_argument("--output-jsonl", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument("--cutoff", type=float, default=5.0)
    args = parser.parse_args()

    parsed = load_parsed_4eul(args.parsed_jsonl)
    residues = parse_pdb_atoms(args.pdb)
    ca_map = build_ca_index_map(parsed, residues)

    chromo_atoms = [
        coord
        for residue in residues.values()
        if residue["resname"] in CHROMO_RESNAMES
        for _, coord in residue["heavy_atoms"]
    ]
    chromo_labels = [
        residue_label(residue)
        for residue in residues.values()
        if residue["resname"] in CHROMO_RESNAMES
    ]
    if not chromo_atoms:
        raise ValueError(f"No chromophore residue found in {args.pdb}")

    cutoff2 = args.cutoff * args.cutoff
    shell_positions = []
    shell_labels = []
    for key, residue in residues.items():
        if not is_shell_candidate(residue):
            continue
        if any(dist2(atom_coord, chromo_coord) <= cutoff2 for _, atom_coord in residue["heavy_atoms"] for chromo_coord in chromo_atoms):
            if key in ca_map:
                shell_positions.append(ca_map[key])
                shell_labels.append(residue_label(residue))

    shell_positions = sorted(set(shell_positions))
    positions = sorted(set(LITERATURE_FIXED) | set(shell_positions))
    sequence = parsed[f"seq_chain_{CHAIN}"]

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_jsonl.write_text(json.dumps({TARGET: {CHAIN: positions}}) + "\n")

    row = {
        "target": TARGET,
        "chain": CHAIN,
        "pdb_file": args.pdb.as_posix(),
        "chromophore": ";".join(chromo_labels),
        "cutoff_A": args.cutoff,
        "seq_len": len(sequence),
        "literature_fixed_count": len(LITERATURE_FIXED),
        "literature_fixed_positions": " ".join(map(str, LITERATURE_FIXED)),
        "shell_count": len(shell_positions),
        "shell_positions": " ".join(map(str, shell_positions)),
        "shell_labels": ";".join(shell_labels),
        "total_count": len(positions),
        "positions": " ".join(map(str, positions)),
        "aas": "".join(sequence[pos - 1] for pos in positions),
    }
    with args.summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)

    print(f"4EUL fixed positions: {len(positions)}")
    print(f"Wrote {args.output_jsonl}")
    print(f"Wrote {args.summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
