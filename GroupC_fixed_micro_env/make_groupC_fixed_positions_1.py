#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate GroupC fixed positions for chromophore-aware ProteinMPNN design.

GroupC = GroupB fixed positions + chromophore-centered shell residues.

Input:
  1. parsed_pdbs jsonl from chromophore-aware parser
  2. original single-chain PDB folder

Output:
  1. ProteinMPNN fixed_positions_jsonl
  2. summary CSV for checking fixed residues

Important:
  Output positions are ProteinMPNN parsed chain indices, 1-based.
"""

import argparse
import csv
import json
import math
from pathlib import Path


CHROMO_RESNAMES = {"CRO", "CR2", "NRQ"}
WATER_RESNAMES = {"HOH", "WAT", "DOD"}

# Common ions / small molecules to exclude from protein-shell candidates
EXCLUDE_RESNAMES = {
    "HOH", "WAT", "DOD",
    "NA", "K", "CL", "MG", "MN", "CA", "ZN", "FE", "CU", "CO", "NI",
    "SO4", "PO4", "GOL", "EDO"
}

# GroupB fixed positions: ProteinMPNN parsed indices, 1-based
GROUP_B_FIXED = {
    "1EMB": {
        "A": [64, 65, 66, 68, 95, 147, 202, 204, 221]
    },
    "1HUY": {
        "A": [68, 69, 70, 72, 99, 151, 206, 208, 225]
    },
    "1MYW": {
        "A": [47, 65, 66, 67, 68, 70, 97, 149, 154, 164, 176, 204, 206, 223]
    },
    "3BXB_A": {
        "A": [60, 61, 62, 140, 171, 194]
    },
    "4EUL": {
        "A": [62, 63, 64, 65, 67, 94, 146, 201, 203, 220]
    },
    "4OQW_A": {
        "A": [28, 41, 63, 64, 65, 143, 158, 197]
    }
}


def infer_element(atom_name: str, element_field: str) -> str:
    element = element_field.strip()
    if element:
        return element.upper()

    # Infer from atom name if element column is empty
    s = atom_name.strip()
    if not s:
        return ""
    # Atom names like CA1 / C1 / N2 / OXT
    if len(s) >= 2 and s[0].isdigit():
        return s[1].upper()
    return s[0].upper()


def dist2(a, b):
    return (
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


def parse_pdb_atoms(pdb_path: Path):
    """
    Return residues dict:
      key = (chain, resseq, icode, resname)
      value = {
        "chain": ...,
        "resseq": int,
        "icode": str,
        "resname": str,
        "atoms": {atom_name: (x,y,z)},
        "heavy_atoms": [(atom_name, (x,y,z)), ...],
        "records": set(["ATOM", "HETATM"])
      }
    """
    residues = {}

    with pdb_path.open() as f:
        for line in f:
            record = line[0:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue

            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip() or "_"
            resseq_str = line[22:26].strip()
            icode = line[26].strip()

            if not resseq_str:
                continue

            try:
                resseq = int(resseq_str)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            element = infer_element(atom_name, line[76:78] if len(line) >= 78 else "")
            key = (chain, resseq, icode, resname)

            if key not in residues:
                residues[key] = {
                    "chain": chain,
                    "resseq": resseq,
                    "icode": icode,
                    "resname": resname,
                    "atoms": {},
                    "heavy_atoms": [],
                    "records": set(),
                }

            residues[key]["records"].add(record)
            residues[key]["atoms"][atom_name] = (x, y, z)

            if element != "H":
                residues[key]["heavy_atoms"].append((atom_name, (x, y, z)))

    return residues


def load_parsed_jsonl(parsed_path: Path):
    parsed = {}
    for line in parsed_path.open():
        line = line.strip()
        if not line:
            continue
        obj = json.loads(line)
        parsed[obj["name"]] = obj
    return parsed


def find_pdb_for_name(pdb_dir: Path, name: str):
    candidates = []

    # exact
    candidates.append(pdb_dir / f"{name}.pdb")
    candidates.append(pdb_dir / f"{name.upper()}.pdb")
    candidates.append(pdb_dir / f"{name.lower()}.pdb")

    # for names like 3BXB_A, try 3BXB.pdb
    if "_" in name:
        base = name.split("_")[0]
        candidates.append(pdb_dir / f"{base}.pdb")
        candidates.append(pdb_dir / f"{base.upper()}.pdb")
        candidates.append(pdb_dir / f"{base.lower()}.pdb")

    for p in candidates:
        if p.exists():
            return p

    # fallback: case-insensitive search
    lower_map = {p.name.lower(): p for p in pdb_dir.glob("*.pdb")}
    for c in candidates:
        hit = lower_map.get(c.name.lower())
        if hit:
            return hit

    return None


def is_nan_coord(xyz):
    return any(isinstance(v, float) and math.isnan(v) for v in xyz)


def build_ca_index_map(parsed_obj, chain: str, residues, tol=0.15):
    """
    Map PDB residue key -> ProteinMPNN parsed index by CA coordinate matching.
    This avoids errors caused by different PDB numbering or missing residues.
    """
    coords = parsed_obj.get(f"coords_chain_{chain}")
    if coords is None:
        raise ValueError(f"coords_chain_{chain} not found in parsed jsonl for {parsed_obj['name']}")

    ca_list = coords.get(f"CA_chain_{chain}")
    if ca_list is None:
        raise ValueError(f"CA_chain_{chain} not found in parsed jsonl for {parsed_obj['name']}")

    parsed_ca = []
    for i, xyz in enumerate(ca_list, start=1):
        if xyz is None or is_nan_coord(xyz):
            continue
        parsed_ca.append((i, tuple(float(v) for v in xyz)))

    mapping = {}
    tol2 = tol * tol

    for key, res in residues.items():
        if res["chain"] != chain:
            continue
        if "CA" not in res["atoms"]:
            continue

        ca = res["atoms"]["CA"]

        best_i = None
        best_d2 = None

        for idx, pca in parsed_ca:
            d2 = dist2(ca, pca)
            if best_d2 is None or d2 < best_d2:
                best_d2 = d2
                best_i = idx

        if best_d2 is not None and best_d2 <= tol2:
            mapping[key] = best_i

    return mapping


def residue_label(res):
    icode = res["icode"] if res["icode"] else ""
    return f"{res['chain']}:{res['resname']}{res['resseq']}{icode}"


def get_chromophore_atoms(residues, chain=None):
    chromo_atoms = []
    chromo_labels = []

    for key, res in residues.items():
        if chain is not None and res["chain"] != chain:
            continue
        if res["resname"] in CHROMO_RESNAMES:
            chromo_labels.append(residue_label(res))
            for atom_name, xyz in res["heavy_atoms"]:
                chromo_atoms.append(xyz)

    return chromo_atoms, chromo_labels


def is_protein_like_candidate(res):
    """
    Candidate shell residue:
      - not chromophore
      - not water/ion/small molecule
      - has CA and at least N or C atom
    This keeps standard residues and common modified residues like MLY if parser handled them.
    """
    resname = res["resname"]

    if resname in CHROMO_RESNAMES:
        return False
    if resname in EXCLUDE_RESNAMES:
        return False
    if "CA" not in res["atoms"]:
        return False
    if "N" not in res["atoms"] and "C" not in res["atoms"]:
        return False

    return True


def find_shell_residues(residues, chromo_atoms, cutoff, chain=None):
    cutoff2 = cutoff * cutoff
    shell_keys = []

    for key, res in residues.items():
        if chain is not None and res["chain"] != chain:
            continue
        if not is_protein_like_candidate(res):
            continue

        hit = False
        for _, atom_xyz in res["heavy_atoms"]:
            for c_xyz in chromo_atoms:
                if dist2(atom_xyz, c_xyz) <= cutoff2:
                    hit = True
                    break
            if hit:
                break

        if hit:
            shell_keys.append(key)

    return shell_keys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", required=True, help="Folder containing original single-chain PDB files")
    parser.add_argument("--parsed_jsonl", required=True, help="Parsed jsonl generated by chromophore-aware parser")
    parser.add_argument("--output_jsonl", required=True, help="Output fixed positions jsonl for ProteinMPNN")
    parser.add_argument("--summary_csv", required=True, help="Output summary CSV")
    parser.add_argument("--cutoff", type=float, default=5.0, help="Chromophore shell cutoff in Angstrom")
    parser.add_argument("--chain", default="A", help="Chain ID, default A")
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    parsed_path = Path(args.parsed_jsonl)
    output_jsonl = Path(args.output_jsonl)
    summary_csv = Path(args.summary_csv)

    parsed = load_parsed_jsonl(parsed_path)

    fixed_groupC = {}
    summary_rows = []

    print("========== make_groupC_fixed_positions ==========")
    print(f"[INFO] pdb_dir       = {pdb_dir}")
    print(f"[INFO] parsed_jsonl  = {parsed_path}")
    print(f"[INFO] output_jsonl  = {output_jsonl}")
    print(f"[INFO] summary_csv   = {summary_csv}")
    print(f"[INFO] cutoff        = {args.cutoff} Å")
    print(f"[INFO] chain         = {args.chain}")

    for name, obj in parsed.items():
        chain = args.chain

        if name not in GROUP_B_FIXED:
            raise KeyError(f"{name} not found in GROUP_B_FIXED. Please add it.")

        if chain not in GROUP_B_FIXED[name]:
            raise KeyError(f"{name} chain {chain} not found in GROUP_B_FIXED.")

        pdb_path = find_pdb_for_name(pdb_dir, name)
        if pdb_path is None:
            raise FileNotFoundError(f"Cannot find PDB for parsed name: {name}")

        residues = parse_pdb_atoms(pdb_path)
        chromo_atoms, chromo_labels = get_chromophore_atoms(residues, chain=chain)

        if not chromo_atoms:
            raise RuntimeError(f"No chromophore residue {CHROMO_RESNAMES} found for {name} chain {chain}")

        ca_map = build_ca_index_map(obj, chain, residues)

        shell_keys = find_shell_residues(
            residues=residues,
            chromo_atoms=chromo_atoms,
            cutoff=args.cutoff,
            chain=chain,
        )

        shell_positions = []
        unmapped_shell = []

        for key in shell_keys:
            if key in ca_map:
                shell_positions.append(ca_map[key])
            else:
                unmapped_shell.append(residue_label(residues[key]))

        shell_positions = sorted(set(shell_positions))

        groupB_positions = sorted(set(GROUP_B_FIXED[name][chain]))
        groupC_positions = sorted(set(groupB_positions) | set(shell_positions))
        shell_extra_positions = sorted(set(shell_positions) - set(groupB_positions))

        fixed_groupC[name] = {chain: groupC_positions}

        seq = obj.get(f"seq_chain_{chain}", "")
        seq_len = len(seq)

        # Validate positions
        bad = [p for p in groupC_positions if p < 1 or p > seq_len]
        if bad:
            raise ValueError(f"{name}: out-of-range GroupC positions {bad}, seq_len={seq_len}")

        groupB_aas = "".join(seq[p - 1] for p in groupB_positions)
        shell_extra_aas = "".join(seq[p - 1] for p in shell_extra_positions)
        groupC_aas = "".join(seq[p - 1] for p in groupC_positions)

        shell_labels = []
        shell_extra_labels = []

        pos_to_label = {}
        for key in shell_keys:
            if key in ca_map:
                pos_to_label[ca_map[key]] = residue_label(residues[key])

        for p in shell_positions:
            shell_labels.append(pos_to_label.get(p, f"parsed_pos_{p}"))

        for p in shell_extra_positions:
            shell_extra_labels.append(pos_to_label.get(p, f"parsed_pos_{p}"))

        print()
        print(f"[TARGET] {name}")
        print(f"  pdb_path              = {pdb_path}")
        print(f"  chromophore           = {', '.join(chromo_labels)}")
        print(f"  GroupB count          = {len(groupB_positions)}")
        print(f"  shell count           = {len(shell_positions)}")
        print(f"  shell extra count     = {len(shell_extra_positions)}")
        print(f"  GroupC total count    = {len(groupC_positions)}")
        print(f"  shell extra positions = {shell_extra_positions}")
        if unmapped_shell:
            print(f"  [WARN] unmapped shell residues = {unmapped_shell}")

        summary_rows.append({
            "target": name,
            "chain": chain,
            "pdb_file": str(pdb_path),
            "chromophore": ";".join(chromo_labels),
            "cutoff_A": args.cutoff,
            "seq_len": seq_len,
            "groupB_count": len(groupB_positions),
            "groupB_positions": " ".join(map(str, groupB_positions)),
            "groupB_aas": groupB_aas,
            "shell_count": len(shell_positions),
            "shell_positions": " ".join(map(str, shell_positions)),
            "shell_labels": ";".join(shell_labels),
            "shell_extra_count": len(shell_extra_positions),
            "shell_extra_positions": " ".join(map(str, shell_extra_positions)),
            "shell_extra_labels": ";".join(shell_extra_labels),
            "shell_extra_aas": shell_extra_aas,
            "groupC_total_count": len(groupC_positions),
            "groupC_positions": " ".join(map(str, groupC_positions)),
            "groupC_aas": groupC_aas,
            "unmapped_shell_residues": ";".join(unmapped_shell),
        })

    output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_jsonl.open("w") as f:
        f.write(json.dumps(fixed_groupC) + "\n")

    fieldnames = list(summary_rows[0].keys())
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print()
    print("========== DONE ==========")
    print(f"[OK] wrote fixed positions jsonl: {output_jsonl}")
    print(f"[OK] wrote summary csv          : {summary_csv}")


if __name__ == "__main__":
    main()