#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse
import csv
import json
import math
import re
import numpy as np

from Bio.PDB import PDBParser, Superimposer


# Default paths for the 4EUL structure-metric workflow.
DESIGN_PDB_DIR = Path("results/example_models")
TARGET_PDB_DIR = Path("data/raw")
CONFIDENCE_DIR = Path("outputs/af3_stage2")
WT_SCORE_FASTA = Path("data/proteinmpnn/4EUL_designs.fa")

OUT_CSV = Path("results/metrics/4EUL_structure_metrics_with_plddt_pae.csv")

# Design PDB file format: {DESIGN_PREFIX}{pdb_id}_{rank}_{score}_model.pdb
DESIGN_PREFIX = ""
DESIGN_SUFFIX = "_model.pdb"

# Preferred chain for CA extraction and alignment.
PREFERRED_CHAIN_ID = "A"

# 4EUL-only: exclude the chromophore-forming tripeptide from the design CA list.
EXCLUDE_DESIGN_POSITIONS = {
    "4EUL": [63, 64, 65],
}

# If lengths still differ after exclusion, truncate to the shared length.
ALLOW_TRUNCATE_IF_LENGTH_MISMATCH = True


def get_first_model(structure):
    return next(structure.get_models())


def choose_chain(model, preferred_chain_id="A"):
    chains = list(model.get_chains())
    if not chains:
        raise ValueError("No chains found")

    for chain in chains:
        if chain.id == preferred_chain_id:
            return chain

    return chains[0]


def get_ca_coords(pdb_path: Path, preferred_chain_id="A"):
    """
    Extract CA coordinates from the preferred chain in the first PDB model.

    Water, ions, and common small molecules are skipped. The parser avoids a
    strict standard-amino-acid whitelist so residue order stays compatible with
    chromophore-aware preprocessing.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = get_first_model(structure)
    chain = choose_chain(model, preferred_chain_id)

    coords = []
    residue_ids = []

    SKIP_RESNAMES = {
        "HOH", "WAT", "DOD",
        "CA", "MG", "MN", "NA", "K", "CL", "ZN", "FE", "CU", "CO", "NI",
        "SO4", "PO4", "GOL", "EDO"
    }

    for residue in chain.get_residues():
        hetflag, resseq, icode = residue.id
        resname = residue.get_resname().strip()

        if resname in SKIP_RESNAMES:
            continue

        if "CA" not in residue:
            continue

        ca = residue["CA"]
        coords.append(ca.get_coord())
        residue_ids.append((hetflag, resseq, icode.strip(), resname))

    coords = np.array(coords, dtype=float)

    return {
        "chain_id": chain.id,
        "coords": coords,
        "residue_ids": residue_ids,
    }


def remove_positions_by_1based(coords: np.ndarray, remove_positions: list[int]):
    """Remove 1-based sequence positions from a CA coordinate array."""
    if len(coords) == 0:
        return coords

    remove_set = set(remove_positions)
    keep_indices = [
        i for i in range(len(coords))
        if (i + 1) not in remove_set
    ]

    return coords[keep_indices]


def kabsch_superpose_rmsd(mobile_coords, target_coords):
    """Superpose mobile coordinates onto target coordinates and return RMSD."""

    if len(mobile_coords) != len(target_coords):
        raise ValueError(
            f"mobile and target length mismatch: "
            f"{len(mobile_coords)} vs {len(target_coords)}"
        )

    if len(mobile_coords) < 3:
        raise ValueError("Need at least 3 CA atoms for superposition")

    class FakeAtom:
        def __init__(self, coord):
            self.coord = np.array(coord, dtype=float)

        def get_coord(self):
            return self.coord

        def transform(self, rot, tran):
            self.coord = np.dot(self.coord, rot) + tran

    fixed_atoms = [FakeAtom(c) for c in target_coords]
    moving_atoms = [FakeAtom(c) for c in mobile_coords]

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    sup.apply(moving_atoms)

    transformed = np.array([a.get_coord() for a in moving_atoms], dtype=float)

    return float(sup.rms), transformed


def approximate_tm_score(transformed_mobile_coords, target_coords, norm_len):
    """
    Approximate a CA TM-score after sequential residue alignment and Kabsch
    superposition. This is not a replacement for official TMscore/US-align.
    """
    if len(target_coords) == 0:
        return float("nan")

    L = max(int(norm_len), 1)

    if L > 15:
        d0 = 1.24 * ((L - 15) ** (1.0 / 3.0)) - 1.8
        d0 = max(d0, 0.5)
    else:
        d0 = 0.5

    dists = np.linalg.norm(transformed_mobile_coords - target_coords, axis=1)
    tm_terms = 1.0 / (1.0 + (dists / d0) ** 2)

    return float(np.sum(tm_terms) / L)


def parse_design_filename(design_path: Path, target_ids: list[str]):
    """
    Parse target ID, rank, and MPNN score from a design PDB filename.

    The matching is case-insensitive against target PDB IDs.
    """
    name = design_path.name

    if not name.startswith(DESIGN_PREFIX) or not name.endswith(DESIGN_SUFFIX):
        return None

    core = name[len(DESIGN_PREFIX):-len(DESIGN_SUFFIX)]
    # core = {pdb_id}_{rank}_{score}

    target_id_map = {tid.lower(): tid for tid in target_ids}
    target_ids_lower = sorted(target_id_map.keys(), key=len, reverse=True)

    matched_lower = None
    rest = None

    for tid_lower in target_ids_lower:
        prefix = tid_lower + "_"
        if core.lower().startswith(prefix):
            matched_lower = tid_lower
            rest = core[len(prefix):]
            break

    if matched_lower is None:
        return None

    matched_pdb_id = target_id_map[matched_lower]

    rank = ""
    score = ""

    for part in rest.split("_"):
        part_lower = part.lower()
        if part_lower.startswith("rank"):
            rank = part
        elif part_lower.startswith("score"):
            score = part[5:]

    return {
        "pdb_id": matched_pdb_id,
        "rank": rank,
        "mpnn_score_from_name": score,
        "core": core,
    }


def compare_design_to_target_excluding_chromophore(
    design_pdb: Path,
    target_pdb: Path,
    pdb_id: str,
    preferred_chain_id: str,
):
    """
    Remove CA positions corresponding to the chromophore-forming tripeptide
    from the design structure, then compare sequentially to the target PDB.
    """
    target = get_ca_coords(target_pdb, preferred_chain_id)
    design = get_ca_coords(design_pdb, preferred_chain_id)

    target_coords_all = target["coords"]
    design_coords_all = design["coords"]

    target_ca_len = len(target_coords_all)
    design_ca_len = len(design_coords_all)

    exclude_positions = EXCLUDE_DESIGN_POSITIONS.get(pdb_id, [])
    design_coords_excluded = remove_positions_by_1based(
        design_coords_all,
        exclude_positions,
    )

    design_ca_len_after_exclusion = len(design_coords_excluded)

    length_match_after_exclusion = (
        target_ca_len == design_ca_len_after_exclusion
    )

    warning = ""

    if not length_match_after_exclusion:
        msg = (
            f"Length mismatch after exclusion for {pdb_id}: "
            f"target={target_ca_len}, "
            f"design_after_exclusion={design_ca_len_after_exclusion}"
        )

        if not ALLOW_TRUNCATE_IF_LENGTH_MISMATCH:
            raise ValueError(msg)

        common_len = min(target_ca_len, design_ca_len_after_exclusion)
        warning = msg + f"; truncated_to={common_len}"

        target_coords = target_coords_all[:common_len]
        design_coords = design_coords_excluded[:common_len]
    else:
        common_len = target_ca_len
        target_coords = target_coords_all
        design_coords = design_coords_excluded

    rmsd, design_superposed = kabsch_superpose_rmsd(
        mobile_coords=design_coords,
        target_coords=target_coords,
    )

    tm_score = approximate_tm_score(
        transformed_mobile_coords=design_superposed,
        target_coords=target_coords,
        norm_len=target_ca_len,
    )

    return {
        "target_chain": target["chain_id"],
        "design_chain": design["chain_id"],

        "target_ca_len": target_ca_len,
        "design_ca_len": design_ca_len,
        "excluded_design_positions": " ".join(map(str, exclude_positions)),
        "design_ca_len_after_exclusion": design_ca_len_after_exclusion,
        "length_match_after_exclusion": length_match_after_exclusion,
        "common_ca_len": common_len,

        "ca_rmsd_excluding_chromophore": rmsd,
        "ca_tm_score_approx_excluding_chromophore": tm_score,

        "warning": warning,
    }


def read_wt_mpnn_score(fasta_path: Path) -> str:
    if not fasta_path.exists():
        return ""

    with fasta_path.open() as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                continue
            match = re.search(r"score=([0-9.]+)", line)
            return match.group(1) if match else ""
    return ""


def parse_rank_num(value) -> str:
    match = re.search(r"(\d+)", str(value))
    return str(int(match.group(1))) if match else str(value)


def parse_confidence_rank(path: Path):
    match = re.search(r"4eul_rank(\d+)_score[0-9.]+_confidences\.json$", path.name, re.IGNORECASE)
    if not match:
        return None
    return parse_rank_num(match.group(1))


def extract_confidence_metrics(conf_json: Path):
    with conf_json.open("r", encoding="utf-8") as f:
        data = json.load(f)

    atom_plddt_mean = ""
    pae_mean = ""

    if "atom_plddts" in data:
        values = np.asarray(data["atom_plddts"], dtype=float)
        values = values[np.isfinite(values)]
        atom_plddt_mean = float(np.mean(values)) if values.size else ""

    if "pae" in data:
        values = np.asarray(data["pae"], dtype=float)
        values = values[np.isfinite(values)]
        pae_mean = float(np.mean(values)) if values.size else ""

    score_match = re.search(r"score([0-9.]+)_confidences\.json$", conf_json.name, re.IGNORECASE)
    score_from_conf_name = score_match.group(1) if score_match else ""

    return {
        "atom_plddt_mean": atom_plddt_mean,
        "pae_mean": pae_mean,
        "confidence_json": conf_json.as_posix(),
        "score_from_conf_name": score_from_conf_name,
        "confidence_matched": 1,
    }


def build_confidence_map(confidence_dir: Path):
    if not confidence_dir.exists():
        print(f"[WARN] Confidence dir not found, pLDDT/PAE will be blank: {confidence_dir}")
        return {}

    conf_map = {}
    for conf_json in sorted(confidence_dir.rglob("*_confidences.json")):
        rank = parse_confidence_rank(conf_json)
        if rank is None:
            continue
        try:
            conf_map[rank] = extract_confidence_metrics(conf_json)
        except Exception as exc:
            print(f"[WARN] Failed to parse confidence JSON {conf_json}: {exc}")

    print(f"[INFO] Matched AF3 confidence JSON files: {len(conf_map)}")
    return conf_map


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate 4EUL structure metrics after excluding the chromophore-forming tripeptide."
    )
    parser.add_argument("--design-pdb-dir", type=Path, default=DESIGN_PDB_DIR)
    parser.add_argument("--target-pdb-dir", type=Path, default=TARGET_PDB_DIR)
    parser.add_argument("--confidence-dir", type=Path, default=CONFIDENCE_DIR)
    parser.add_argument("--wt-score-fasta", type=Path, default=WT_SCORE_FASTA)
    parser.add_argument("--out-csv", type=Path, default=OUT_CSV)
    parser.add_argument("--preferred-chain", default=PREFERRED_CHAIN_ID)
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.design_pdb_dir.exists():
        raise FileNotFoundError(f"Design PDB dir not found: {args.design_pdb_dir}")

    if not args.target_pdb_dir.exists():
        raise FileNotFoundError(f"Target PDB dir not found: {args.target_pdb_dir}")

    target_pdbs = {
        p.stem: p
        for p in sorted(args.target_pdb_dir.glob("*.pdb"))
    }

    target_ids = sorted(target_pdbs.keys())

    if not target_pdbs:
        raise FileNotFoundError(f"No target PDB files found in: {args.target_pdb_dir}")

    design_pdbs = sorted(args.design_pdb_dir.glob("*.pdb"))

    print("========== INPUT ==========")
    print(f"[INFO] Design PDB dir : {args.design_pdb_dir}")
    print(f"[INFO] Target PDB dir : {args.target_pdb_dir}")
    print(f"[INFO] Output CSV     : {args.out_csv}")
    print(f"[INFO] Target PDBs    : {len(target_pdbs)}")
    print(f"[INFO] Design PDBs    : {len(design_pdbs)}")
    print(f"[INFO] Method         : exclude chromophore tripeptide from design CA, then sequential CA alignment")

    wt_mpnn_score = read_wt_mpnn_score(args.wt_score_fasta)
    confidence_by_rank = build_confidence_map(args.confidence_dir)

    rows = []

    n_ok = 0
    n_fail = 0
    n_skip_name = 0
    n_missing_target = 0

    for design_pdb in design_pdbs:
        parsed = parse_design_filename(design_pdb, target_ids)

        if parsed is None:
            n_skip_name += 1
            print(f"[SKIP] filename not matched: {design_pdb.name}")
            continue

        pdb_id = parsed["pdb_id"]
        target_pdb = target_pdbs.get(pdb_id)

        if target_pdb is None:
            n_missing_target += 1
            print(f"[WARN] missing target for {design_pdb.name}: pdb_id={pdb_id}")
            continue

        try:
            metrics = compare_design_to_target_excluding_chromophore(
                design_pdb=design_pdb,
                target_pdb=target_pdb,
                pdb_id=pdb_id,
                preferred_chain_id=args.preferred_chain,
            )

            row = {
                "pdb_id": pdb_id,
                "wt_mpnn_score": wt_mpnn_score,
                "rank": parsed["rank"],
                "mpnn_score_from_name": parsed["mpnn_score_from_name"],
                **metrics,
                "design_file": design_pdb.as_posix(),
                "target_file": target_pdb.as_posix(),
            }
            row.update(
                confidence_by_rank.get(
                    parse_rank_num(parsed["rank"]),
                    {
                        "atom_plddt_mean": "",
                        "pae_mean": "",
                        "confidence_json": "",
                        "score_from_conf_name": "",
                        "confidence_matched": 0,
                    },
                )
            )

            rows.append(row)
            n_ok += 1

            print(
                f"[OK] {design_pdb.name} vs {target_pdb.name} | "
                f"RMSD={metrics['ca_rmsd_excluding_chromophore']:.4f} | "
                f"TM~={metrics['ca_tm_score_approx_excluding_chromophore']:.4f} | "
                f"target_CA={metrics['target_ca_len']} | "
                f"design_CA={metrics['design_ca_len']} -> "
                f"{metrics['design_ca_len_after_exclusion']}"
            )

            if metrics["warning"]:
                print(f"     [WARN] {metrics['warning']}")

        except Exception as e:
            n_fail += 1
            print(f"[FAIL] {design_pdb.name}: {e}")

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "pdb_id",
        "wt_mpnn_score",
        "rank",
        "mpnn_score_from_name",

        "ca_rmsd_excluding_chromophore",
        "ca_tm_score_approx_excluding_chromophore",

        "common_ca_len",
        "target_ca_len",
        "design_ca_len",
        "excluded_design_positions",
        "design_ca_len_after_exclusion",
        "length_match_after_exclusion",

        "target_chain",
        "design_chain",
        "warning",

        "design_file",
        "target_file",

        "atom_plddt_mean",
        "pae_mean",
        "confidence_json",
        "score_from_conf_name",
        "confidence_matched",
    ]

    with args.out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("========== DONE ==========")
    print(f"[INFO] OK             : {n_ok}")
    print(f"[INFO] Failed         : {n_fail}")
    print(f"[INFO] Skip name      : {n_skip_name}")
    print(f"[INFO] Missing target : {n_missing_target}")
    print(f"[INFO] Output CSV     : {args.out_csv}")


if __name__ == "__main__":
    main()
