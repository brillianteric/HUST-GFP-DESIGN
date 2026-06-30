#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Convert AlphaFold 3 mmCIF model outputs to PDB files."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO


DEFAULT_AF3_OUTPUT_DIR = Path("outputs/af3_stage2")
DEFAULT_OUT_PDB_DIR = Path("results/example_models")
DEFAULT_CIF_PATTERN = "*_model.cif"
DEFAULT_OUTPUT_PREFIX = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--af3-output-dir", type=Path, default=DEFAULT_AF3_OUTPUT_DIR)
    parser.add_argument("--out-pdb-dir", type=Path, default=DEFAULT_OUT_PDB_DIR)
    parser.add_argument("--cif-pattern", default=DEFAULT_CIF_PATTERN)
    parser.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing PDB files.")
    return parser.parse_args()


def find_cif_files(af3_output_dir: Path, cif_pattern: str) -> list[Path]:
    """Recursively find AF3 model mmCIF files."""
    if not af3_output_dir.exists():
        raise FileNotFoundError(f"AF3 output directory not found: {af3_output_dir}")

    cif_files = sorted(af3_output_dir.rglob(cif_pattern))

    print("========== CIF SCAN ==========")
    print(f"[INFO] AF3 output dir : {af3_output_dir}")
    print(f"[INFO] CIF pattern    : {cif_pattern}")
    print(f"[INFO] Found CIF files: {len(cif_files)}")

    return cif_files


def convert_cif_to_pdb(cif_path: Path, out_pdb_path: Path, parser: MMCIFParser, io: PDBIO) -> None:
    """Convert one mmCIF structure to PDB using Biopython."""
    structure_id = cif_path.stem
    structure = parser.get_structure(structure_id, str(cif_path))
    io.set_structure(structure)
    io.save(str(out_pdb_path))


def main() -> int:
    args = parse_args()
    cif_files = find_cif_files(args.af3_output_dir, args.cif_pattern)
    args.out_pdb_dir.mkdir(parents=True, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    converted = 0
    skipped_existing = 0
    failed = 0

    print("========== CONVERT ==========")

    for cif_path in cif_files:
        sequence_id = cif_path.name.replace("_model.cif", "")
        out_pdb_name = f"{args.output_prefix}{sequence_id}_model.pdb"
        out_pdb_path = args.out_pdb_dir / out_pdb_name

        if not args.overwrite and out_pdb_path.exists():
            skipped_existing += 1
            print(f"[SKIP] exists: {out_pdb_path}")
            continue

        try:
            convert_cif_to_pdb(cif_path, out_pdb_path, parser, io)
            converted += 1
            print(f"[OK] {sequence_id}: {cif_path} -> {out_pdb_path}")
        except Exception as exc:  # noqa: BLE001 - command-line utility should continue converting
            failed += 1
            print(f"[WARN] Failed: sequence_id={sequence_id} | {cif_path} | {exc}")

    print("========== DONE ==========")
    print(f"[INFO] Found CIF files : {len(cif_files)}")
    print(f"[INFO] Converted       : {converted}")
    print(f"[INFO] Skipped existing: {skipped_existing}")
    print(f"[INFO] Failed          : {failed}")
    print(f"[INFO] Output dir      : {args.out_pdb_dir}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
