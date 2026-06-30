#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO


# =========================
# 鐢ㄦ埛鍙慨鏀瑰弬鏁?# =========================
AF3_OUTPUT_DIR = Path("outputs/af3_stage2")
OUT_PDB_DIR = Path("results/example_models")

# 鍙浆鎹?AF3 涓昏緭鍑虹殑 *_model.cif
CIF_PATTERN = "*_model.cif"

# 杈撳嚭鏂囦欢鍚嶅墠缂€锛涗笉鎯冲姞鍓嶇紑灏辫涓?""
OUTPUT_PREFIX = ""

# 濡傛灉杈撳嚭 PDB 宸插瓨鍦紝鏄惁璺宠繃
SKIP_EXISTING = True


def find_cif_files(af3_output_dir: Path) -> list[Path]:
    """
    閫掑綊鎵弿 AF3 杈撳嚭鐩綍锛屾壘鍒版墍鏈?*_model.cif銆?    """
    if not af3_output_dir.exists():
        raise FileNotFoundError(f"AF3 output directory not found: {af3_output_dir}")

    cif_files = sorted(af3_output_dir.rglob(CIF_PATTERN))

    print("========== CIF SCAN ==========")
    print(f"[INFO] AF3 output dir : {af3_output_dir}")
    print(f"[INFO] CIF pattern    : {CIF_PATTERN}")
    print(f"[INFO] Found CIF files: {len(cif_files)}")

    return cif_files


def convert_cif_to_pdb(cif_path: Path, out_pdb_path: Path, parser: MMCIFParser, io: PDBIO):
    """
    Convert one mmCIF file to PDB using Biopython.
    """
    structure_id = cif_path.stem
    structure = parser.get_structure(structure_id, str(cif_path))
    io.set_structure(structure)
    io.save(str(out_pdb_path))


def main():
    cif_files = find_cif_files(AF3_OUTPUT_DIR)

    OUT_PDB_DIR.mkdir(parents=True, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    converted = 0
    skipped_existing = 0
    failed = 0

    print("========== CONVERT ==========")

    for cif_path in cif_files:
        # 渚嬪 seq_1_model.cif -> seq_1
        sequence_id = cif_path.name.replace("_model.cif", "")

        out_pdb_name = f"{OUTPUT_PREFIX}{sequence_id}_model.pdb"
        out_pdb_path = OUT_PDB_DIR / out_pdb_name

        if SKIP_EXISTING and out_pdb_path.exists():
            skipped_existing += 1
            print(f"[SKIP] exists: {out_pdb_path}")
            continue

        try:
            convert_cif_to_pdb(cif_path, out_pdb_path, parser, io)
            converted += 1
            print(f"[OK] {sequence_id}: {cif_path} -> {out_pdb_path}")
        except Exception as e:
            failed += 1
            print(f"[WARN] Failed: sequence_id={sequence_id} | {cif_path} | {e}")

    print("========== DONE ==========")
    print(f"[INFO] Found CIF files : {len(cif_files)}")
    print(f"[INFO] Converted       : {converted}")
    print(f"[INFO] Skipped existing: {skipped_existing}")
    print(f"[INFO] Failed          : {failed}")
    print(f"[INFO] Output dir      : {OUT_PDB_DIR}")


if __name__ == "__main__":
    main()
