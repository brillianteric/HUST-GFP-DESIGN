#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
璁＄畻鑽у厜铔嬬櫧璁捐缁撴瀯鐨勫彂鑹插洟寰幆澧冩寚鏍囷拷?
鏍稿績鎸囨爣锟?1. rmsd_shell5_ca
2. rmsd_shell8_ca
3. chromo_clash_count_2A
4. chromo_contact_F1
5. shell5_mutation_count / shell5_mutation_rate

杈撳叆锟?- 鍘熷 PDB锛氬繀椤讳繚鐣欐垚鐔熷彂鑹插洟 CRO / CR2 / NRQ
- AF3 棰勬祴 PDB锛氬凡杞垚 PDB 鐨勯娴嬬粨锟?"""

import csv
import math
import re
from pathlib import Path

import numpy as np


# ============================================================
# 鐢ㄦ埛鍙傛暟鍖猴細鍙渶瑕佹敼杩欓噷
# ============================================================

# User parameters adapted for this repository.
NATIVE_DIR = Path("data/raw")
PRED_DIR = Path("results/example_models")
OUT_CSV = Path("results/metrics/4EUL_microenv_metrics.csv")
DESIGN_SET = "4EUL"
RECURSIVE = True
CHAIN_ID = "A"
SHELL5_CUTOFF = 5.0
SHELL8_CUTOFF = 8.0
CLASH_CUTOFF = 2.0
CONTACT_CUTOFF = 4.0


# ============================================================
# 椤圭洰鍥哄畾淇℃伅
# ============================================================

CHROMO_RESNAMES = {"CRO", "CR2", "NRQ"}

TARGET_ALIASES = {
    "4eul": "4EUL",
}

# 锟?chromophore-aware parser 瑙ｆ瀽鍚庣殑鍙戣壊鍥㈠墠浣撲笁鑲戒綅缃紝1-based
CHROMO_POSITIONS = {
    "4EUL": [63, 64, 65],
}

# Native mature chromophore restored to the parent tripeptide used by parser.
CHROMO_PRECURSOR = {
    "4EUL": "TYG",
}

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MLY": "K",
}

EXCLUDE_RESNAMES = {
    "HOH", "WAT", "DOD",
    "NA", "K", "CL", "MG", "MN", "CA", "ZN", "FE", "CU", "CO", "NI",
    "SO4", "PO4", "GOL", "EDO",
}


# ============================================================
# 鍩虹鍑芥暟
# ============================================================

def infer_element(atom_name, element_field):
    element = element_field.strip()
    if element:
        return element.upper()

    atom_name = atom_name.strip()
    if not atom_name:
        return ""

    if atom_name[0].isdigit() and len(atom_name) > 1:
        return atom_name[1].upper()

    return atom_name[0].upper()


def dist2(a, b):
    d = a - b
    return float(np.dot(d, d))


def rmsd(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if len(a) == 0:
        return float("nan")

    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def kabsch(moving, fixed):
    """
    杩斿洖 R, t锛屼娇 moving_aligned = moving @ R + t
    """
    moving = np.asarray(moving, dtype=float)
    fixed = np.asarray(fixed, dtype=float)

    mc = moving.mean(axis=0)
    fc = fixed.mean(axis=0)

    X = moving - mc
    Y = fixed - fc

    H = X.T @ Y
    U, S, Vt = np.linalg.svd(H)
    R = U @ Vt

    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = U @ Vt

    t = fc - mc @ R
    return R, t


def transform(coord, R, t):
    return coord @ R + t


def clean_float(x):
    if isinstance(x, float) and math.isnan(x):
        return ""
    return x


# ============================================================
# PDB 瑙ｆ瀽
# ============================================================

def parse_pdb_residues(pdb_path, chain_id="A"):
    residues = {}
    order = []

    with open(pdb_path) as f:
        for line in f:
            record = line[0:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue

            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip() or "_"

            if chain_id and chain != chain_id:
                continue

            try:
                resseq = int(line[22:26].strip())
                icode = line[26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            key = (chain, resseq, icode, resname)

            if key not in residues:
                residues[key] = {
                    "chain": chain,
                    "resseq": resseq,
                    "icode": icode,
                    "resname": resname,
                    "atoms": {},
                    "heavy_atoms": [],
                }
                order.append(key)

            coord = np.array([x, y, z], dtype=float)
            residues[key]["atoms"][atom_name] = coord

            element = infer_element(atom_name, line[76:78] if len(line) >= 78 else "")
            if element != "H":
                residues[key]["heavy_atoms"].append((atom_name, coord))

    return [residues[k] for k in order]


def is_protein_like(res):
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


def residue_label(res):
    icode = res["icode"] if res["icode"] else ""
    return f"{res['chain']}:{res['resname']}{res['resseq']}{icode}"


def build_native_model(native_pdb, target, chain_id="A"):
    """
    鍘熷 PDB 涓彂鑹插洟锟?CRO/CR2/NRQ锟?    杩欓噷鎶婃垚鐔熷彂鑹插洟杩樺師锟?3 锟?parsed positions锛屼互瀵归綈 ProteinMPNN/AF3 搴忓垪锟?    """
    residues = parse_pdb_residues(native_pdb, chain_id=chain_id)

    parsed = {}
    seq = []
    chromo_atoms = []
    chromo_labels = []

    idx = 0

    for res in residues:
        resname = res["resname"]

        if resname in CHROMO_RESNAMES:
            chromo_labels.append(residue_label(res))

            for _, coord in res["heavy_atoms"]:
                chromo_atoms.append(coord)

            precursor = CHROMO_PRECURSOR[target]

            for j, aa in enumerate(precursor, start=1):
                idx += 1
                atoms = {}
                heavy_atoms = []

                for base in ["N", "CA", "C", "O"]:
                    atom_name = f"{base}{j}"
                    if atom_name in res["atoms"]:
                        atoms[base] = res["atoms"][atom_name]
                        heavy_atoms.append((base, res["atoms"][atom_name]))

                parsed[idx] = {
                    "aa": aa,
                    "label": f"{residue_label(res)}_{j}",
                    "atoms": atoms,
                    "heavy_atoms": heavy_atoms,
                    "is_chromo": True,
                }
                seq.append(aa)

            continue

        if not is_protein_like(res):
            continue

        idx += 1
        aa = AA3_TO_AA1.get(resname, "X")

        parsed[idx] = {
            "aa": aa,
            "label": residue_label(res),
            "atoms": res["atoms"],
            "heavy_atoms": res["heavy_atoms"],
            "is_chromo": False,
        }
        seq.append(aa)

    if not chromo_atoms:
        raise RuntimeError(f"No chromophore CRO/CR2/NRQ found in {native_pdb}")

    return {
        "parsed": parsed,
        "sequence": "".join(seq),
        "chromo_atoms": chromo_atoms,
        "chromo_labels": chromo_labels,
        "chromo_positions": CHROMO_POSITIONS[target],
    }


def build_pred_model(pred_pdb, chain_id="A"):
    residues = parse_pdb_residues(pred_pdb, chain_id=chain_id)

    parsed = {}
    seq = []

    idx = 0

    for res in residues:
        if not is_protein_like(res):
            continue

        idx += 1
        resname = res["resname"]
        aa = AA3_TO_AA1.get(resname, "X")

        parsed[idx] = {
            "aa": aa,
            "label": residue_label(res),
            "atoms": res["atoms"],
            "heavy_atoms": res["heavy_atoms"],
        }
        seq.append(aa)

    return {
        "parsed": parsed,
        "sequence": "".join(seq),
    }


# ============================================================
# 鎸囨爣璁＄畻
# ============================================================

def get_shell_positions(native_model, cutoff):
    parsed = native_model["parsed"]
    chromo_atoms = native_model["chromo_atoms"]
    chromo_positions = set(native_model["chromo_positions"])

    cutoff2 = cutoff * cutoff
    shell_positions = []

    for pos, res in parsed.items():
        if pos in chromo_positions:
            continue

        hit = False

        for _, atom_coord in res["heavy_atoms"]:
            for chromo_coord in chromo_atoms:
                if dist2(atom_coord, chromo_coord) <= cutoff2:
                    hit = True
                    break

            if hit:
                break

        if hit:
            shell_positions.append(pos)

    return sorted(shell_positions)


def get_alignment(native_model, pred_model):
    native_parsed = native_model["parsed"]
    pred_parsed = pred_model["parsed"]
    chromo_positions = set(native_model["chromo_positions"])

    moving = []
    fixed = []
    used_positions = []

    max_pos = min(max(native_parsed), max(pred_parsed))

    for pos in range(1, max_pos + 1):
        if pos in chromo_positions:
            continue

        nres = native_parsed.get(pos)
        pres = pred_parsed.get(pos)

        if nres is None or pres is None:
            continue

        if "CA" not in nres["atoms"] or "CA" not in pres["atoms"]:
            continue

        moving.append(pres["atoms"]["CA"])
        fixed.append(nres["atoms"]["CA"])
        used_positions.append(pos)

    if len(moving) < 3:
        raise RuntimeError("Too few CA pairs for alignment")

    moving = np.asarray(moving)
    fixed = np.asarray(fixed)

    R, t = kabsch(moving, fixed)
    return R, t, used_positions


def calc_shell_ca_rmsd(native_model, pred_model, positions, R, t):
    native_coords = []
    pred_coords = []

    for pos in positions:
        nres = native_model["parsed"].get(pos)
        pres = pred_model["parsed"].get(pos)

        if nres is None or pres is None:
            continue

        if "CA" not in nres["atoms"] or "CA" not in pres["atoms"]:
            continue

        native_coords.append(nres["atoms"]["CA"])
        pred_coords.append(transform(pres["atoms"]["CA"], R, t))

    return rmsd(native_coords, pred_coords), len(native_coords)


def get_pred_heavy_atoms_excluding_chromo(pred_model, chromo_positions, R, t):
    atoms = []

    for pos, res in pred_model["parsed"].items():
        if pos in chromo_positions:
            continue

        for atom_name, coord in res["heavy_atoms"]:
            atoms.append((pos, atom_name, transform(coord, R, t)))

    return atoms


def calc_chromo_clash(native_model, pred_model, R, t, clash_cutoff=2.0):
    chromo_atoms = native_model["chromo_atoms"]
    chromo_positions = set(native_model["chromo_positions"])

    pred_atoms = get_pred_heavy_atoms_excluding_chromo(
        pred_model=pred_model,
        chromo_positions=chromo_positions,
        R=R,
        t=t,
    )

    cutoff2 = clash_cutoff * clash_cutoff
    clash_count = 0
    min_dist = float("inf")

    for c in chromo_atoms:
        for _, _, p in pred_atoms:
            d2 = dist2(c, p)
            d = math.sqrt(d2)

            if d < min_dist:
                min_dist = d

            if d2 < cutoff2:
                clash_count += 1

    if min_dist == float("inf"):
        min_dist = float("nan")

    return clash_count, min_dist


def contact_positions_to_chromophore(model, chromo_atoms, cutoff=4.0, R=None, t=None, exclude_positions=None):
    exclude_positions = set(exclude_positions or [])
    cutoff2 = cutoff * cutoff
    contacts = set()

    for pos, res in model["parsed"].items():
        if pos in exclude_positions:
            continue

        hit = False

        for _, coord in res["heavy_atoms"]:
            coord2 = transform(coord, R, t) if R is not None and t is not None else coord

            for c in chromo_atoms:
                if dist2(coord2, c) <= cutoff2:
                    hit = True
                    break

            if hit:
                break

        if hit:
            contacts.add(pos)

    return contacts


def calc_contact_f1(native_model, pred_model, R, t, contact_cutoff=4.0):
    chromo_atoms = native_model["chromo_atoms"]
    chromo_positions = set(native_model["chromo_positions"])

    native_contacts = contact_positions_to_chromophore(
        model=native_model,
        chromo_atoms=chromo_atoms,
        cutoff=contact_cutoff,
        R=None,
        t=None,
        exclude_positions=chromo_positions,
    )

    pred_contacts = contact_positions_to_chromophore(
        model=pred_model,
        chromo_atoms=chromo_atoms,
        cutoff=contact_cutoff,
        R=R,
        t=t,
        exclude_positions=chromo_positions,
    )

    inter = native_contacts & pred_contacts

    precision = len(inter) / len(pred_contacts) if pred_contacts else float("nan")
    recall = len(inter) / len(native_contacts) if native_contacts else float("nan")

    if pred_contacts and native_contacts and precision + recall > 0:
        f1 = 2 * precision * recall / (precision + recall)
    else:
        f1 = float("nan")

    return precision, recall, f1, len(native_contacts), len(pred_contacts)


def calc_shell_mutation(native_model, pred_model, positions):
    total = 0
    mutation = 0

    for pos in positions:
        nres = native_model["parsed"].get(pos)
        pres = pred_model["parsed"].get(pos)

        if nres is None or pres is None:
            continue

        total += 1

        if nres["aa"] != pres["aa"]:
            mutation += 1

    rate = mutation / total if total else float("nan")
    return mutation, total, rate


# ============================================================
# 鏂囦欢璇嗗埆
# ============================================================

def infer_target_from_filename(pdb_path):
    stem = pdb_path.stem.lower()

    for key in sorted(TARGET_ALIASES.keys(), key=len, reverse=True):
        pattern = r"(^|_)" + re.escape(key) + r"(_|$)"
        if re.search(pattern, stem):
            return TARGET_ALIASES[key]

    return None


def find_native_pdb(native_dir, target):
    candidates = [
        native_dir / f"{target}.pdb",
        native_dir / f"{target.upper()}.pdb",
        native_dir / f"{target.lower()}.pdb",
    ]

    if "_" in target:
        base = target.split("_")[0]
        candidates.extend([
            native_dir / f"{base}.pdb",
            native_dir / f"{base.upper()}.pdb",
            native_dir / f"{base.lower()}.pdb",
        ])

    lower_map = {p.name.lower(): p for p in native_dir.glob("*.pdb")}

    for p in candidates:
        if p.exists():
            return p

        hit = lower_map.get(p.name.lower())
        if hit is not None:
            return hit

    return None


def parse_rank_score(pdb_path):
    stem = pdb_path.stem

    rank = ""
    score = ""

    m = re.search(r"rank(\d+)", stem, flags=re.IGNORECASE)
    if m:
        rank = int(m.group(1))

    m = re.search(r"score([0-9]+(?:\.[0-9]+)?)", stem, flags=re.IGNORECASE)
    if m:
        score = float(m.group(1))

    return rank, score


# ============================================================
# main
# ============================================================

def main():
    pred_files = sorted(PRED_DIR.rglob("*.pdb") if RECURSIVE else PRED_DIR.glob("*.pdb"))

    if not pred_files:
        raise FileNotFoundError(f"No pdb files found in {PRED_DIR}")

    print("========== compute_microenv_metrics ==========")
    print(f"[INFO] DESIGN_SET = {DESIGN_SET}")
    print(f"[INFO] NATIVE_DIR = {NATIVE_DIR}")
    print(f"[INFO] PRED_DIR   = {PRED_DIR}")
    print(f"[INFO] OUT_CSV    = {OUT_CSV}")
    print(f"[INFO] RECURSIVE  = {RECURSIVE}")
    print(f"[INFO] pred files = {len(pred_files)}")

    native_cache = {}
    rows = []

    ok = 0
    skipped = 0
    failed = 0

    for pred_pdb in pred_files:
        target = infer_target_from_filename(pred_pdb)

        if target is None:
            print(f"[SKIP] cannot infer target: {pred_pdb.name}")
            skipped += 1
            continue

        try:
            if target not in native_cache:
                native_pdb = find_native_pdb(NATIVE_DIR, target)

                if native_pdb is None:
                    raise FileNotFoundError(f"Cannot find native PDB for target={target}")

                native_model = build_native_model(native_pdb, target, chain_id=CHAIN_ID)

                shell5 = get_shell_positions(native_model, SHELL5_CUTOFF)
                shell8 = get_shell_positions(native_model, SHELL8_CUTOFF)

                native_cache[target] = {
                    "native_pdb": native_pdb,
                    "native_model": native_model,
                    "shell5": shell5,
                    "shell8": shell8,
                }

                print(
                    f"[NATIVE] {target}: "
                    f"native={native_pdb.name}, "
                    f"seq_len={len(native_model['sequence'])}, "
                    f"chromo={native_model['chromo_labels']}, "
                    f"chromo_pos={native_model['chromo_positions']}, "
                    f"shell5={len(shell5)}, shell8={len(shell8)}"
                )

            native_item = native_cache[target]
            native_model = native_item["native_model"]
            native_pdb = native_item["native_pdb"]
            shell5 = native_item["shell5"]
            shell8 = native_item["shell8"]

            pred_model = build_pred_model(pred_pdb, chain_id=CHAIN_ID)

            R, t, align_positions = get_alignment(native_model, pred_model)

            rmsd_shell5_ca, shell5_ca_count = calc_shell_ca_rmsd(
                native_model, pred_model, shell5, R, t
            )

            rmsd_shell8_ca, shell8_ca_count = calc_shell_ca_rmsd(
                native_model, pred_model, shell8, R, t
            )

            chromo_clash_count_2A, min_chromo_protein_dist = calc_chromo_clash(
                native_model, pred_model, R, t, clash_cutoff=CLASH_CUTOFF
            )

            contact_precision, contact_recall, contact_f1, native_contact_n, pred_contact_n = calc_contact_f1(
                native_model, pred_model, R, t, contact_cutoff=CONTACT_CUTOFF
            )

            shell5_mutation_count, shell5_mutation_total, shell5_mutation_rate = calc_shell_mutation(
                native_model, pred_model, shell5
            )

            rank, score_from_name = parse_rank_score(pred_pdb)

            row = {
                "design_set": DESIGN_SET,
                "target": target,
                "pred_file": str(pred_pdb),
                "native_file": str(native_pdb),
                "rank": rank,
                "score_from_filename": score_from_name,

                "native_seq_len": len(native_model["sequence"]),
                "pred_seq_len": len(pred_model["sequence"]),
                "alignment_ca_count": len(align_positions),

                "chromo_labels": ";".join(native_model["chromo_labels"]),
                "chromo_positions": " ".join(map(str, native_model["chromo_positions"])),

                "shell5_count": len(shell5),
                "shell5_positions": " ".join(map(str, shell5)),
                "shell8_count": len(shell8),
                "shell8_positions": " ".join(map(str, shell8)),

                "rmsd_shell5_ca": rmsd_shell5_ca,
                "rmsd_shell5_ca_count": shell5_ca_count,
                "rmsd_shell8_ca": rmsd_shell8_ca,
                "rmsd_shell8_ca_count": shell8_ca_count,

                "chromo_clash_count_2A": chromo_clash_count_2A,
                "min_chromo_protein_dist": min_chromo_protein_dist,

                "native_chromo_contact_count": native_contact_n,
                "pred_chromo_contact_count": pred_contact_n,
                "chromo_contact_precision": contact_precision,
                "chromo_contact_recall": contact_recall,
                "chromo_contact_F1": contact_f1,

                "shell5_mutation_count": shell5_mutation_count,
                "shell5_mutation_total": shell5_mutation_total,
                "shell5_mutation_rate": shell5_mutation_rate,
            }

            rows.append(row)
            ok += 1
            print(f"[OK] {pred_pdb.name} -> {target}")

        except Exception as e:
            failed += 1
            print(f"[FAIL] {pred_pdb.name}: {e}")

    if not rows:
        raise RuntimeError("No successful results. Please check file names and paths.")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows[0].keys())

    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow({k: clean_float(v) for k, v in row.items()})

    print()
    print("========== DONE ==========")
    print(f"[INFO] OK      : {ok}")
    print(f"[INFO] SKIPPED : {skipped}")
    print(f"[INFO] FAILED  : {failed}")
    print(f"[OK] output CSV: {OUT_CSV}")


if __name__ == "__main__":
    main()
