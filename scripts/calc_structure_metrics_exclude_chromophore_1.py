#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import csv
import json
import math
import re
import numpy as np

from Bio.PDB import PDBParser, Superimposer


# =========================
# 鐢ㄦ埛鍙慨鏀瑰弬鏁?# =========================
DESIGN_PDB_DIR = Path("results/example_models")
TARGET_PDB_DIR = Path("data/raw")
CONFIDENCE_DIR = Path("outputs/af3_stage2")
WT_SCORE_FASTA = Path("data/proteinmpnn/4EUL_designs.fa")

OUT_CSV = Path("results/metrics/4EUL_structure_metrics_with_plddt_pae.csv")

# 璁捐缁撴瀯鏂囦欢鏍煎紡锛?# $DESIGN_PREFIX_{pdb_id}_{rank}_{score}_model.pdb
DESIGN_PREFIX = ""
DESIGN_SUFFIX = "_model.pdb"

# 浼樺厛姣旇緝鐨勯摼
PREFERRED_CHAIN_ID = "A"

# 杩欎簺浣嶇疆鏄璁＄粨鏋?/ MPNN 杈撳嚭搴忓垪涓殑 1-based 浣嶇疆
# 4EUL-only: exclude the chromophore-forming tripeptide from the design CA list.
EXCLUDE_DESIGN_POSITIONS = {
    "4EUL": [63, 64, 65],
}

# 濡傛灉鎺掗櫎鍚庝笁鑰呴暱搴︿粛涓嶄竴鑷达紝鏄惁鐢?min length 鎴柇缁х画璁＄畻
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
    鎻愬彇鎸囧畾 PDB 绗竴妯″瀷銆佹寚瀹氶摼涓殑 CA 鍧愭爣銆?    娉ㄦ剰锛?    - 璺宠繃姘淬€佺瀛愩€佸皬鍒嗗瓙锛?    - 灏ゅ叾璺宠繃 4EUL 涓殑閽欑瀛?CA锛?    - 涓嶄娇鐢?STANDARD_AA 鐧藉悕鍗曪紝閬垮厤鐮村潖鍘熸湁椤哄簭瀵归綈閫昏緫銆?    """
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

        # 鎺掗櫎姘淬€佺瀛愩€佸皬鍒嗗瓙锛屽挨鍏舵槸 4EUL 閲岀殑閽欑瀛?CA
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
    """
    浠?design CA 鍧愭爣涓垹闄ゆ寚瀹?1-based 搴忓垪浣嶇疆銆?    """
    if len(coords) == 0:
        return coords

    remove_set = set(remove_positions)
    keep_indices = [
        i for i in range(len(coords))
        if (i + 1) not in remove_set
    ]

    return coords[keep_indices]


def kabsch_superpose_rmsd(mobile_coords, target_coords):
    """
    瀵?mobile -> target 鍋?Kabsch superposition锛岃繑鍥?RMSD 鍜屽彉鎹㈠悗鐨?mobile 鍧愭爣銆?    """

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
    杩戜技 C伪 TM-score銆?
    娉ㄦ剰锛?    杩欐槸鍦ㄥ凡鏈?residue 椤哄簭瀵归綈鍜?Kabsch superposition 鍚庤绠楃殑 approximate TM-score锛?    涓嶆槸 US-align / TMscore 瀹樻柟宸ュ叿鐨勫畬鏁寸粨鏋勬瘮瀵圭粨鏋溿€?    """
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
    浠?fixed_{pdb_id}_{rank}_{score}_model.pdb 瑙ｆ瀽 pdb_id / rank / score銆?
    鏀寔灏忓啓璁捐鏂囦欢鍚嶅尮閰嶅ぇ鍐?target PDB锛?      fixed_4oqw_a_rank01_score0.7287_model.pdb
      -> 4EUL.pdb
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
):
    """
    璺嚎 B锛?    浠?design 缁撴瀯涓垹闄?chromophore-forming tripeptide 瀵瑰簲 CA锛?    鍐嶅拰鍘熷 target PDB 鐨?CA 鎸夐『搴忔瘮瀵广€?    """
    target = get_ca_coords(target_pdb, PREFERRED_CHAIN_ID)
    design = get_ca_coords(design_pdb, PREFERRED_CHAIN_ID)

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
        "confidence_json": str(conf_json),
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


def main():
    if not DESIGN_PDB_DIR.exists():
        raise FileNotFoundError(f"Design PDB dir not found: {DESIGN_PDB_DIR}")

    if not TARGET_PDB_DIR.exists():
        raise FileNotFoundError(f"Target PDB dir not found: {TARGET_PDB_DIR}")

    target_pdbs = {
        p.stem: p
        for p in sorted(TARGET_PDB_DIR.glob("*.pdb"))
    }

    target_ids = sorted(target_pdbs.keys())

    if not target_pdbs:
        raise FileNotFoundError(f"No target PDB files found in: {TARGET_PDB_DIR}")

    design_pdbs = sorted(DESIGN_PDB_DIR.glob("*.pdb"))

    print("========== INPUT ==========")
    print(f"[INFO] Design PDB dir : {DESIGN_PDB_DIR}")
    print(f"[INFO] Target PDB dir : {TARGET_PDB_DIR}")
    print(f"[INFO] Output CSV     : {OUT_CSV}")
    print(f"[INFO] Target PDBs    : {len(target_pdbs)}")
    print(f"[INFO] Design PDBs    : {len(design_pdbs)}")
    print(f"[INFO] Method         : exclude chromophore tripeptide from design CA, then sequential CA alignment")

    wt_mpnn_score = read_wt_mpnn_score(WT_SCORE_FASTA)
    confidence_by_rank = build_confidence_map(CONFIDENCE_DIR)

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
            )

            row = {
                "pdb_id": pdb_id,
                "wt_mpnn_score": wt_mpnn_score,
                "rank": parsed["rank"],
                "mpnn_score_from_name": parsed["mpnn_score_from_name"],
                **metrics,
                "design_file": str(design_pdb),
                "target_file": str(target_pdb),
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

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

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

    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("========== DONE ==========")
    print(f"[INFO] OK             : {n_ok}")
    print(f"[INFO] Failed         : {n_fail}")
    print(f"[INFO] Skip name      : {n_skip_name}")
    print(f"[INFO] Missing target : {n_missing_target}")
    print(f"[INFO] Output CSV     : {OUT_CSV}")


if __name__ == "__main__":
    main()
