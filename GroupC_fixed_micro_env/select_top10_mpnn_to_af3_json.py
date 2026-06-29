#!/usr/bin/env python3
import json
import re
from pathlib import Path

# =========================
# 直接在这里改路径
# =========================
INPUT_DIR = Path("/work/home/mayongze/ProteinMPNN_project/Graduation_Design/GroupC_fixed_micro_env/output_GroupC/seqs")
OUTPUT_DIR = Path("/work/home/mayongze/ProteinMPNN_project/Graduation_Design/GroupC_fixed_micro_env/output_GroupC/af3_json_top10")

TOP_N = 10
MODEL_SEED = 1

# 你的 fa 中有 XXX；如果不想让 X 进入 AF3，可设为 False
ALLOW_X = True

score_re = re.compile(r"score=([0-9.]+)")


def read_fasta(fa):
    records = []
    header = None
    seq = []

    with open(fa) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

    if header is not None:
        records.append((header, "".join(seq)))

    return records


def get_score(header):
    m = score_re.search(header)
    if not m:
        return None
    return float(m.group(1))


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def make_af3_json(name, seq):
    return {
        "name": name,
        "modelSeeds": [MODEL_SEED],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 1
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fa_files = sorted(list(INPUT_DIR.glob("*.fa")) + list(INPUT_DIR.glob("*.fasta")))

    if not fa_files:
        raise FileNotFoundError(f"No fa/fasta files found in {INPUT_DIR}")

    list_file = OUTPUT_DIR / "af3_json_files.list"
    manifest_file = OUTPUT_DIR / "selected_top10_manifest.tsv"

    all_json_paths = []
    manifest_rows = []

    for fa in fa_files:
        records = read_fasta(fa)

        if len(records) <= 1:
            print(f"[SKIP] {fa.name}: only {len(records)} record")
            continue

        # 跳过第一条 WT/input
        design_records = records[1:]

        candidates = []
        for idx, (header, seq) in enumerate(design_records, start=2):
            score = get_score(header)
            if score is None:
                print(f"[SKIP] {fa.name} record {idx}: no score")
                continue

            seq = seq.upper().replace(" ", "")

            if not ALLOW_X and "X" in seq:
                print(f"[SKIP] {fa.name} record {idx}: contains X")
                continue

            candidates.append({
                "idx": idx,
                "header": header,
                "seq": seq,
                "score": score
            })

        candidates.sort(key=lambda x: x["score"])
        selected = candidates[:TOP_N]

        out_subdir = OUTPUT_DIR / safe_name(fa.stem)
        out_subdir.mkdir(parents=True, exist_ok=True)

        for rank, item in enumerate(selected, start=1):
            name = f"{safe_name(fa.stem)}_rank{rank:02d}_score{item['score']:.4f}"
            json_path = out_subdir / f"{name}.json"

            af3_json = make_af3_json(name, item["seq"])

            with open(json_path, "w") as f:
                json.dump(af3_json, f, indent=2)

            all_json_paths.append(json_path)

            manifest_rows.append(
                f"{fa}\t{json_path}\t{rank}\t{item['idx']}\t{item['score']}\t{item['header']}"
            )

        print(f"[OK] {fa.name}: selected {len(selected)} sequences")

    with open(list_file, "w") as f:
        for p in all_json_paths:
            f.write(str(p) + "\n")

    with open(manifest_file, "w") as f:
        f.write("source_fa\tjson_path\trank\trecord_index\tscore\theader\n")
        for row in manifest_rows:
            f.write(row + "\n")

    print("\n========== DONE ==========")
    print(f"Input dir:  {INPUT_DIR}")
    print(f"Output dir: {OUTPUT_DIR}")
    print(f"JSON count: {len(all_json_paths)}")
    print(f"List file:  {list_file}")
    print(f"Manifest:   {manifest_file}")


if __name__ == "__main__":
    main()