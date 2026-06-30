#!/usr/bin/env python3
"""Select low-scoring ProteinMPNN designs and write AlphaFold 3 input JSON."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


SCORE_RE = re.compile(r"score=([0-9.]+)")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq: list[str] = []
    with path.open() as handle:
        for line in handle:
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


def score_from_header(header: str) -> float | None:
    match = SCORE_RE.search(header)
    return float(match.group(1)) if match else None


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def make_af3_json(name: str, sequence: str, seed: int) -> dict:
    return {
        "name": name,
        "modelSeeds": [seed],
        "sequences": [{"protein": {"id": "A", "sequence": sequence}}],
        "dialect": "alphafold3",
        "version": 1,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--model-seed", type=int, default=1)
    parser.add_argument("--allow-x", action="store_true")
    args = parser.parse_args()

    records = read_fasta(args.input_fasta)
    if len(records) <= 1:
        raise ValueError("Expected one WT record followed by ProteinMPNN design records")

    candidates = []
    for record_index, (header, sequence) in enumerate(records[1:], start=2):
        score = score_from_header(header)
        sequence = sequence.upper().replace(" ", "")
        if score is None:
            continue
        if not args.allow_x and "X" in sequence:
            continue
        candidates.append((score, record_index, header, sequence))

    candidates.sort(key=lambda item: item[0])
    selected = candidates[: args.top_n]

    target = safe_name(args.input_fasta.stem.split("_")[0])
    target_dir = args.output_dir / target
    target_dir.mkdir(parents=True, exist_ok=True)

    list_file = args.output_dir / "af3_json_files.list"
    manifest_file = args.output_dir / "selected_top10_manifest.tsv"

    json_paths = []
    manifest_rows = ["source_fa\tjson_path\trank\trecord_index\tscore\theader"]
    for rank, (score, record_index, header, sequence) in enumerate(selected, start=1):
        name = f"{target}_rank{rank:02d}_score{score:.4f}"
        json_path = target_dir / f"{name}.json"
        json_path.write_text(json.dumps(make_af3_json(name, sequence, args.model_seed), indent=2) + "\n")
        json_paths.append(json_path)
        manifest_rows.append(
            f"{args.input_fasta}\t{json_path}\t{rank}\t{record_index}\t{score}\t{header}"
        )

    list_file.write_text("".join(f"{path}\n" for path in json_paths))
    manifest_file.write_text("\n".join(manifest_rows) + "\n")

    print(f"Selected {len(json_paths)} designs")
    print(f"JSON list: {list_file}")
    print(f"Manifest: {manifest_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
