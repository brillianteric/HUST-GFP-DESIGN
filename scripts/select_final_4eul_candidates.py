#!/usr/bin/env python3
"""
Reproduce the final 4EUL GFP candidate selection.

This script implements the conservative filtering logic used to choose six
4EUL-derived designs for the 2026 Protein Design Synbio Challenge.

Inputs
------
1. structure_metrics_with_plddt_pae.csv
2. microenv_metrics.csv

Outputs
-------
1. selected_4EUL_candidate_ids.csv
   Final selected candidates, sorted by ProteinMPNN score.
2. selected_4EUL_audit.csv
   All target candidates with pass/fail status for each filter.

Example
-------
python select_final_4eul_candidates.py \
  --structure structure_metrics_with_plddt_pae.csv \
  --microenv microenv_metrics.csv \
  --out selected_4EUL_candidate_ids.csv \
  --audit-out selected_4EUL_audit.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import Dict, Iterable

import pandas as pd


DEFAULT_THRESHOLDS = {
    # Global structure preservation
    "ca_rmsd_max": 0.325,
    "tm_score_min": 0.9965,
    "plddt_min": 94.8,
    "pae_max": 2.17,
    # Chromophore microenvironment preservation
    "shell5_rmsd_max": 0.325,
    "shell8_rmsd_max": 0.29,
    "contact_f1_min": 0.84,
    "clash_2a_max": 5,
    "shell5_mutation_max": 0,
}


def parse_rank(value: object) -> int:
    """Convert values such as 'rank01' or 1 into integer rank numbers."""
    match = re.search(r"(\d+)", str(value))
    if not match:
        raise ValueError(f"Cannot parse rank from value: {value!r}")
    return int(match.group(1))


def require_columns(df: pd.DataFrame, columns: Iterable[str], table_name: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {missing}")


def normalize_bool(series: pd.Series) -> pd.Series:
    """Handle bool columns read either as bools or strings."""
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def make_candidate_id(row: pd.Series, candidate_prefix: str) -> str:
    return f"{candidate_prefix}_rank{int(row['rank_num']):02d}_score{float(row['mpnn_score_from_name']):.4f}"


def build_filter_columns(df: pd.DataFrame, thresholds: Dict[str, float]) -> pd.DataFrame:
    """Add one boolean column per filter and a final all-pass column."""
    df = df.copy()

    filter_exprs = {
        "pass_length_match_after_exclusion": normalize_bool(df["length_match_after_exclusion"]),
        "pass_shell5_mutation_count": df["shell5_mutation_count"].le(thresholds["shell5_mutation_max"]),
        "pass_ca_rmsd": df["ca_rmsd_excluding_chromophore"].le(thresholds["ca_rmsd_max"]),
        "pass_tm_score": df["ca_tm_score_approx_excluding_chromophore"].ge(thresholds["tm_score_min"]),
        "pass_plddt": df["atom_plddt_mean"].ge(thresholds["plddt_min"]),
        "pass_pae": df["pae_mean"].le(thresholds["pae_max"]),
        "pass_shell5_rmsd": df["rmsd_shell5_ca"].le(thresholds["shell5_rmsd_max"]),
        "pass_shell8_rmsd": df["rmsd_shell8_ca"].le(thresholds["shell8_rmsd_max"]),
        "pass_contact_f1": df["chromo_contact_F1"].ge(thresholds["contact_f1_min"]),
        "pass_clash_2a": df["chromo_clash_count_2A"].le(thresholds["clash_2a_max"]),
    }

    pass_cols = []
    for col, values in filter_exprs.items():
        df[col] = values.fillna(False).astype(bool)
        pass_cols.append(col)

    df["passes_all_filters"] = df[pass_cols].all(axis=1)

    def failed_reasons(row: pd.Series) -> str:
        failed = [col.replace("pass_", "") for col in pass_cols if not bool(row[col])]
        return ";".join(failed)

    df["failed_filters"] = df.apply(failed_reasons, axis=1)
    return df


def load_and_merge(structure_path: Path, microenv_path: Path, target: str) -> pd.DataFrame:
    structure = pd.read_csv(structure_path)
    microenv = pd.read_csv(microenv_path)

    require_columns(
        structure,
        [
            "pdb_id",
            "rank",
            "mpnn_score_from_name",
            "wt_mpnn_score",
            "ca_rmsd_excluding_chromophore",
            "ca_tm_score_approx_excluding_chromophore",
            "design_ca_len",
            "length_match_after_exclusion",
            "design_file",
            "atom_plddt_mean",
            "pae_mean",
        ],
        "structure metrics CSV",
    )
    require_columns(
        microenv,
        [
            "target",
            "rank",
            "rmsd_shell5_ca",
            "rmsd_shell8_ca",
            "chromo_clash_count_2A",
            "chromo_contact_F1",
            "shell5_mutation_count",
        ],
        "microenvironment metrics CSV",
    )

    structure = structure.copy()
    structure["target"] = structure["pdb_id"].astype(str).str.upper()
    structure["rank_num"] = structure["rank"].apply(parse_rank)

    microenv = microenv.copy()
    microenv["target"] = microenv["target"].astype(str).str.upper()
    microenv["rank_num"] = microenv["rank"].apply(parse_rank)

    merged = structure.merge(
        microenv,
        on=["target", "rank_num"],
        suffixes=("_structure", "_microenv"),
        how="inner",
        validate="one_to_one",
    )
    merged = merged[merged["target"].eq(target.upper())].copy()
    if merged.empty:
        raise ValueError(f"No merged candidates found for target={target!r}")
    return merged


def build_selected_table(df: pd.DataFrame, candidate_prefix: str, top_n: int) -> pd.DataFrame:
    passed = df[df["passes_all_filters"]].copy()
    passed = passed.sort_values(["mpnn_score_from_name", "rank_num"], ascending=[True, True]).head(top_n)

    selected = pd.DataFrame(
        {
            "Seq_ID": range(1, len(passed) + 1),
            "Candidate_ID": passed.apply(lambda row: make_candidate_id(row, candidate_prefix), axis=1),
            "Target": passed["target"],
            "Rank": passed["rank_num"].astype(int),
            "Predicted_Length_without_added_M": passed["design_ca_len"].astype(int),
            "MPNN_score": passed["mpnn_score_from_name"],
            "WT_MPNN_score": passed["wt_mpnn_score"],
            "CA_RMSD_excluding_chromophore": passed["ca_rmsd_excluding_chromophore"],
            "TM_score_approx_excluding_chromophore": passed["ca_tm_score_approx_excluding_chromophore"],
            "pLDDT_mean": passed["atom_plddt_mean"],
            "PAE_mean": passed["pae_mean"],
            "Shell5_RMSD": passed["rmsd_shell5_ca"],
            "Shell8_RMSD": passed["rmsd_shell8_ca"],
            "Contact_F1": passed["chromo_contact_F1"],
            "Clash_2A": passed["chromo_clash_count_2A"],
            "Shell5_mutation_count": passed["shell5_mutation_count"],
            "Design_File": passed["design_file"],
        }
    )
    return selected


def build_audit_table(df: pd.DataFrame, candidate_prefix: str) -> pd.DataFrame:
    audit = df.copy().sort_values("rank_num")
    audit.insert(0, "Candidate_ID", audit.apply(lambda row: make_candidate_id(row, candidate_prefix), axis=1))

    keep_cols = [
        "Candidate_ID",
        "target",
        "rank_num",
        "mpnn_score_from_name",
        "wt_mpnn_score",
        "design_ca_len",
        "length_match_after_exclusion",
        "ca_rmsd_excluding_chromophore",
        "ca_tm_score_approx_excluding_chromophore",
        "atom_plddt_mean",
        "pae_mean",
        "rmsd_shell5_ca",
        "rmsd_shell8_ca",
        "chromo_contact_F1",
        "chromo_clash_count_2A",
        "shell5_mutation_count",
        "passes_all_filters",
        "failed_filters",
        "design_file",
    ]
    keep_cols += [c for c in audit.columns if c.startswith("pass_")]
    return audit[keep_cols]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Reproduce final 4EUL GFP candidate selection from metrics CSV files."
    )
    parser.add_argument("--structure", required=True, type=Path, help="4EUL structure metrics CSV")
    parser.add_argument("--microenv", required=True, type=Path, help="4EUL microenvironment metrics CSV")
    parser.add_argument("--target", default="4EUL", help="Target PDB ID to select from; default: 4EUL")
    parser.add_argument("--candidate-prefix", default="4EUL", help="Prefix used for Candidate_ID; default: 4EUL")
    parser.add_argument("--top-n", default=6, type=int, help="Number of final candidates to output; default: 6")
    parser.add_argument("--out", required=True, type=Path, help="Output selected candidate CSV")
    parser.add_argument("--audit-out", required=True, type=Path, help="Output audit CSV with pass/fail details")

    # Thresholds are exposed as CLI arguments so the selection rule is transparent.
    for name, default in DEFAULT_THRESHOLDS.items():
        parser.add_argument(f"--{name.replace('_', '-')}", type=float, default=default)

    args = parser.parse_args()
    thresholds = {name: getattr(args, name) for name in DEFAULT_THRESHOLDS}

    try:
        df = load_and_merge(args.structure, args.microenv, args.target)
        df = build_filter_columns(df, thresholds)
        selected = build_selected_table(df, args.candidate_prefix, args.top_n)
        audit = build_audit_table(df, args.candidate_prefix)
    except Exception as exc:  # noqa: BLE001 - command-line script should report cleanly
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.audit_out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out, index=False)
    audit.to_csv(args.audit_out, index=False)

    print(f"Target: {args.target.upper()}")
    print(f"Total candidates considered: {len(df)}")
    print(f"Candidates passing all filters: {int(df['passes_all_filters'].sum())}")
    print("Selected ranks:", ", ".join(f"rank{r:02d}" for r in selected["Rank"].astype(int)))
    print(f"Selected output: {args.out}")
    print(f"Audit output: {args.audit_out}")

    failed = audit[~audit["passes_all_filters"]]
    if not failed.empty:
        print("\nRejected candidates:")
        for _, row in failed.iterrows():
            print(f"  rank{int(row['rank_num']):02d}: {row['failed_filters']}")

    if len(selected) < args.top_n:
        print(
            f"WARNING: only {len(selected)} candidates passed all filters, fewer than requested top_n={args.top_n}.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
