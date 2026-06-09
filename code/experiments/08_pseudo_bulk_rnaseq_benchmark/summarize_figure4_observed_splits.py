#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the final two-table Figure 4 non-default allele evidence audit.

Supplementary Table S18 uses unfiltered MiXCR-derived observed-region coverage
for method comparison. Supplementary Table S19 uses PanTCR-retained evidence for
the pseudo-bulk recovered allele mechanism audit.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd

CODE_EXPERIMENT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = CODE_EXPERIMENT_DIR.parents[2]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / CODE_EXPERIMENT_DIR.name

sys.path.append(str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))
from common_pantcr_io import (  # noqa: E402
    defining_positions,
    evidence_coverage_by_gene,
    load_default_refs,
    load_trb_index,
    trim_sequence,
)


SC_ARTIFACTS = {
    ("SC-2", "TRBV11-1", "TRBV11-1*01"),
    ("SC-2", "TRBV6-9", "TRBV6-9*01"),
    ("SC-2", "TRBV7-3", "TRBV7-3*01"),
    ("SC-3", "TRBV7-4", "TRBV7-4*01"),
    ("SC-4", "TRBV7-4", "TRBV7-4*01"),
    ("SC-6", "TRBV7-4", "TRBV7-4*01"),
    ("SC-8", "TRBV5-3", "TRBV5-3*01"),
}


def load_trimmed_default_refs(ref_root: Path) -> Dict[str, str]:
    index_df = load_trb_index(ref_root / "TRB_index.csv")
    refs = load_default_refs(ref_root / "IMGT_TRB_default.csv")
    return {gene: trim_sequence(seq, gene, "V", index_df) for gene, seq in refs.items()}


def observed_status(def_pos: List[int], covered_positions: Set[int]) -> Tuple[str, int]:
    n_covered = sum(1 for p in def_pos if p in covered_positions)
    if n_covered == len(def_pos):
        return "All defining sites observed", n_covered
    if n_covered == 0:
        return "No defining site observed", n_covered
    return "Partially observed defining sites", n_covered


def build_truth_key(row: pd.Series, source: str) -> str:
    if source == "semi":
        fields = [
            row["dataset"],
            row["sample_id"],
            row["gene"],
            row["truth_allele"],
            row["truth_seq"],
        ]
    else:
        fields = [
            row["SC_ID"],
            row["DatasetID"],
            row["gene"],
            row["truth_allele"],
            row["truth_seq"],
        ]
    return "||".join(str(x) for x in fields)


def add_observed_strata(
    df: pd.DataFrame,
    ref_root: Path,
    source: str,
) -> pd.DataFrame:
    default_refs = load_trimmed_default_refs(ref_root)
    coverage_cache: Dict[str, Dict[str, Set[int]]] = {}
    rows = []
    for _, row in df.iterrows():
        gene = str(row["gene"])
        truth_seq = str(row["truth_seq"]).strip().upper()
        ref_seq = default_refs.get(gene, "")
        if not ref_seq:
            raise ValueError(f"Missing default reference for {gene}")
        def_pos = defining_positions(truth_seq, ref_seq)
        if not def_pos:
            continue
        if source == "sc" and (str(row["SC_ID"]), gene, str(row["truth_allele"])) in SC_ARTIFACTS:
            continue
        evidence_file = Path(str(row["evidence_file"]))
        if not evidence_file.is_absolute():
            evidence_file = Path.cwd() / evidence_file
        cache_key = str(evidence_file)
        if cache_key not in coverage_cache:
            # Method-comparison table uses observed regions from MiXCR-derived
            # rearrangement evidence before PanTCR-specific support filtering.
            coverage_cache[cache_key] = evidence_coverage_by_gene(evidence_file, min_naive=None)
        obs_status, n_obs = observed_status(def_pos, coverage_cache[cache_key].get(gene, set()))
        out = row.to_dict()
        out["truth_key"] = build_truth_key(row, source)
        out["No. of default-relative defining sites"] = len(def_pos)
        out["No. of observed defining sites"] = n_obs
        out["Observed defining-site stratum"] = obs_status
        rows.append(out)
    return pd.DataFrame(rows)


def summarize_methods(stratified: pd.DataFrame, method_map: Dict[str, str], dataset_col: str, benchmark: str, dataset_order: List[str]) -> pd.DataFrame:
    base = stratified[stratified["method"].isin(method_map.keys())].copy()
    truth_meta = (
        base.sort_values(["truth_key", "method"])
        .drop_duplicates("truth_key")[
            [
                "truth_key",
                dataset_col,
                "Observed defining-site stratum",
                "No. of default-relative defining sites",
                "No. of observed defining sites",
            ]
        ]
        .copy()
    )
    rows = []
    strata_order = [
        "All defining sites observed",
        "Partially observed defining sites",
        "No defining site observed",
    ]
    for dataset in dataset_order:
        d_truth = truth_meta[truth_meta[dataset_col] == dataset]
        for stratum in strata_order:
            keys = set(d_truth[d_truth["Observed defining-site stratum"] == stratum]["truth_key"])
            if not keys:
                continue
            row = {
                "Benchmark": benchmark,
                "Dataset": dataset,
                "Observed Defining-site Stratum": stratum,
                "No. of Non-default Truth Alleles": len(keys),
            }
            for raw_method, display in method_map.items():
                m = base[(base["method"] == raw_method) & (base["truth_key"].isin(keys))]
                row[f"{display} Exact Recovered"] = int((m["status"] == "exact_tp").sum())
            rows.append(row)
    out = pd.DataFrame(rows)
    return out


def add_overall(rows: pd.DataFrame, benchmark: str, datasets: List[str]) -> pd.DataFrame:
    sub = rows[(rows["Benchmark"] == benchmark) & (rows["Dataset"].isin(datasets))].copy()
    if sub.empty:
        return rows
    numeric_cols = [c for c in sub.columns if c not in {"Benchmark", "Dataset", "Observed Defining-site Stratum"}]
    overall = (
        sub.groupby("Observed Defining-site Stratum", as_index=False)[numeric_cols]
        .sum()
        .assign(Benchmark=benchmark, Dataset="Overall")
    )
    overall = overall[rows.columns]
    return pd.concat([rows, overall], ignore_index=True)


def build_method_comparison(args: argparse.Namespace) -> pd.DataFrame:
    semi = pd.read_csv(args.semi_truth_csv)
    semi = semi[(semi["gene_type"] == "V") & (semi["method"].isin(["Bayes", "MiXCR", "FindAlleles"]))].copy()
    semi = add_observed_strata(semi, args.semi_ref_root, source="semi")
    semi_summary = summarize_methods(
        semi,
        {"Bayes": "PanTCR", "MiXCR": "MiXCR-default", "FindAlleles": "FindAlleles"},
        "dataset",
        "Semi-synthetic AIRR-seq",
        ["AIRR-1", "AIRR-2"],
    )

    sc = pd.read_csv(args.sc_truth_csv)
    sc = sc[(sc["gene_type"] == "V") & (sc["method"].isin(["PanTCR.semi", "MiXCR", "FindAlleles"]))].copy()
    sc = add_observed_strata(sc, args.sc_ref_root, source="sc")
    sc_summary = summarize_methods(
        sc,
        {"PanTCR.semi": "PanTCR", "MiXCR": "MiXCR-default", "FindAlleles": "FindAlleles"},
        "SC_ID",
        "Pseudo-bulk RNA-seq",
        [f"SC-{i}" for i in range(1, 9)],
    )
    sc_summary = add_overall(sc_summary, "Pseudo-bulk RNA-seq", [f"SC-{i}" for i in range(1, 9)])
    return pd.concat([semi_summary, sc_summary], ignore_index=True)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--semi-truth-csv",
        type=Path,
        default=PACKAGE_ROOT
        / "experiments"
        / "07_semi_synthetic_airr_benchmark"
        / "generated"
        / "evidence_analysis"
        / "per_truth_call_status.csv",
    )
    ap.add_argument("--semi-ref-root", type=Path, default=PACKAGE_ROOT / "data" / "ref")
    ap.add_argument(
        "--sc-truth-csv",
        type=Path,
        default=EXPERIMENT_DIR / "generated" / "evidence_analysis" / "per_truth_call_status.csv",
    )
    ap.add_argument("--sc-ref-root", type=Path, default=PACKAGE_ROOT / "data" / "ref")
    ap.add_argument(
        "--examples-csv",
        type=Path,
        default=EXPERIMENT_DIR
        / "generated"
        / "fig4e_nonreference_variant_audit"
        / "fig4e_recovered_nonreference_alleles_for_table.csv",
    )
    ap.add_argument("--out-dir", type=Path, default=EXPERIMENT_DIR / "generated" / "final_two_table_audit")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    method_comparison = build_method_comparison(args)
    recovered_examples = pd.read_csv(args.examples_csv)
    method_comparison.to_csv(out_dir / "method_comparison_by_observed_defining_site.csv", index=False, encoding="utf-8-sig")
    recovered_examples.to_csv(out_dir / "pantcr_recovered_scbulk_examples.csv", index=False, encoding="utf-8-sig")
    with (out_dir / "workbook_data.json").open("w", encoding="utf-8") as fh:
        json.dump(
            {
                "method_comparison": method_comparison.fillna("").to_dict(orient="records"),
                "recovered_examples": recovered_examples.fillna("").to_dict(orient="records"),
            },
            fh,
            ensure_ascii=False,
        )
    print(f"Wrote final two-table audit data to {out_dir}")
    print(method_comparison.to_string(index=False))


if __name__ == "__main__":
    main()
