#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build S15 source CSVs from generated threshold-sensitivity eval summaries."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


METHOD_RE = re.compile(r"^(?P<kind>Bayes_fixedGraph|BayesNoPrior)_inferMinNaive(?P<threshold>\d+)$")
METHOD_LABELS = {
    "Bayes_fixedGraph": "PanTCR_fixed_graph",
    "BayesNoPrior": "PanTCR_NP_no_prior",
}
METRICS = ["Precision", "Recall", "F1-score"]


def formatted_metric(row: pd.Series) -> str:
    if isinstance(row.get("formatted"), str) and row["formatted"]:
        return row["formatted"]
    mean = float(row["mean_of_population_means"])
    std = float(row["std_of_population_means"])
    return f"{mean:.3f}±{std:.3f}"


def parse_method(method: str) -> tuple[str, int] | None:
    m = METHOD_RE.match(str(method))
    if not m:
        return None
    return METHOD_LABELS[m.group("kind")], int(m.group("threshold"))


def build_sources(eval_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    means = {}
    for _, row in eval_summary.iterrows():
        parsed = parse_method(str(row["method"]))
        if parsed is None:
            continue
        label, threshold = parsed
        metric = str(row["metric"])
        if metric not in METRICS:
            continue
        key = (str(row["expr"]), str(row["gene"]), label, threshold)
        means[(key, metric)] = float(row["mean_of_population_means"])
        rows.append(
            {
                "expr": str(row["expr"]),
                "gene": str(row["gene"]),
                "method_label": label,
                "inference_min_naive": threshold,
                "metric": metric,
                "formatted": formatted_metric(row),
            }
        )

    long = pd.DataFrame(rows)
    if long.empty:
        raise ValueError("No threshold-sensitivity methods found in population_mean_eval_summary.csv")
    metrics = (
        long.pivot_table(
            index=["gene", "method_label", "inference_min_naive"],
            columns="metric",
            values="formatted",
            aggfunc="first",
        )
        .reset_index()
    )
    for metric in METRICS:
        if metric not in metrics.columns:
            metrics[metric] = ""
    metrics = metrics[["gene", "method_label", "inference_min_naive", *METRICS]]
    metrics = metrics.sort_values(["gene", "method_label", "inference_min_naive"]).reset_index(drop=True)

    delta_rows = []
    keys = sorted({key for key, metric in means if metric == "F1-score"})
    for expr, gene, _label, threshold in keys:
        pan = means.get(((expr, gene, "PanTCR_fixed_graph", threshold), "F1-score"))
        np = means.get(((expr, gene, "PanTCR_NP_no_prior", threshold), "F1-score"))
        if pan is None or np is None:
            continue
        delta_rows.append(
            {
                "expr": expr,
                "inference_min_naive": threshold,
                "gene": gene,
                "PanTCR_NP_no_prior": np,
                "PanTCR_fixed_graph": pan,
                "PanTCR_minus_NP_F1": pan - np,
            }
        )
    delta = pd.DataFrame(delta_rows).sort_values(["expr", "inference_min_naive", "gene"]).reset_index(drop=True)
    return metrics, delta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval-summary-dir", required=True, help="Directory containing population_mean_eval_summary.csv")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    eval_summary_dir = Path(args.eval_summary_dir)
    out_dir = Path(args.out_dir)
    summary_csv = eval_summary_dir / "population_mean_eval_summary.csv"
    if not summary_csv.exists():
        raise FileNotFoundError(f"Missing evaluation summary CSV: {summary_csv}")
    metrics, delta = build_sources(pd.read_csv(summary_csv))
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out_dir / "threshold_precision_recall_f1_by_gene.csv", index=False, encoding="utf-8")
    delta.to_csv(out_dir / "threshold_f1_delta_summary.csv", index=False, encoding="utf-8")
    print(f"Wrote threshold source CSVs to {out_dir}")


if __name__ == "__main__":
    main()
