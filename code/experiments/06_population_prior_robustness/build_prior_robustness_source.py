#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the S14 source CSV from generated population-prior eval summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_EXPR_LABELS = {"expr_ScenarioA": "Scenario A", "expr_ScenarioB": "Scenario B"}
MODE_ORDER = [
    ("No graph prior", ("BayesNoPrior_no_prior", "BayesNoPrior_review")),
    ("EUR-mismatched prior", ("Bayes_priorMismatch", "Bayes_priorMismatch_AFR_to_EUR")),
    ("Global prior", ("Bayes_priorGlobal", "Bayes_priorGlobal_AFR")),
    ("AFR-matched prior", ("Bayes_priorMatched", "Bayes_priorMatched_AFR")),
]


def parse_expr_labels(text: str) -> dict[str, str]:
    labels = dict(DEFAULT_EXPR_LABELS)
    for item in text.replace(",", " ").split():
        if ":" not in item:
            continue
        expr, label = item.split(":", 1)
        labels[expr] = label.replace("_", " ")
    return labels


def match_method(method: str, aliases: tuple[str, ...]) -> bool:
    return any(method == alias or method.startswith(alias + "_") for alias in aliases)


def metric_fmt(values: pd.Series) -> str:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return "nan±nan"
    return f"{numeric.mean():.3f}±{numeric.std(ddof=1):.3f}"


def build_source(per_sample: pd.DataFrame, expr_labels: dict[str, str], gene: str, population: str) -> pd.DataFrame:
    df = per_sample.copy()
    df = df[(df["gene"].astype(str) == gene) & (df["population"].astype(str) == population)].copy()
    if df.empty:
        raise ValueError(f"No per-sample eval rows found for gene={gene}, population={population}")
    rows = []
    for expr, setting in expr_labels.items():
        expr_df = df[df["expr"].astype(str) == expr].copy()
        if expr_df.empty:
            continue
        for mode, aliases in MODE_ORDER:
            sub = expr_df[expr_df["method"].astype(str).map(lambda m: match_method(m, aliases))].copy()
            if sub.empty:
                continue
            for col in ["NO. of Truth Alleles", "NO. of Predicted Alleles", "TP", "FP", "FN"]:
                sub[col] = pd.to_numeric(sub[col], errors="coerce").fillna(0).astype(int)
            rows.append(
                {
                    "Setting": setting,
                    "Mode": mode,
                    "NO. of Truth Alleles": int(sub["NO. of Truth Alleles"].sum()),
                    "NO. of Predicted Alleles": int(sub["NO. of Predicted Alleles"].sum()),
                    "TP": int(sub["TP"].sum()),
                    "FP": int(sub["FP"].sum()),
                    "FN": int(sub["FN"].sum()),
                    "Precision": metric_fmt(sub["Precision"]),
                    "Recall": metric_fmt(sub["Recall"]),
                    "F1-score": metric_fmt(sub["F1-score"]),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("No prior-robustness method rows matched the requested aliases.")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--eval-summary-dir", required=True, help="Directory containing per_sample_eval_metrics.csv")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--expr-labels", default="expr_ScenarioA:Scenario_A expr_ScenarioB:Scenario_B")
    ap.add_argument("--gene", default="V")
    ap.add_argument("--population", default="AFR")
    args = ap.parse_args()

    per_sample_path = Path(args.eval_summary_dir) / "per_sample_eval_metrics.csv"
    if not per_sample_path.exists():
        raise FileNotFoundError(f"Missing per-sample evaluation metrics: {per_sample_path}")
    source = build_source(
        pd.read_csv(per_sample_path),
        parse_expr_labels(args.expr_labels),
        gene=args.gene,
        population=args.population,
    )
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    source.to_csv(out_csv, index=False, encoding="utf-8")
    print(f"Wrote prior robustness source CSV to {out_csv}")


if __name__ == "__main__":
    main()
