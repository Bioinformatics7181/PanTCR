#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize PanTCR evaluation CSV files under a generated results tree.

The input contract is a package-local or explicitly supplied results root with
files laid out as:

    results/eval/<expr>/<method>/<population>/<gene>/eval_*_summary.csv
    results/eval/<expr>/<method>/<population>/<gene>/eval_*_detailed.csv

This utility is the common post-processing step for downstream experiments that
reuse generated inference and evaluation outputs from an upstream simulation.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Iterable

import pandas as pd


METRICS = ["Precision", "Recall", "F1-score"]
SUMMARY_RE = re.compile(
    r"^eval_(?P<method>.+?)_(?P<pop>[A-Za-z]+)_(?P<seed>seed\d+)_(?P<gene>[VJ])_summary\.csv$"
)


def readable_path(path: Path) -> str | Path:
    if os.name != "nt":
        return path
    resolved = str(path.resolve())
    if resolved.startswith("\\\\?\\"):
        return resolved
    return "\\\\?\\" + resolved


def split_words(value: str | Iterable[str] | None, default: list[str]) -> list[str]:
    if value is None:
        return list(default)
    if isinstance(value, str):
        out = [x for x in value.replace(",", " ").split() if x]
        return out or list(default)
    return [str(x) for x in value]


def read_summary_metrics(summary_file: Path) -> dict[str, float]:
    table = pd.read_csv(readable_path(summary_file))
    if not {"Metric", "Value"}.issubset(table.columns):
        raise ValueError(f"Evaluation summary missing Metric/Value columns: {summary_file}")
    values = dict(zip(table["Metric"].astype(str), table["Value"]))
    missing = [metric for metric in METRICS if metric not in values]
    if missing:
        raise ValueError(f"Evaluation summary missing metrics {missing}: {summary_file}")
    return {metric: float(values.get(metric, "nan")) for metric in METRICS}


def read_detailed_counts(summary_file: Path) -> dict[str, int | str]:
    detailed = Path(str(summary_file).replace("_summary.csv", "_detailed.csv"))
    if not detailed.exists():
        return {
            "detailed_file": "",
            "TP": pd.NA,
            "FP": pd.NA,
            "FN": pd.NA,
            "NO. of Truth Alleles": pd.NA,
            "NO. of Predicted Alleles": pd.NA,
        }
    table = pd.read_csv(readable_path(detailed))
    if "Status" not in table.columns:
        raise ValueError(f"Evaluation detailed file missing Status column: {detailed}")
    counts = table["Status"].astype(str).value_counts().to_dict()
    unknown = sorted(set(counts) - {"TP", "FP", "FN"})
    if unknown:
        raise ValueError(f"Evaluation detailed file has unsupported Status values {unknown}: {detailed}")
    tp = int(counts.get("TP", 0))
    fp = int(counts.get("FP", 0))
    fn = int(counts.get("FN", 0))
    return {
        "detailed_file": str(detailed),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "NO. of Truth Alleles": tp + fn,
        "NO. of Predicted Alleles": tp + fp,
    }


def collect_eval_summaries(
    results_root: Path,
    exprs: list[str],
    methods: list[str],
    genes: list[str],
    populations: list[str],
) -> pd.DataFrame:
    rows: list[dict] = []
    method_set = set(methods) if methods else None
    gene_set = {g.upper() for g in genes} if genes else None
    pop_set = set(populations) if populations else None
    for expr in exprs:
        base = results_root / "eval" / expr
        if not base.exists():
            continue
        for summary_file in sorted(base.glob("*/*/[VJ]/eval_*_summary.csv")):
            rel = summary_file.relative_to(base)
            method = rel.parts[0]
            pop = rel.parts[1]
            gene = rel.parts[2]
            if method_set is not None and method not in method_set:
                continue
            if gene_set is not None and gene not in gene_set:
                continue
            if pop_set is not None and pop not in pop_set:
                continue
            m = SUMMARY_RE.match(summary_file.name)
            seed = m.group("seed") if m else ""
            row = {
                "expr": expr,
                "method": method,
                "population": pop,
                "gene": gene,
                "seed": seed,
                "summary_file": str(summary_file),
            }
            row.update(read_summary_metrics(summary_file))
            row.update(read_detailed_counts(summary_file))
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_population_means(per_sample: pd.DataFrame) -> pd.DataFrame:
    if per_sample.empty:
        return pd.DataFrame()
    melted = per_sample.melt(
        id_vars=["expr", "method", "population", "gene", "seed"],
        value_vars=METRICS,
        var_name="metric",
        value_name="value",
    )
    rows = []
    for keys, sub in melted.groupby(["expr", "method", "gene", "metric"], sort=True):
        expr, method, gene, metric = keys
        pop_means = sub.groupby("population")["value"].mean()
        mean = pop_means.mean()
        std = pop_means.std(ddof=1)
        rows.append(
            {
                "expr": expr,
                "method": method,
                "gene": gene,
                "metric": metric,
                "mean_of_population_means": mean,
                "std_of_population_means": std,
                "n_populations": pop_means.size,
                "n_samples": sub.shape[0],
                "formatted": f"{mean:.3f}±{std:.3f}" if pd.notna(std) else f"{mean:.3f}±nan",
            }
        )
    return pd.DataFrame(rows)


def summarize_counts(per_sample: pd.DataFrame) -> pd.DataFrame:
    if per_sample.empty or "TP" not in per_sample.columns:
        return pd.DataFrame()
    rows = []
    count_cols = ["NO. of Truth Alleles", "NO. of Predicted Alleles", "TP", "FP", "FN"]
    for keys, sub in per_sample.groupby(["expr", "method", "gene"], sort=True):
        expr, method, gene = keys
        row = {"expr": expr, "method": method, "gene": gene, "n_samples": int(len(sub))}
        for col in count_cols:
            row[col] = int(pd.to_numeric(sub[col], errors="coerce").fillna(0).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def write_markdown(summary: pd.DataFrame, out_md: Path) -> None:
    out_md.parent.mkdir(parents=True, exist_ok=True)
    with out_md.open("w", encoding="utf-8") as fh:
        fh.write("# Evaluation Summary\n\n")
        if summary.empty:
            fh.write("No evaluation summary files were found.\n")
            return
        for expr, expr_df in summary.groupby("expr", sort=True):
            fh.write(f"## {expr}\n\n")
            pivot = expr_df.pivot_table(
                index=["method", "metric"],
                columns="gene",
                values="formatted",
                aggfunc="first",
            ).reset_index()
            fh.write("```csv\n")
            fh.write(pivot.to_csv(index=False, lineterminator="\n"))
            fh.write("```\n\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", required=True, help="Directory containing results/eval")
    ap.add_argument("--exprs", required=True)
    ap.add_argument("--methods", default="")
    ap.add_argument("--genes", default="V J")
    ap.add_argument("--populations", default="")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    results_root = Path(args.results_root)
    out_dir = Path(args.out_dir)
    exprs = split_words(args.exprs, [])
    methods = split_words(args.methods, [])
    genes = split_words(args.genes, ["V", "J"])
    populations = split_words(args.populations, [])

    per_sample = collect_eval_summaries(results_root, exprs, methods, genes, populations)
    out_dir.mkdir(parents=True, exist_ok=True)
    per_sample.to_csv(out_dir / "per_sample_eval_metrics.csv", index=False)
    summary = summarize_population_means(per_sample)
    summary.to_csv(out_dir / "population_mean_eval_summary.csv", index=False)
    counts = summarize_counts(per_sample)
    counts.to_csv(out_dir / "aggregate_eval_counts.csv", index=False)
    write_markdown(summary, out_dir / "population_mean_eval_summary.md")
    print(f"Wrote evaluation summaries to {out_dir}")


if __name__ == "__main__":
    main()
