#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build benchmark gene-aware recovery matrices.

Matrix 1 reports exact-match P/R/F1/TP/FP/FN by benchmark, method, and
TRBV gene category. Matrix 2 stratifies truth alleles into exact recovery,
evidence-compatible partial recovery, and not recovered/evidence-conflicting.

This script is a gene-level diagnostic. The final overall S2/S5 primary
sequence-recovery tables are rebuilt by rebuild_sequence_only_primary_tables.py,
where overall exact recovery ignores upstream gene-label assignment.

The second matrix keeps partial recovery separate from the standard P/R/F1
calculation, so observed-region compatibility does not change the exact-match
benchmark metrics.
"""

from __future__ import annotations

import math
import os
import re
import sys
import csv
import io
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


CODE_EXPERIMENT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = CODE_EXPERIMENT_DIR.parents[2]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / CODE_EXPERIMENT_DIR.name
REF_ROOT = PACKAGE_ROOT / "data" / "ref"
DEFAULT_SOURCE_ROOT = EXPERIMENT_DIR / "generated" / "full_provenance_inputs"
SOURCE_ROOT = Path(
    os.environ.get("PANTCR_FULL_PROVENANCE_ROOT", DEFAULT_SOURCE_ROOT)
).resolve()
OUT = EXPERIMENT_DIR / "generated" / "gene_level_recovery_matrices"
sys.path.append(str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))

from common_pantcr_io import (  # noqa: E402
    defining_positions,
    evidence_observed_bases_by_gene,
    load_default_refs,
    load_trb_index,
    trim_sequence,
)


METHOD_MAP = {
    "Bayes": "PanTCR",
    "BayesNoPrior": "PanTCR-NP",
    "PanTCRLeaveout": "PanTCR",
    "PanTCR_fixed_graph": "PanTCR",
    "PanTCR_NP_no_prior": "PanTCR-NP",
    "MiXCR": "MiXCR-default",
    "MiXCR-all": "MiXCR-all",
    "FindAlleles": "FindAlleles",
    "PanTCR.semi": "PanTCR",
}

INSILICO_DIRS = {
    "expr_ScenarioA": "scenario_a",
    "expr_ScenarioB": "scenario_b",
    "expr_ScenarioC": "scenario_c",
    "expr_FullLength": "full_length",
}

SCENARIO_LABEL_TO_EXPR = {
    "Scenario A": "expr_ScenarioA",
    "Scenario B": "expr_ScenarioB",
    "Scenario C": "expr_ScenarioC",
}

STRICT_PARTIAL_DEFINITION = (
    "The complete inferred sequence is not an exact match, but the sample has "
    "physical mutation evidence for at least one target-defining site; every "
    "covered target-defining site supports the target allele without conflict; "
    "and a same-gene prediction matches the target allele at all those covered "
    "target-defining sites. Uncovered positions are not treated as read-supported."
)


@dataclass(frozen=True)
class DatasetSpec:
    benchmark: str
    dataset: str
    truth_path: Path
    pred_path: Path
    sample_cols: tuple[str, ...]
    min_naive: float


def clean_method(value: object) -> str:
    return METHOD_MAP.get(str(value), str(value))


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, low_memory=False)


def safe_div(n: float, d: float) -> float:
    return float(n) / float(d) if d else math.nan


def safe_precision(tp: float, fp: float) -> float:
    denom = tp + fp
    if denom == 0:
        return 0.0
    return float(tp) / float(denom)


def safe_recall(tp: float, fn: float) -> float:
    denom = tp + fn
    if denom == 0:
        return math.nan
    return float(tp) / float(denom)


def f1_score(p: float, r: float) -> float:
    if math.isnan(p) or math.isnan(r):
        return math.nan
    if p + r == 0:
        return 0.0
    return 2 * p * r / (p + r)


def metric_text(value: float) -> str:
    return "not available" if math.isnan(value) else f"{value:.3f}"


def count_text(value: float | int) -> str:
    if pd.isna(value):
        return "not available"
    return str(int(value))


def configure_source_root(path: str | Path) -> None:
    global SOURCE_ROOT
    SOURCE_ROOT = Path(path).resolve()


def source_root_message() -> str:
    return (
        f"Full provenance source root is required but was not found: {SOURCE_ROOT}. "
        "Provide --source-root or set PANTCR_FULL_PROVENANCE_ROOT. The source root must contain "
        "01_count_matrix_and_coverage_strata/, 12_semi_simu_evidence_analysis/, and "
        "13_scbulk_real_evidence_analysis/."
    )


def relative_source_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(SOURCE_ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def sanitize_path_value(value: object) -> object:
    if value is None or pd.isna(value):
        return value
    text = str(value)
    if not text:
        return text
    text_path = Path(text)
    if text_path.is_absolute():
        for root in [SOURCE_ROOT, PACKAGE_ROOT]:
            try:
                return str(text_path.resolve().relative_to(root)).replace("\\", "/")
            except ValueError:
                pass
        return text
    normalized_parts = Path(text.replace("\\", "/")).parts
    source_anchor = SOURCE_ROOT.name
    if source_anchor and source_anchor in normalized_parts:
        anchor_idx = normalized_parts.index(source_anchor)
        return "/".join(normalized_parts[anchor_idx + 1:])
    return text.replace("\\", "/")


def sanitize_source_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    path_like = [
        col for col in out.columns
        if col.endswith("_file")
        or col.endswith("_path")
        or col.lower() in {"source truth file", "source prediction file"}
    ]
    for col in path_like:
        out[col] = out[col].map(sanitize_path_value)
    return out


def sample_key(row: pd.Series, cols: Iterable[str]) -> tuple[str, ...]:
    return tuple(str(row.get(col, "")) for col in cols)


def build_dataset_specs() -> list[DatasetSpec]:
    specs: list[DatasetSpec] = []
    base = SOURCE_ROOT / "01_count_matrix_and_coverage_strata"
    for expr, subdir in INSILICO_DIRS.items():
        d = base / subdir
        specs.append(
            DatasetSpec(
                benchmark="In silico simulation",
                dataset=expr,
                truth_path=d / "per_truth_call_status.csv",
                pred_path=d / "per_prediction_status.csv",
                sample_cols=("expr", "method_clean", "population", "seed", "gene_type"),
                min_naive=1.0,
            )
        )

    specs.append(
        DatasetSpec(
            benchmark="Semi-synthetic AIRR-seq",
            dataset="AIRR pooled",
            truth_path=SOURCE_ROOT / "12_semi_simu_evidence_analysis" / "results" / "per_truth_call_status.csv",
            pred_path=SOURCE_ROOT / "12_semi_simu_evidence_analysis" / "results" / "per_prediction_status.csv",
            sample_cols=("expr", "dataset", "method_clean", "sample_id", "gene_type"),
            min_naive=2.0,
        )
    )
    specs.append(
        DatasetSpec(
            benchmark="Pseudo-bulk RNA-seq",
            dataset="SC pooled",
            truth_path=SOURCE_ROOT / "13_scbulk_real_evidence_analysis" / "results" / "per_truth_call_status.csv",
            pred_path=SOURCE_ROOT / "13_scbulk_real_evidence_analysis" / "results" / "per_prediction_status.csv",
            sample_cols=("SC_ID", "DatasetID", "method_clean", "gene_type"),
            min_naive=2.0,
        )
    )
    return specs


def normalize_tables(spec: DatasetSpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    truth = read_csv(spec.truth_path)
    pred = read_csv(spec.pred_path)

    truth = truth[truth["gene_type"].eq("V")].copy()
    pred = pred[pred["gene_type"].eq("V")].copy()

    truth["Benchmark"] = spec.benchmark
    pred["Benchmark"] = spec.benchmark
    if spec.benchmark == "Semi-synthetic AIRR-seq":
        truth["Dataset"] = truth["dataset"].astype(str)
        pred["Dataset"] = pred["dataset"].astype(str)
    elif spec.benchmark == "Pseudo-bulk RNA-seq":
        truth["Dataset"] = truth["SC_ID"].astype(str)
        pred["Dataset"] = pred["SC_ID"].astype(str)
    else:
        truth["Dataset"] = spec.dataset
        pred["Dataset"] = spec.dataset

    truth["method_clean"] = truth["method"].map(clean_method)
    pred["method_clean"] = pred["method"].map(clean_method)
    truth["Method"] = truth["method_clean"]
    pred["Method"] = pred["method_clean"]
    truth["Gene category"] = truth["gene"].astype(str)
    pred["Gene category"] = pred["gene"].astype(str)
    truth["Source truth file"] = relative_source_path(spec.truth_path)
    pred["Source prediction file"] = relative_source_path(spec.pred_path)
    truth["Min naive evidence filter"] = spec.min_naive
    return truth, pred


def add_all_gene_rows(df: pd.DataFrame, gene_col: str = "Gene category") -> pd.DataFrame:
    all_df = df.copy()
    all_df[gene_col] = "All TRBV"
    return pd.concat([df, all_df], ignore_index=True, sort=False)


def build_standard_matrix(truth_all: pd.DataFrame, pred_all: pd.DataFrame) -> pd.DataFrame:
    truth_for_group = add_all_gene_rows(truth_all)
    pred_for_group = add_all_gene_rows(pred_all)
    group_cols = ["Benchmark", "Dataset", "Method", "Gene category"]

    truth_counts = (
        truth_for_group.groupby(group_cols, dropna=False)
        .agg(
            Truth_count=("status", "size"),
            TP=("status", lambda s: int((s == "exact_tp").sum())),
        )
        .reset_index()
    )
    truth_counts["FN"] = truth_counts["Truth_count"] - truth_counts["TP"]

    pred_counts = (
        pred_for_group.groupby(group_cols, dropna=False)
        .agg(
            Prediction_count=("status", "size"),
            FP=("status", lambda s: int((s == "fp").sum())),
        )
        .reset_index()
    )

    out = truth_counts.merge(pred_counts, on=group_cols, how="outer").fillna(0)
    for col in ["Truth_count", "Prediction_count", "TP", "FP", "FN"]:
        out[col] = out[col].astype(int)
    out["Precision"] = [safe_precision(tp, fp) for tp, fp in zip(out["TP"], out["FP"])]
    out["Recall"] = [safe_recall(tp, fn) for tp, fn in zip(out["TP"], out["FN"])]
    out["F1_score"] = [f1_score(p, r) for p, r in zip(out["Precision"], out["Recall"])]
    out["Metric source"] = "Count-derived exact-match matrix; partial recovery is not included in TP."
    out["Precision"] = out["Precision"].map(metric_text)
    out["Recall"] = out["Recall"].map(metric_text)
    out["F1 score"] = out["F1_score"].map(metric_text)
    out = out.drop(columns=["F1_score"])
    return out.sort_values(group_cols).reset_index(drop=True)


def count_metric_row(row: pd.Series, gene_category: str, source: str) -> dict[str, object]:
    truth = int(row["NO. of Truth Alleles"])
    pred = int(row["NO. of Predicted Alleles"])
    tp = int(row["TP"])
    fp = int(row["FP"])
    fn = int(row["FN"])
    precision = safe_precision(tp, fp)
    recall = safe_recall(tp, fn)
    return {
        "Benchmark": "In silico simulation",
        "Dataset": SCENARIO_LABEL_TO_EXPR.get(str(row["Dataset"]), str(row["Dataset"])),
        "Method": "MiXCR-all",
        "Gene category": gene_category,
        "Truth_count": truth,
        "TP": tp,
        "FN": fn,
        "Prediction_count": pred,
        "FP": fp,
        "Precision": metric_text(precision),
        "Recall": metric_text(recall),
        "F1 score": metric_text(f1_score(precision, recall)),
        "Metric source": source,
    }


def mixcr_all_regenerated_rows() -> pd.DataFrame:
    return pd.DataFrame()


def keep_required_scope(df: pd.DataFrame, include_mixcr_all: bool) -> pd.DataFrame:
    """Keep only benchmark/method combinations requested for gene-level recovery."""
    required = {
        "In silico simulation": {"MiXCR-default", "FindAlleles", "PanTCR-NP", "PanTCR"},
        "Semi-synthetic AIRR-seq": {"MiXCR-default", "FindAlleles", "PanTCR"},
        "Pseudo-bulk RNA-seq": {"MiXCR-default", "FindAlleles", "PanTCR"},
    }
    if include_mixcr_all:
        required["In silico simulation"] = set(required["In silico simulation"]) | {"MiXCR-all"}
    mask = [
        str(row["Method"]) in required.get(str(row["Benchmark"]), set())
        for _, row in df.iterrows()
    ]
    return df.loc[mask].copy()


def load_population_mean_metrics(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8")
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    current_expr: str | None = None
    lines = text.splitlines()
    idx = 0
    while idx < len(lines):
        line = lines[idx]
        expr_match = re.match(r"^##\s+(expr_[A-Za-z0-9_]+)", line.strip())
        if expr_match:
            current_expr = expr_match.group(1)
            idx += 1
            continue
        if current_expr and line.strip() == "```csv":
            idx += 1
            csv_lines = []
            while idx < len(lines) and lines[idx].strip() != "```":
                csv_lines.append(lines[idx])
                idx += 1
            reader = csv.DictReader(io.StringIO("\n".join(csv_lines)))
            for row in reader:
                method = clean_method(row.get("method", ""))
                metric = str(row.get("metric", ""))
                value = str(row.get("V", ""))
                if method and metric in {"Precision", "Recall", "F1-score"} and value:
                    lookup.setdefault((current_expr, method), {})[metric] = value
        idx += 1
    return lookup


def load_full_length_metrics(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 6:
        return {}
    method_row = rows[2]
    metric_rows = rows[3:6]
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for col_idx in range(1, min(5, len(method_row))):
        method = clean_method(method_row[col_idx])
        for metric_row in metric_rows:
            if not metric_row:
                continue
            metric = metric_row[0]
            if metric in {"Precision", "Recall", "F1-score"} and col_idx < len(metric_row):
                lookup.setdefault(("expr_FullLength", method), {})[metric] = metric_row[col_idx]
    return lookup


def insilico_official_metrics() -> dict[tuple[str, str], dict[str, str]]:
    return {}


def semi_official_metrics() -> dict[tuple[str, str], dict[str, str]]:
    p = SOURCE_ROOT / "12_semi_simu_evidence_analysis" / "complete_metrics" / "method_level_complete_summary.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p)
    df = df[df["gene_type"].eq("V")].copy()
    df["Method"] = df["method"].map(clean_method)
    df = df[df["Method"].isin({"MiXCR-default", "FindAlleles", "PanTCR"})]
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for _, r in df.iterrows():
        lookup[(str(r["dataset"]), str(r["Method"]))] = {
            "Precision": f"{float(r['precision']):.3f}±{float(r.get('precision_sd', 0)):.3f}",
            "Recall": f"{float(r['recall']):.3f}±{float(r.get('recall_sd', 0)):.3f}",
            "F1-score": f"{float(r['f1']):.3f}±{float(r.get('f1_sd', 0)):.3f}",
        }
    return lookup


def pseudo_official_metrics() -> dict[tuple[str, str], dict[str, str]]:
    p = SOURCE_ROOT / "13_scbulk_real_evidence_analysis" / "complete_metrics" / "method_level_complete_summary.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p)
    df = df[df["gene_type"].eq("V")].copy()
    df["Method"] = df["method"].map(clean_method)
    df = df[df["Method"].isin({"MiXCR-default", "FindAlleles", "PanTCR"})]
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for _, r in df.iterrows():
        lookup[(str(r["SC_ID"]), str(r["Method"]))] = {
            "Precision": f"{float(r['precision']):.3f}",
            "Recall": f"{float(r['recall']):.3f}",
            "F1-score": f"{float(r['f1']):.3f}",
        }
    return lookup


def apply_official_total_metrics(standard: pd.DataFrame) -> pd.DataFrame:
    out = standard.copy()
    sim = insilico_official_metrics()
    semi = semi_official_metrics()
    pseudo = pseudo_official_metrics()
    for idx, row in out.iterrows():
        if str(row.get("Gene category")) != "All TRBV" or str(row.get("Method")) == "MiXCR-all":
            continue
        key = (str(row.get("Dataset")), str(row.get("Method")))
        vals = None
        if row.get("Benchmark") == "In silico simulation":
            vals = sim.get(key)
        elif row.get("Benchmark") == "Semi-synthetic AIRR-seq":
            vals = semi.get(key)
        elif row.get("Benchmark") == "Pseudo-bulk RNA-seq":
            vals = pseudo.get(key)
        if not vals:
            continue
        out.at[idx, "Precision"] = vals["Precision"]
        out.at[idx, "Recall"] = vals["Recall"]
        out.at[idx, "F1 score"] = vals["F1-score"]
        out.at[idx, "Metric source"] = (
            "All-TRBV P/R/F1 from regenerated source summaries; TP/FP/FN are exact-match total counts. "
            "Partial recovery is not included in TP."
        )
    return out


def build_prediction_index(pred_all: pd.DataFrame, sample_cols: tuple[str, ...]) -> dict[tuple[str, ...], list[str]]:
    out: dict[tuple[str, ...], list[str]] = {}
    for _, row in pred_all.iterrows():
        key = sample_key(row, sample_cols) + (str(row.get("gene", "")),)
        out.setdefault(key, []).append(str(row.get("pred_seq", "")).strip().upper())
    return out


def clean_supported_positions(row: pd.Series, observed_cache: dict, index_df, default_refs_raw) -> list[int]:
    gene = str(row.get("gene", ""))
    truth_seq = str(row.get("truth_seq", "")).strip().upper()
    gene_type = str(row.get("gene_type", ""))
    ref = trim_sequence(default_refs_raw.get(gene, ""), gene, gene_type, index_df)
    if not truth_seq or not ref:
        return []
    def_pos = defining_positions(truth_seq, ref)
    evidence_path = Path(str(row.get("evidence_file", "")))
    if not evidence_path.exists():
        return []
    min_naive = float(row.get("Min naive evidence filter", 2.0))
    cache_key = (str(evidence_path), min_naive)
    if cache_key not in observed_cache:
        observed_cache[cache_key] = evidence_observed_bases_by_gene(evidence_path, min_naive=min_naive)
    gene_observed = observed_cache[cache_key].get(gene, {})
    positions = []
    for pos in def_pos:
        bases = {str(base).upper() for base in gene_observed.get(pos, set()) if str(base).strip()}
        truth_base = truth_seq[pos] if 0 <= pos < len(truth_seq) else ""
        if bases and truth_base and bases == {truth_base}:
            positions.append(pos)
    return positions


def pred_matches_positions(pred_seq: str, truth_seq: str, positions: list[int]) -> bool:
    if not positions:
        return False
    pred_seq = str(pred_seq).strip().upper()
    truth_seq = str(truth_seq).strip().upper()
    return all(pos < len(pred_seq) and pos < len(truth_seq) and pred_seq[pos] == truth_seq[pos] for pos in positions)


def classify_truth_recovery_for_spec(spec: DatasetSpec, truth: pd.DataFrame, pred: pd.DataFrame) -> pd.DataFrame:
    pred_index = build_prediction_index(pred, spec.sample_cols)
    index_df = load_trb_index(REF_ROOT / "TRB_index.csv")
    default_refs_raw = load_default_refs(REF_ROOT / "IMGT_TRB_default.csv")
    observed_cache = {}

    categories = []
    partial_candidate_counts = []
    for _, row in truth.iterrows():
        if str(row.get("status", "")) == "exact_tp":
            categories.append("Exact recovery")
            partial_candidate_counts.append(0)
            continue

        same_gene_preds = pred_index.get(sample_key(row, spec.sample_cols) + (str(row.get("gene", "")),), [])
        if not same_gene_preds:
            categories.append("Not recovered or evidence-conflicting")
            partial_candidate_counts.append(0)
            continue

        covered = int(row.get("covered_defining_variant_count", 0) or 0)
        matched = int(row.get("matched_covered_defining_variant_count", 0) or 0)
        conflicting = int(row.get("conflicting_covered_defining_variant_count", 0) or 0)
        mixed = int(row.get("mixed_covered_defining_variant_count", 0) or 0)
        evidence_compatible = covered > 0 and matched == covered and conflicting == 0 and mixed == 0
        if not evidence_compatible:
            categories.append("Not recovered or evidence-conflicting")
            partial_candidate_counts.append(0)
            continue

        positions = clean_supported_positions(row, observed_cache, index_df, default_refs_raw)
        truth_seq = str(row.get("truth_seq", "")).strip().upper()
        candidates = [seq for seq in same_gene_preds if pred_matches_positions(seq, truth_seq, positions)]
        if candidates:
            categories.append("Evidence-compatible partial recovery")
            partial_candidate_counts.append(len(candidates))
        else:
            categories.append("Not recovered or evidence-conflicting")
            partial_candidate_counts.append(0)

    out = truth.copy()
    out["Truth recovery class"] = categories
    out["Strict partial candidate count"] = partial_candidate_counts
    out["Strict partial definition"] = STRICT_PARTIAL_DEFINITION
    return out


def build_truth_recovery_matrix(classified_truth: pd.DataFrame) -> pd.DataFrame:
    truth_for_group = add_all_gene_rows(classified_truth)
    group_cols = ["Benchmark", "Dataset", "Method", "Gene category"]
    matrix = (
        truth_for_group.groupby(group_cols + ["Truth recovery class"], dropna=False)
        .size()
        .reset_index(name="n")
        .pivot_table(index=group_cols, columns="Truth recovery class", values="n", fill_value=0, aggfunc="sum")
        .reset_index()
    )
    matrix.columns.name = None
    for col in [
        "Exact recovery",
        "Evidence-compatible partial recovery",
        "Not recovered or evidence-conflicting",
    ]:
        if col not in matrix.columns:
            matrix[col] = 0
    matrix["Truth count"] = matrix[
        [
            "Exact recovery",
            "Evidence-compatible partial recovery",
            "Not recovered or evidence-conflicting",
        ]
    ].sum(axis=1)
    matrix["Exact recovery rate"] = [
        metric_text(safe_div(n, d))
        for n, d in zip(matrix["Exact recovery"], matrix["Truth count"])
    ]
    matrix["Exact or partial recovery rate"] = [
        metric_text(safe_div(e + p, d))
        for e, p, d in zip(
            matrix["Exact recovery"],
            matrix["Evidence-compatible partial recovery"],
            matrix["Truth count"],
        )
    ]
    matrix["Definition note"] = STRICT_PARTIAL_DEFINITION
    return matrix.sort_values(group_cols).reset_index(drop=True)


def write_markdown(standard: pd.DataFrame, recovery: pd.DataFrame, out_path: Path) -> None:
    def markdown_table(df: pd.DataFrame) -> str:
        if df.empty:
            return "_No rows._"
        text_df = df.copy().where(pd.notna(df), "")
        cols = [str(c) for c in text_df.columns]
        rows = [[str(v) for v in row] for row in text_df.to_numpy().tolist()]
        lines = [
            "| " + " | ".join(cols) + " |",
            "| " + " | ".join(["---"] * len(cols)) + " |",
        ]
        for row in rows:
            lines.append("| " + " | ".join(cell.replace("|", "\\|") for cell in row) + " |")
        return "\n".join(lines)

    all_standard = standard[standard["Gene category"].eq("All TRBV")].copy()
    all_recovery = recovery[recovery["Gene category"].eq("All TRBV")].copy()
    lines = [
        "# gene-level recovery two-matrix implementation",
        "",
        "## Matrix 1: gene-aware exact-match performance",
        "",
        "This matrix reports P/R/F1/TP/FP/FN under a gene-aware exact-match diagnostic. Partial recovery is not counted as TP.",
        "",
        markdown_table(all_standard.head(30)),
        "",
        "## Matrix 2: truth recovery evidence stratification",
        "",
        "This matrix decomposes truth alleles into exact recovery, evidence-compatible partial recovery, and not recovered/evidence-conflicting. It does not redefine the P/R/F1 metrics.",
        "",
        markdown_table(all_recovery.head(30)),
        "",
        "## Strict partial definition",
        "",
        STRICT_PARTIAL_DEFINITION,
        "",
        "Complete per-TRBV gene rows are available in the accompanying Excel workbook and CSV files.",
        "",
    ]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        default=None,
        help=(
            "Full provenance root containing the normalized per-truth/per-prediction inputs. "
            "If omitted, PANTCR_FULL_PROVENANCE_ROOT is used; otherwise the package-local "
            "generated/full_provenance_inputs directory is expected."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.source_root:
        configure_source_root(args.source_root)
    if not SOURCE_ROOT.exists():
        raise FileNotFoundError(source_root_message())

    OUT.mkdir(parents=True, exist_ok=True)
    truth_tables = []
    pred_tables = []
    classified_tables = []

    for spec in build_dataset_specs():
        truth, pred = normalize_tables(spec)
        truth_tables.append(truth)
        pred_tables.append(pred)
        classified_tables.append(classify_truth_recovery_for_spec(spec, truth, pred))

    truth_all = pd.concat(truth_tables, ignore_index=True, sort=False)
    pred_all = pd.concat(pred_tables, ignore_index=True, sort=False)
    classified_truth = pd.concat(classified_tables, ignore_index=True, sort=False)

    standard = build_standard_matrix(truth_all, pred_all)
    mixcr_all = mixcr_all_regenerated_rows()
    if not mixcr_all.empty:
        standard = pd.concat([standard, mixcr_all], ignore_index=True, sort=False)
        standard = standard.sort_values(["Benchmark", "Dataset", "Method", "Gene category"]).reset_index(drop=True)
    standard = keep_required_scope(standard, include_mixcr_all=True)
    standard = apply_official_total_metrics(standard)

    recovery = build_truth_recovery_matrix(classified_truth)
    recovery = keep_required_scope(recovery, include_mixcr_all=False)

    standard_csv = OUT / "standard_performance_matrix_by_trbv_gene.csv"
    recovery_csv = OUT / "truth_recovery_matrix_by_trbv_gene.csv"
    classified_csv = OUT / "strict_partial_per_truth_classification.csv"
    xlsx = OUT / "gene_level_recovery_matrices.xlsx"
    md = OUT / "gene_level_recovery_matrices_summary.md"

    standard.to_csv(standard_csv, index=False)
    recovery.to_csv(recovery_csv, index=False)
    sanitize_source_columns(classified_truth).to_csv(classified_csv, index=False)
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        standard.to_excel(writer, sheet_name="Standard performance matrix", index=False)
        recovery.to_excel(writer, sheet_name="Truth recovery matrix", index=False)
        pd.DataFrame(
            [
                {
                    "Term": "Evidence-compatible partial recovery",
                    "Definition": STRICT_PARTIAL_DEFINITION,
                    "Scope note": "Not included in TP/P/R/F1; used only in the truth recovery evidence matrix.",
                },
                {
                    "Term": "MiXCR-all limitation",
                    "Definition": "MiXCR-all rows require regenerated MiXCR-all evaluation details.",
                    "Scope note": "Provide regenerated source inputs before adding MiXCR-all rows to this diagnostic workbook.",
                },
            ]
        ).to_excel(writer, sheet_name="Definitions and notes", index=False)
    write_markdown(standard, recovery, md)

    print(standard_csv)
    print(recovery_csv)
    print(classified_csv)
    print(xlsx)
    print(md)


if __name__ == "__main__":
    main()
