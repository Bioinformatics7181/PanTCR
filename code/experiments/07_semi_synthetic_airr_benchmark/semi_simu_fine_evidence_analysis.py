#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fine evidence-support analysis for manuscript semi-synthetic AIRR-1/AIRR-2 results."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[3]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / Path(__file__).resolve().parent.name
sys.path.append(str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))
from common_pantcr_io import (  # noqa: E402
    compatible_match,
    compatible_match_score,
    coverage_stratum,
    defining_positions,
    evidence_coverage_by_gene,
    evidence_observed_bases_by_gene,
    evidence_support_counts,
    evidence_support_stratum,
    load_default_refs,
    load_pmtr_sequences,
    load_trb_index,
    prediction_records,
    sequence_length_difference,
    trim_sequence,
    truth_records,
)


EXPR_LABELS = {"expr_AIRR1": "AIRR-1", "expr_AIRR2": "AIRR-2"}
METHODS = ["MiXCR", "FindAlleles", "Bayes", "BayesNoPrior"]
GENES = ["V", "J"]


def copy_refs(ref_root: Path, workspace_root: Path) -> Path:
    dst = workspace_root / "ref"
    dst.mkdir(parents=True, exist_ok=True)
    for name in ["TRB_index.csv", "IMGT_TRB_default.csv", "pmTR_TRB_V_J_cleaned.csv"]:
        src = ref_root / name
        if src.exists():
            shutil.copy2(src, dst / name)
    return dst


def sample_id_from_infer(infer_csv: Path, gene: str) -> str:
    name = infer_csv.name
    prefix = "infer_"
    suffix = f".{gene}.csv"
    if name.startswith(prefix) and name.endswith(suffix):
        return name[len(prefix):-len(suffix)]
    return name


def validate_infer_schema(infer_csv: Path) -> int:
    df = pd.read_csv(infer_csv)
    required = {"gene"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Inference CSV missing required columns {missing}: {infer_csv}")
    seq_cols = [c for c in ["seq_A", "seq_B"] if c in df.columns]
    if not seq_cols:
        raise ValueError(f"Inference CSV has no sequence columns seq_A/seq_B: {infer_csv}")
    return int(len(df))


def find_validation_fold(results_root: Path, expr: str, gene: str, evidence_name: str) -> Optional[int]:
    base = results_root / "validation" / expr / gene
    if not base.exists():
        return None
    for fold_dir in sorted(base.glob("fold_*")):
        if (fold_dir / evidence_name).exists():
            try:
                return int(fold_dir.name.split("_")[-1])
            except Exception:
                return None
    return None


def match_truth_and_predictions(truth, preds, gene_type):
    truth_rows_by_index = {}
    pred_rows = []
    used_preds = set()

    for ti, t in enumerate(truth):
        for pi, p in enumerate(preds):
            if pi in used_preds or p["gene"] != t["gene"]:
                continue
            if p["seq"] == t["seq"]:
                used_preds.add(pi)
                truth_rows_by_index[ti] = (ti, pi, "exact_tp", 0, 0, 0)
                break

    for ti, t in enumerate(truth):
        if ti in truth_rows_by_index:
            continue

        partial = None
        best_tie_count = 0
        compatible_candidates = []
        for pi, p in enumerate(preds):
            if pi in used_preds or p["gene"] != t["gene"]:
                continue
            score = compatible_match_score(t["seq"], p["seq"], gene_type)
            if score > 0:
                compatible_candidates.append((score, len(p["seq"]), pi))
        if compatible_candidates:
            compatible_candidates = sorted(compatible_candidates, reverse=True)
            best_score = compatible_candidates[0][0]
            best_tie_count = sum(1 for score, _length, _pi in compatible_candidates if score == best_score)
            partial = compatible_candidates[0][2]
        if partial is not None:
            used_preds.add(partial)
            truth_rows_by_index[ti] = (ti, partial, "partial_recovery", len(compatible_candidates), compatible_match_score(t["seq"], preds[partial]["seq"], gene_type), best_tie_count)
        else:
            same_gene_preds = [p for p in preds if p["gene"] == t["gene"]]
            truth_rows_by_index[ti] = (ti, None, "no_call" if not same_gene_preds else "fn", 0, 0, 0)

    for pi, _p in enumerate(preds):
        pred_rows.append((pi, "matched" if pi in used_preds else "fp"))
    truth_rows = [truth_rows_by_index[ti] for ti in range(len(truth))]
    return truth_rows, pred_rows


def metric_from_eval_details(results_root: Path, expr: str) -> pd.DataFrame:
    rows = []
    required_metric_keys = {"Precision", "Recall", "F1-score"}
    for method in METHODS:
        for gene in GENES:
            eval_dir = results_root / "eval" / expr / method / "MIX" / gene
            if not eval_dir.exists():
                raise FileNotFoundError(f"Missing evaluation directory for {expr}/{method}/{gene}: {eval_dir}")
            per_file_metrics = []
            status_counts: Dict[str, int] = {}
            n_files = 0
            n_summary_files = 0
            for summary in sorted(eval_dir.glob("*_summary.csv")):
                df = pd.read_csv(summary)
                if not {"Metric", "Value"}.issubset(df.columns):
                    raise ValueError(f"Evaluation summary missing Metric/Value columns: {summary}")
                vals = dict(zip(df["Metric"].astype(str), df["Value"]))
                missing = sorted(required_metric_keys - set(vals))
                if missing:
                    raise ValueError(f"Evaluation summary missing required metrics {missing}: {summary}")
                n_summary_files += 1
                per_file_metrics.append((
                    float(vals["Precision"]),
                    float(vals["Recall"]),
                    float(vals["F1-score"]),
                ))
            for detailed in sorted(eval_dir.glob("*_detailed.csv")):
                df = pd.read_csv(detailed)
                if "Status" not in df.columns:
                    raise ValueError(f"Evaluation detailed file missing Status column: {detailed}")
                unknown_status = sorted(set(df["Status"].dropna().astype(str)) - {"TP", "FP", "FN"})
                if unknown_status:
                    raise ValueError(f"Evaluation detailed file has unsupported Status values {unknown_status}: {detailed}")
                n_files += 1
                for status, n in df["Status"].astype(str).value_counts().items():
                    status_counts[status] = status_counts.get(status, 0) + int(n)
            if n_files == 0 or n_summary_files == 0:
                raise RuntimeError(f"No evaluation files found for {expr}/{method}/{gene}: {eval_dir}")
            if n_files != n_summary_files:
                raise RuntimeError(f"Detailed/summary evaluation file count mismatch for {expr}/{method}/{gene}: detailed={n_files}, summary={n_summary_files}")
            tp = status_counts.get("TP", 0)
            fp = status_counts.get("FP", 0)
            fn = status_counts.get("FN", 0)
            precision = tp / (tp + fp) if (tp + fp) else 0.0
            recall = tp / (tp + fn) if (tp + fn) else 0.0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
            if per_file_metrics:
                metric_df = pd.DataFrame(per_file_metrics, columns=["precision", "recall", "f1"])
                mean_precision = float(metric_df["precision"].mean())
                mean_recall = float(metric_df["recall"].mean())
                mean_f1 = float(metric_df["f1"].mean())
                sd_precision = float(metric_df["precision"].std(ddof=0))
                sd_recall = float(metric_df["recall"].std(ddof=0))
                sd_f1 = float(metric_df["f1"].std(ddof=0))
            else:
                mean_precision = precision
                mean_recall = recall
                mean_f1 = f1
                sd_precision = sd_recall = sd_f1 = 0.0
            rows.append({
                "expr": expr,
                "dataset": EXPR_LABELS.get(expr, expr),
                "method": method,
                "gene_type": gene,
                "n_eval_files": n_files,
                "n_summary_files": n_summary_files,
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "precision": mean_precision,
                "recall": mean_recall,
                "f1": mean_f1,
                "precision_sd": sd_precision,
                "recall_sd": sd_recall,
                "f1_sd": sd_f1,
                "aggregate_precision_from_counts": precision,
                "aggregate_recall_from_counts": recall,
                "aggregate_f1_from_counts": f1,
            })
    return pd.DataFrame(rows)


def analyze_expr(results_root: Path, ref_dir: Path, expr: str, min_naive: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    index_df = load_trb_index(ref_dir / "TRB_index.csv")
    pmtr_map = load_pmtr_sequences(ref_dir / "pmTR_TRB_V_J_cleaned.csv")
    default_refs = load_default_refs(ref_dir / "IMGT_TRB_default.csv")

    truth_rows = []
    pred_rows = []

    for method in METHODS:
        for gene_type in GENES:
            infer_dir = results_root / "infer" / expr / method / "MIX" / gene_type
            if not infer_dir.exists():
                raise FileNotFoundError(f"Missing inference directory for {expr}/{method}/{gene_type}: {infer_dir}")
            for infer_csv in sorted(infer_dir.glob("*.csv")):
                infer_table_rows = validate_infer_schema(infer_csv)
                sample_id = sample_id_from_infer(infer_csv, gene_type)
                label_csv = results_root / "labels" / expr / "MIX" / f"genotype_{sample_id}.csv"
                evidence_csv = results_root / "mutations" / expr / "MIX" / gene_type / f"{sample_id}.{gene_type}_sequences.csv"
                if not label_csv.exists():
                    raise FileNotFoundError(f"Missing genotype label for {expr}/{gene_type}/{sample_id}: {label_csv}")
                if not evidence_csv.exists():
                    raise FileNotFoundError(f"Missing mutation evidence for {expr}/{gene_type}/{sample_id}: {evidence_csv}")
                coverage = evidence_coverage_by_gene(evidence_csv, min_naive=min_naive)
                observed = evidence_observed_bases_by_gene(evidence_csv, min_naive=min_naive)
                fold = find_validation_fold(results_root, expr, gene_type, evidence_csv.name) if evidence_csv.exists() else None

                truth = truth_records(label_csv, gene_type, index_df, pmtr_map, prefer_label_sequence=True)
                preds = prediction_records(infer_csv, gene_type, index_df)
                matched_truth, matched_preds = match_truth_and_predictions(truth, preds, gene_type)

                for ti, pi, status, compatible_candidate_count, compatible_match_len, compatible_best_tie_count in matched_truth:
                    t = truth[ti]
                    gene = t["gene"]
                    ref = trim_sequence(default_refs.get(gene, ""), gene, gene_type, index_df)
                    if not ref:
                        raise ValueError(f"Missing default reference sequence for {expr}/{gene_type}/{gene}")
                    def_pos = defining_positions(t["seq"], ref)
                    length_diff = sequence_length_difference(t["seq"], ref)
                    covered = coverage.get(gene, set())
                    support_counts = evidence_support_counts(t["seq"], def_pos, observed.get(gene, {}))
                    truth_rows.append({
                        "expr": expr,
                        "dataset": EXPR_LABELS.get(expr, expr),
                        "truth_unit": "unique allele sequence per gene per sample in evaluation scope",
                        "partial_status_definition": "prefix/suffix sequence-compatible partial, not conservative defining-site recovery",
                        "partial_matching_rule": "best compatible prediction by prefix/suffix overlap length after exact matches",
                        "method": method,
                        "sample_id": sample_id,
                        "population": "MIX",
                        "fold": fold,
                        "gene_type": gene_type,
                        "gene": gene,
                        "truth_allele": t["allele"],
                        "truth_seq": t["seq"],
                        "truth_seq_full_catalog": t.get("seq_full_catalog", t["seq"]),
                        "truth_sequence_scope": t.get("truth_sequence_scope", "eval_trimmed"),
                        "matched_pred_seq": preds[pi]["seq"] if pi is not None else "",
                        "compatible_partial_candidate_count": compatible_candidate_count,
                        "compatible_partial_match_length": compatible_match_len,
                        "compatible_partial_best_tie_count": compatible_best_tie_count,
                        "status": status,
                        "coverage_stratum": coverage_stratum(def_pos, covered, length_difference=length_diff),
                        "evidence_support_stratum": evidence_support_stratum(
                            def_pos,
                            covered,
                            support_counts,
                            length_difference=length_diff,
                        ),
                        "n_defining_positions": len(def_pos),
                        "n_defining_positions_covered": sum(1 for p in def_pos if p in covered),
                        "eval_sequence_length_difference_vs_default_ref": length_diff,
                        "has_unresolved_eval_length_difference": bool(length_diff != 0),
                        **support_counts,
                        "single_site_compatible": (
                            support_counts["covered_defining_variant_count"] == 1
                            and support_counts["matched_covered_defining_variant_count"] == 1
                            and support_counts["conflicting_covered_defining_variant_count"] == 0
                        ),
                        "infer_file": str(infer_csv),
                        "infer_table_rows": infer_table_rows,
                        "label_file": str(label_csv),
                        "evidence_file": str(evidence_csv),
                    })

                for pi, status in matched_preds:
                    p = preds[pi]
                    pred_rows.append({
                        "expr": expr,
                        "dataset": EXPR_LABELS.get(expr, expr),
                        "prediction_evidence_annotation_scope": "prediction status only; false-positive prediction site-level evidence is not annotated in this experiment",
                        "method": method,
                        "sample_id": sample_id,
                        "gene_type": gene_type,
                        "gene": p["gene"],
                        "pred_allele": p["allele"],
                        "pred_seq": p["seq"],
                        "status": status,
                        "infer_file": str(infer_csv),
                    })

    return pd.DataFrame(truth_rows), pd.DataFrame(pred_rows)


def summarize_truth(truth_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    out = {}
    if truth_df.empty:
        return out
    truth_df = truth_df.copy()
    truth_df["recovered_exact_or_partial"] = truth_df["status"].isin(["exact_tp", "partial_recovery"])
    keys = ["expr", "dataset", "method", "gene_type", "evidence_support_stratum"]
    out["summary_by_evidence_support"] = (
        truth_df.groupby(keys)
        .agg(
            n_truth_alleles=("status", "size"),
            exact_recovery_n=("status", lambda s: int((s == "exact_tp").sum())),
            partial_recovery_n=("status", lambda s: int((s == "partial_recovery").sum())),
            recovered_exact_or_partial_n=("recovered_exact_or_partial", "sum"),
            recovered_exact_or_partial_rate=("recovered_exact_or_partial", "mean"),
            no_call_n=("status", lambda s: int((s == "no_call").sum())),
            fn_n=("status", lambda s: int((s == "fn").sum())),
            mixed_covered_defining_variant_count=("mixed_covered_defining_variant_count", "sum"),
            conflicting_covered_defining_variant_count=("conflicting_covered_defining_variant_count", "sum"),
            clean_matched_covered_defining_variant_count=("clean_matched_covered_defining_variant_count", "sum"),
        )
        .reset_index()
    )
    keys2 = ["expr", "dataset", "method", "gene_type", "coverage_stratum"]
    out["summary_by_coverage"] = (
        truth_df.groupby(keys2)
        .agg(
            n_truth_alleles=("status", "size"),
            recovered_exact_or_partial_n=("recovered_exact_or_partial", "sum"),
            recovered_exact_or_partial_rate=("recovered_exact_or_partial", "mean"),
        )
        .reset_index()
    )
    partial = truth_df[truth_df["status"].eq("partial_recovery")].copy()
    if partial.empty:
        out["partial_match_diagnostics"] = pd.DataFrame(
            columns=[
                "dataset",
                "method",
                "gene_type",
                "partial_recovery_n",
                "max_compatible_partial_candidate_count",
                "mean_compatible_partial_candidate_count",
                "max_compatible_partial_match_length",
                "max_compatible_partial_best_tie_count",
                "partial_rows_with_ties",
            ]
        )
    else:
        out["partial_match_diagnostics"] = (
            partial.groupby(["dataset", "method", "gene_type"])
            .agg(
                partial_recovery_n=("status", "size"),
                max_compatible_partial_candidate_count=("compatible_partial_candidate_count", "max"),
                mean_compatible_partial_candidate_count=("compatible_partial_candidate_count", "mean"),
                max_compatible_partial_match_length=("compatible_partial_match_length", "max"),
                max_compatible_partial_best_tie_count=("compatible_partial_best_tie_count", "max"),
                partial_rows_with_ties=("compatible_partial_best_tie_count", lambda s: int((pd.to_numeric(s, errors="coerce").fillna(0) > 1).sum())),
            )
            .reset_index()
        )
    return out


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |"]
    lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
    for _, row in df.iterrows():
        vals = ["" if pd.isna(row[c]) else str(row[c]) for c in cols]
        vals = [v.replace("\n", " ") for v in vals]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_report(out_dir: Path, metrics: pd.DataFrame, truth_df: pd.DataFrame, summaries: Dict[str, pd.DataFrame], min_naive: float) -> None:
    lines = []
    lines.append("# Semi-synthetic AIRR evidence-support analysis\n")
    lines.append("## Scope\n")
    lines.append(f"This report covers the corresponding semi-synthetic AIRR-1 (`expr_AIRR1`) and AIRR-2 (`expr_AIRR2`) result directories. The analysis reads labels, inference outputs, mutation evidence, and evaluation outputs generated by the upstream semi-synthetic pipeline. Mutation-evidence rows were filtered with `NaiveDiversityIndex >= {min_naive:g}` to match the real/semi-synthetic inference setting.\n")
    lines.append("## Local metric extraction from existing evaluation details\n")
    lines.append("Values below are per-sample means with descriptive per-sample standard deviations from generated evaluation CSVs. The SD values expose sample-level heterogeneity across rerun samples.\n")
    metric_show = metrics[metrics["method"].isin(["Bayes", "BayesNoPrior", "MiXCR", "FindAlleles"])].copy()
    metric_show["precision"] = metric_show.apply(lambda r: f"{r['precision']:.3f}±{r['precision_sd']:.3f}", axis=1)
    metric_show["recall"] = metric_show.apply(lambda r: f"{r['recall']:.3f}±{r['recall_sd']:.3f}", axis=1)
    metric_show["f1"] = metric_show.apply(lambda r: f"{r['f1']:.3f}±{r['f1_sd']:.3f}", axis=1)
    for col in ["aggregate_precision_from_counts", "aggregate_recall_from_counts", "aggregate_f1_from_counts"]:
        metric_show[col] = metric_show[col].map(lambda x: f"{x:.3f}")
    lines.append("The `n_eval_files` column counts detailed evaluation files used for aggregate TP/FP/FN counts, while `n_summary_files` counts per-sample summary files used for the precision/recall/F1 means. The `TP/FP/FN` columns are aggregate detailed-status counts. The `precision/recall/f1` columns are per-sample means, while `aggregate_*_from_counts` gives the corresponding count-pooled metric to avoid mixing denominators.\n")
    lines.append(markdown_table(metric_show[[
        "dataset", "method", "gene_type", "n_eval_files", "n_summary_files",
        "precision", "recall", "f1",
        "aggregate_precision_from_counts", "aggregate_recall_from_counts", "aggregate_f1_from_counts",
        "TP", "FP", "FN",
    ]]))
    lines.append("\n## Evidence-supported recovery strata\n")
    ev = summaries.get("summary_by_evidence_support", pd.DataFrame())
    bayes_v = ev[(ev["method"] == "Bayes") & (ev["gene_type"] == "V")].copy() if not ev.empty else ev
    if not bayes_v.empty:
        bayes_v["exact_or_sequence_compatible_partial_rate"] = bayes_v["recovered_exact_or_partial_rate"].map(lambda x: f"{x:.3f}")
        lines.append("Table scope: PanTCR Bayes V only; full method/gene strata are retained in `summary_by_evidence_support.csv`.\n")
        lines.append(markdown_table(bayes_v[["dataset", "evidence_support_stratum", "n_truth_alleles", "exact_recovery_n", "partial_recovery_n", "recovered_exact_or_partial_n", "exact_or_sequence_compatible_partial_rate", "mixed_covered_defining_variant_count", "conflicting_covered_defining_variant_count", "no_call_n", "fn_n"]]))
        lines.append("Conflict columns are summed defining-site counts across truth alleles, not allele counts.\n")
    lines.append("\n## Sequence-compatible partial-match diagnostics\n")
    partial_diag = summaries.get("partial_match_diagnostics", pd.DataFrame()).copy()
    if not partial_diag.empty:
        partial_diag["mean_compatible_partial_candidate_count"] = partial_diag["mean_compatible_partial_candidate_count"].map(lambda x: f"{float(x):.3f}")
    lines.append(markdown_table(partial_diag))
    lines.append("These diagnostics audit the sequence-compatible partial category. The current outputs expose candidate counts, match lengths, and tie counts so partial recovery is not confused with exact full-sequence recovery.\n")
    lines.append("\n## Interpretation\n")
    lines.append("- In semi-synthetic AIRR-1/AIRR-2, PanTCR performance is evaluated against known sample-level genotypes, so the evidence-support strata separate allele-defining sites observed in filtered PanTCR mutation-detection outputs from graph/prior-assisted truth-allele recovery contexts.\n")
    lines.append("- Truth alleles are counted as unique allele sequences per gene per sample, not allele-copy genotype counts.\n")
    lines.append("- The evidence-strata table above is shown for PanTCR Bayes V because the fragmentation concern is V-region allele recovery under fragmented evidence. Full method/gene strata are available in the CSV outputs; this experiment does not annotate false-positive prediction sites.\n")
    lines.append("- `fully_evidence_supported` means all truth-defining variant sites are covered and compatible after the same evidence filter used by inference. `partially_evidence_supported` means the covered defining sites are compatible but not all defining sites are covered. `prior_only_imputation_dependent` means defining sites are not covered in the mutation evidence and should not be described as directly observed.\n")
    lines.append("- If both the target base and an alternative base are observed at the same defining site, the site is counted as mixed/conflicting rather than cleanly supported, because such evidence should not be used as unqualified support for a candidate allele. The `conflicting_covered_evidence` stratum is therefore a cautionary category, not clean read support.\n")
    lines.append("- These tables complement the main Precision/Recall/F1 results by identifying contexts where fragmented coverage coincides with limited recovery; this is a stratified association, not a causal proof that coverage alone explains every miss.\n")
    lines.append("\n## Output files\n")
    for name in [
        "per_truth_call_status.csv",
        "per_prediction_status.csv",
        "metric_check_from_eval_details.csv",
        "summary_by_evidence_support.csv",
        "summary_by_coverage.csv",
        "partial_match_diagnostics.csv",
    ]:
        lines.append(f"- `{name}`")
    (out_dir / "semi_simu_fine_evidence_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    default_results_root = (
        PACKAGE_ROOT
        / "experiments"
        / "01_in_silico_trbv_benchmarks"
        / "benchmark_pipelines"
        / "02_semi_synthetic_airr"
        / "generated"
        / "results"
    )
    ap.add_argument(
        "--results-root",
        default=str(default_results_root),
        help="Expanded semi-synthetic benchmark results root generated by the benchmark pipeline.",
    )
    ap.add_argument("--ref-root", default=str(PACKAGE_ROOT / "data" / "ref"))
    ap.add_argument("--exprs", default="expr_AIRR1 expr_AIRR2")
    ap.add_argument("--workspace", default=str(EXPERIMENT_DIR / "generated" / "workspace"))
    ap.add_argument("--out-dir", default=str(EXPERIMENT_DIR / "generated" / "evidence_analysis"))
    ap.add_argument("--min-naive", type=float, default=2)
    args = ap.parse_args()

    exprs = args.exprs.split()
    workspace = Path(args.workspace)
    candidate_results_root = Path(args.results_root)
    if candidate_results_root.exists():
        results_root = candidate_results_root
    else:
        raise FileNotFoundError(
            "Missing semi-synthetic expanded results root. Expected the benchmark "
            f"pipeline output at {candidate_results_root}."
        )
    ref_dir = copy_refs(Path(args.ref_root), workspace)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_truth = []
    all_pred = []
    all_metrics = []
    for expr in exprs:
        truth_df, pred_df = analyze_expr(results_root, ref_dir, expr, min_naive=args.min_naive)
        all_truth.append(truth_df)
        all_pred.append(pred_df)
        all_metrics.append(metric_from_eval_details(results_root, expr))

    truth_df = pd.concat(all_truth, ignore_index=True) if all_truth else pd.DataFrame()
    pred_df = pd.concat(all_pred, ignore_index=True) if all_pred else pd.DataFrame()
    metrics = pd.concat(all_metrics, ignore_index=True) if all_metrics else pd.DataFrame()
    summaries = summarize_truth(truth_df)

    truth_df.to_csv(out_dir / "per_truth_call_status.csv", index=False, encoding="utf-8-sig")
    pred_df.to_csv(out_dir / "per_prediction_status.csv", index=False, encoding="utf-8-sig")
    metrics.to_csv(out_dir / "metric_check_from_eval_details.csv", index=False, encoding="utf-8-sig")
    for name, df in summaries.items():
        df.to_csv(out_dir / f"{name}.csv", index=False, encoding="utf-8-sig")
    write_report(out_dir, metrics, truth_df, summaries, min_naive=args.min_naive)
    print(f"Wrote semi-synthetic evidence analysis to {out_dir}")
    print(f"truth rows={len(truth_df)} prediction rows={len(pred_df)} metric rows={len(metrics)}")


if __name__ == "__main__":
    main()
