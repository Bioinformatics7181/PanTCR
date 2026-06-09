#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Evaluate inferred TRB allele calls against genotype labels.

Default ``sequence`` mode is the common evaluator for PanTCR, MiXCR-default,
FindAlleles, semi-synthetic, pseudo-bulk, and public-bulk benchmark outputs.
The genotype CSV is the source of truth: truth sequences are read directly from
``seq_A`` and ``seq_B`` and trimmed to the mature evaluation interval with
``TRB_index.csv``. No external allele database is used to re-derive truth.

``custom``/``mixcr-all`` mode is for the expanded-reference MiXCR-all baseline,
whose retained output is a MiXCR clones TSV. It evaluates exact allele-name
hits from ``allVHitsWithScore``/``allJHitsWithScore`` against genotype
``allele_A``/``allele_B``. By default it folds only mature-region-equivalent
allele labels that are indistinguishable in the final clean reference:
``TRBV7-7*01``/``TRBV7-7*02`` and ``TRBV15*01``/``TRBV15*05``.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from common_pantcr_io import (
    compatible_match_score,
    load_trb_index,
    prediction_records,
    truth_records,
)


PACKAGE_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INDEX = PACKAGE_ROOT / "data" / "ref" / "TRB_index.csv"

MATURE_EQUIVALENT_LABELS = {
    # These two pairs have identical sequences over the mature TRBV interval
    # used for evaluation, so allele-name evaluation must collapse them to keep
    # MiXCR-all truth denominators comparable with sequence-based evaluators.
    "TRBV7-7*01": "TRBV7-7*01/02",
    "TRBV7-7*02": "TRBV7-7*01/02",
    "TRBV15*01": "TRBV15*01/05",
    "TRBV15*05": "TRBV15*01/05",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        default="sequence",
        choices=["sequence", "custom", "mixcr-all"],
        help="Evaluation mode. Default: sequence.",
    )
    parser.add_argument("--gt", help="Genotype CSV containing gene, allele_A/B, and seq_A/B columns.")
    parser.add_argument("--infer", help="Inference CSV for sequence mode.")
    parser.add_argument("--index", type=Path, default=DEFAULT_INDEX, help=f"TRB index CSV. Default: {DEFAULT_INDEX}")
    parser.add_argument("--tsv", help="MiXCR clones TSV for custom/mixcr-all mode.")
    parser.add_argument("--csv", help="Alias for --gt in custom/mixcr-all mode.")
    parser.add_argument(
        "--gene_type",
        required=True,
        choices=["V", "J", "all"],
        help="Gene type to evaluate. Sequence mode accepts V or J; custom mode also accepts all.",
    )
    parser.add_argument(
        "--intersect",
        action="store_true",
        help="Sequence mode only: evaluate genes present in both genotype labels and predictions.",
    )
    parser.add_argument(
        "--disable-mature-equivalence-fold",
        action="store_true",
        help="Custom/mixcr-all mode only: do not fold mature-region-equivalent allele labels.",
    )
    parser.add_argument("--out_prefix", default="report", help="Output prefix for summary and detail CSV files.")
    return parser.parse_args()


def clean_value(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.upper() in {"", "NAN", "NONE", "NA"} else text


def canonical_allele_label(label: object, fold_mature_equivalent: bool = True) -> str:
    text = clean_value(label)
    if not text:
        return ""
    if fold_mature_equivalent:
        if text.endswith("_mut"):
            base = text[: -len("_mut")]
            return f"{MATURE_EQUIVALENT_LABELS.get(base, base)}_mut"
        return MATURE_EQUIVALENT_LABELS.get(text, text)
    return text


def custom_truth_label(row: pd.Series, slot: str, fold_mature_equivalent: bool) -> str:
    label = canonical_allele_label(row.get(f"allele_{slot}"), fold_mature_equivalent=fold_mature_equivalent)
    if not label:
        return ""

    mutation_token = clean_value(row.get(f"mut_{slot}"))
    if mutation_token and mutation_token.upper() not in {"FALSE", "0", "NO"}:
        safe_token = re.sub(r"[^A-Za-z0-9_.+-]+", "_", mutation_token)
        label = f"{label}_{safe_token}"
    elif label.endswith("_mut"):
        seq = clean_value(row.get(f"seq_{slot}")).upper()
        if seq:
            seq_hash = hashlib.sha1(seq.encode("utf-8")).hexdigest()[:12]
            label = f"{label}_seq{seq_hash}"
    return label


def gene_type_list(gene_type: str) -> list[str]:
    return ["V", "J"] if gene_type == "all" else [gene_type]


def validate_sequence_args(args: argparse.Namespace) -> tuple[Path, Path]:
    if args.gene_type == "all":
        raise ValueError("sequence mode requires --gene_type V or --gene_type J; use separate calls for V and J.")
    if not args.gt:
        raise ValueError("sequence mode requires --gt.")
    if not args.infer:
        raise ValueError("sequence mode requires --infer.")
    if not args.index.exists():
        raise FileNotFoundError(f"TRB index CSV not found: {args.index}")
    return Path(args.gt), Path(args.infer)


def validate_custom_args(args: argparse.Namespace) -> tuple[Path, Path]:
    gt = args.gt or args.csv
    if not gt:
        raise ValueError("custom/mixcr-all mode requires --gt or --csv.")
    if not args.tsv:
        raise ValueError("custom/mixcr-all mode requires --tsv.")
    return Path(gt), Path(args.tsv)


def records_by_sequence(records: list[dict]) -> dict[str, set[str]]:
    grouped: dict[str, set[str]] = defaultdict(set)
    for row in records:
        seq = str(row.get("seq", "")).strip().upper()
        if not seq:
            continue
        label = clean_value(row.get("allele")) or clean_value(row.get("gene"))
        if label:
            grouped[seq].add(label)
    return dict(grouped)


def labels_for_records(records: list[dict]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    for row in records:
        seq = str(row.get("seq", "")).strip().upper()
        gene = clean_value(row.get("gene"))
        if seq and gene:
            out[seq].add(gene)
    return dict(out)


def legacy_compatible_sequence_matching(
    truth_sequences: list[str],
    pred_sequences: list[str],
    gene_type: str,
) -> tuple[set[int], set[int]]:
    """Use the sequence evaluator's compatible-coverage accounting.

    Each unique predicted sequence is TP if it is compatible with at least one
    unique truth sequence. All compatible truth sequences are marked as covered
    for recall/FN accounting. This keeps the consolidated evaluator aligned
    with the label-sequence truth accounting used by the benchmark workflows.
    """
    matched_pred_indices: set[int] = set()
    matched_truth_indices: set[int] = set()

    for pred_idx, pred_seq in enumerate(pred_sequences):
        pred_matched = False
        for truth_idx, truth_seq in enumerate(truth_sequences):
            score = compatible_match_score(truth_seq, pred_seq, gene_type)
            if score > 0:
                pred_matched = True
                matched_truth_indices.add(truth_idx)
        if pred_matched:
            matched_pred_indices.add(pred_idx)

    return matched_pred_indices, matched_truth_indices


def run_sequence_mode(args: argparse.Namespace) -> None:
    gt_path, infer_path = validate_sequence_args(args)
    index_df = load_trb_index(args.index)

    # Truth comes from the label sequence, not a catalog lookup. The helper
    # still applies the same TRB_index trimming.
    truth_rows = truth_records(gt_path, args.gene_type, index_df, pmtr_map={}, prefer_label_sequence=True)
    pred_rows = prediction_records(infer_path, args.gene_type, index_df)

    if args.intersect:
        truth_genes = {row["gene"] for row in truth_rows}
        pred_genes = {row["gene"] for row in pred_rows}
        common_genes = truth_genes & pred_genes
        truth_rows = [row for row in truth_rows if row["gene"] in common_genes]
        pred_rows = [row for row in pred_rows if row["gene"] in common_genes]

    if not truth_rows:
        raise ValueError(f"No genotype truth records remain for gene type {args.gene_type}.")

    truth = records_by_sequence(truth_rows)
    predictions = labels_for_records(pred_rows)
    truth_sequences = list(truth.keys())
    pred_sequences = list(predictions.keys())

    tp_rows: list[dict[str, str]] = []
    fp_rows: list[dict[str, str]] = []
    fn_rows: list[dict[str, str]] = []
    matched_pred_indices, matched_truth_indices = legacy_compatible_sequence_matching(
        truth_sequences, pred_sequences, args.gene_type
    )

    for pred_idx, pred_seq in enumerate(pred_sequences):
        pred_genes = ",".join(sorted(predictions[pred_seq]))
        if pred_idx in matched_pred_indices:
            tp_rows.append({"Sequence": pred_seq, "GeneType": pred_genes, "Status": "TP"})
        else:
            fp_rows.append({"Sequence": pred_seq, "GeneType": pred_genes, "Status": "FP"})

    for idx, truth_seq in enumerate(truth_sequences):
        if idx not in matched_truth_indices:
            truth_labels = ",".join(sorted(truth[truth_seq]))
            fn_rows.append({"Sequence": truth_seq, "GeneType": truth_labels, "Status": "FN"})

    write_sequence_outputs(args, tp_rows, fp_rows, fn_rows, len(truth_sequences), len(pred_sequences))


def write_sequence_outputs(
    args: argparse.Namespace,
    tp_rows: list[dict[str, str]],
    fp_rows: list[dict[str, str]],
    fn_rows: list[dict[str, str]],
    truth_n: int,
    pred_n: int,
) -> None:
    tp_count = len(tp_rows)
    fp_count = len(fp_rows)
    fn_count = len(fn_rows)
    precision = tp_count / pred_n if pred_n else 0.0
    recall = tp_count / truth_n if truth_n else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall else 0.0

    summary = pd.DataFrame(
        [
            {"Metric": "GeneType", "Value": args.gene_type},
            {"Metric": "EvaluationMode", "Value": "sequence"},
            {"Metric": "IntersectOnly", "Value": args.intersect},
            {"Metric": "TruthSource", "Value": "genotype_seq_A_seq_B_trimmed_by_TRB_index"},
            {"Metric": "TRBIndex", "Value": str(args.index)},
            {"Metric": "Precision", "Value": precision},
            {"Metric": "Recall", "Value": recall},
            {"Metric": "F1-score", "Value": f1_score},
            {"Metric": "TP_Count", "Value": tp_count},
            {"Metric": "FP_Count", "Value": fp_count},
            {"Metric": "FN_Count", "Value": fn_count},
        ]
    )
    summary.to_csv(f"{args.out_prefix}_summary.csv", index=False)

    detailed = pd.DataFrame(tp_rows + fp_rows + fn_rows)
    if detailed.empty:
        detailed = pd.DataFrame(columns=["Sequence", "GeneType", "Status"])
    else:
        detailed = detailed[["Sequence", "GeneType", "Status"]]
    detailed.to_csv(f"{args.out_prefix}_detailed.csv", index=False)

    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1-score:  {f1_score:.4f}")
    print(f"Counts:    TP={tp_count} FP={fp_count} FN={fn_count}")


def parse_tsv_hits(df: pd.DataFrame, col_name: str, fold_mature_equivalent: bool) -> set[str]:
    alleles: set[str] = set()
    if col_name not in df.columns:
        print(f"Warning: Column {col_name} not found in TSV.", file=sys.stderr)
        return alleles

    for hit_str in df[col_name]:
        if pd.isna(hit_str):
            continue
        for hit in str(hit_str).split(","):
            match = re.match(r"([^(]+)", hit.strip())
            if not match:
                continue
            label = canonical_allele_label(match.group(1), fold_mature_equivalent=fold_mature_equivalent)
            if label:
                alleles.add(label)
    return alleles


def parse_csv_truth(df: pd.DataFrame, gene_type: str, fold_mature_equivalent: bool) -> set[str]:
    required = {"gene", "allele_A", "allele_B"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Genotype CSV missing required columns for custom mode: {sorted(missing)}")
    pattern = f"TRB{gene_type}"
    alleles: set[str] = set()
    for _, row in df[df["gene"].astype(str).str.contains(pattern, na=False)].iterrows():
        for slot in ["A", "B"]:
            label = custom_truth_label(row, slot, fold_mature_equivalent=fold_mature_equivalent)
            if label:
                alleles.add(label)
    return alleles


def calculate_metrics(inferred_set: set[str], truth_set: set[str]) -> dict[str, object]:
    tp = inferred_set & truth_set
    fp = inferred_set - truth_set
    fn = truth_set - inferred_set
    precision = len(tp) / (len(tp) + len(fp)) if tp or fp else 0.0
    recall = len(tp) / (len(tp) + len(fn)) if tp or fn else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "Precision": precision,
        "Recall": recall,
        "F1-score": f1_score,
        "TP": sorted(tp),
        "FP": sorted(fp),
        "FN": sorted(fn),
    }


def run_custom_mode(args: argparse.Namespace) -> None:
    gt_path, tsv_path = validate_custom_args(args)
    fold_mature_equivalent = not args.disable_mature_equivalence_fold
    tsv_df = pd.read_csv(tsv_path, sep="\t")
    gt_df = pd.read_csv(gt_path)
    col_map = {"V": "allVHitsWithScore", "J": "allJHitsWithScore"}

    summary_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    for gt in gene_type_list(args.gene_type):
        inferred = parse_tsv_hits(tsv_df, col_map[gt], fold_mature_equivalent)
        truth = parse_csv_truth(gt_df, gt, fold_mature_equivalent)
        metrics = calculate_metrics(inferred, truth)
        summary_rows.extend(
            [
                {"GeneType": gt, "Metric": "GeneType", "Value": gt},
                {"GeneType": gt, "Metric": "EvaluationMode", "Value": "custom_allele_name"},
                {"GeneType": gt, "Metric": "TruthSource", "Value": "genotype_allele_A_allele_B_plus_mutation_token_when_present"},
                {
                    "GeneType": gt,
                    "Metric": "MatureEquivalentFold",
                    "Value": "TRBV7-7*01/02;TRBV15*01/05" if fold_mature_equivalent else "disabled",
                },
                {"GeneType": gt, "Metric": "Precision", "Value": metrics["Precision"]},
                {"GeneType": gt, "Metric": "Recall", "Value": metrics["Recall"]},
                {"GeneType": gt, "Metric": "F1-score", "Value": metrics["F1-score"]},
                {"GeneType": gt, "Metric": "TP_Count", "Value": len(metrics["TP"])},
                {"GeneType": gt, "Metric": "FP_Count", "Value": len(metrics["FP"])},
                {"GeneType": gt, "Metric": "FN_Count", "Value": len(metrics["FN"])},
            ]
        )
        for status in ["TP", "FP", "FN"]:
            for allele in metrics[status]:
                detail_rows.append({"GeneType": gt, "Result": status, "Allele": allele})
        print(
            f"{gt}: Precision={metrics['Precision']:.4f} "
            f"Recall={metrics['Recall']:.4f} F1={metrics['F1-score']:.4f} "
            f"TP={len(metrics['TP'])} FP={len(metrics['FP'])} FN={len(metrics['FN'])}"
        )

    summary_df = pd.DataFrame(summary_rows)
    detail_df = pd.DataFrame(detail_rows, columns=["GeneType", "Result", "Allele"])
    summary_df.to_csv(f"{args.out_prefix}_summary.csv", index=False)
    detail_df.to_csv(f"{args.out_prefix}_details.csv", index=False)
    # Also write the standard evaluator suffix so downstream utilities can read
    # custom-mode outputs without a separate filename branch.
    detail_df.rename(columns={"Result": "Status"}).to_csv(f"{args.out_prefix}_detailed.csv", index=False)


def main() -> None:
    args = parse_args()
    if args.mode == "sequence":
        run_sequence_mode(args)
    else:
        run_custom_mode(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
