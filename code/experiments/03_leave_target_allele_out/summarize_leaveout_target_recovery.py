#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize target-level recovery for the leave-allele-out experiment."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

CODE_EXPERIMENT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = CODE_EXPERIMENT_DIR.parents[2]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / CODE_EXPERIMENT_DIR.name
sys.path.insert(0, str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))

from common_pantcr_io import (
    coverage_stratum,
    evidence_coverage_by_gene,
    load_pmtr_sequences,
    load_trb_index,
    prediction_records,
    read_table_auto,
    trim_sequence,
)


DEFAULT_METHODS = ["MiXCR", "FindAlleles", "PanTCRLeaveout", "BayesNoPrior"]


def gene_type_from_gene(gene: object) -> str:
    text = str(gene)
    if "TRBV" in text:
        return "V"
    if "TRBJ" in text:
        return "J"
    return ""


def split_words(value: str | Iterable[str]) -> list[str]:
    if isinstance(value, str):
        return [x for x in re.split(r"[\s,]+", value.strip()) if x]
    return [str(x) for x in value]


def parse_int_list(value: object) -> list[int]:
    if value is None or pd.isna(value):
        return []
    out: list[int] = []
    for token in str(value).replace(",", ";").split(";"):
        token = token.strip()
        if not token:
            continue
        out.append(int(token))
    return out


def exact_sequence_match(x: str, y: str) -> bool:
    return str(x).strip().upper() == str(y).strip().upper()


def target_sequence_present_in_catalog(catalog: pd.DataFrame, target_gene: str, target_eval_seq: str) -> bool:
    if catalog.empty or "Family" not in catalog.columns or "Sequence" not in catalog.columns:
        return False
    sub = catalog[catalog["Family"].astype(str) == str(target_gene)]
    target = str(target_eval_seq).strip().upper()
    return any(str(seq).strip().upper() == target for seq in sub["Sequence"].dropna())


def covered_target_defining_positions(defining_positions: list[int], covered_positions: set[int]) -> list[int]:
    return [p for p in defining_positions if p in covered_positions]


def prediction_matches_covered_defining_sites(
    pred_seq: str,
    target_seq: str,
    covered_defining: Iterable[int],
) -> bool:
    pred = str(pred_seq).strip().upper()
    target = str(target_seq).strip().upper()
    positions = list(covered_defining)
    if not positions:
        return False
    for pos in positions:
        if pos < 0 or pos >= len(pred) or pos >= len(target):
            return False
        if pred[pos] != target[pos]:
            return False
    return True


def covered_defining_match_counts(
    target_gene_predictions: list[dict],
    target_eval_seq: str,
    covered_defining_positions: Iterable[int],
) -> tuple[int, int, bool]:
    """Return best matched/conflicting counts among target-gene predictions.

    The pair is selected by maximizing matched target-defining sites and then
    minimizing conflicts for that same matched count. This keeps the auxiliary
    fields aligned with the prediction that provides the strongest target
    compatibility signal.
    """
    positions = list(covered_defining_positions)
    if not positions or not target_gene_predictions:
        return 0, 0, False

    target = str(target_eval_seq).strip().upper()
    best_match = 0
    best_conflict = len(positions)
    single_site_compatible = False

    for pred in target_gene_predictions:
        pred_seq = str(pred.get("seq", "")).strip().upper()
        matched = 0
        conflicting = 0
        for pos in positions:
            if pos < 0 or pos >= len(pred_seq) or pos >= len(target):
                conflicting += 1
            elif pred_seq[pos] == target[pos]:
                matched += 1
            else:
                conflicting += 1
        if len(positions) == 1 and matched == 1 and conflicting == 0:
            single_site_compatible = True
        if matched > best_match or (matched == best_match and conflicting < best_conflict):
            best_match = matched
            best_conflict = conflicting

    return best_match, best_conflict, single_site_compatible


def sequence_in_list(seq: str, seqs: Iterable[str]) -> bool:
    seq_norm = str(seq).strip().upper()
    return any(exact_sequence_match(seq_norm, x) for x in seqs)


def prediction_matches_observed_positions(pred_seq: str, target_seq: str, observed_positions: Iterable[int]) -> bool:
    positions = list(observed_positions)
    if not positions:
        return False
    pred = str(pred_seq).strip().upper()
    target = str(target_seq).strip().upper()
    for pos in positions:
        if pos < 0 or pos >= len(pred) or pos >= len(target):
            return False
        if pred[pos] != target[pos]:
            return False
    return True


def classify_recovery_status(
    target_gene_predictions: list[dict],
    target_eval_seq: str,
    other_truth_eval_seqs: list[str],
    covered_defining_positions: set[int],
    defining_positions: list[int],
    observed_region_positions: set[int],
    cross_gene_error: bool,
    known_same_gene_non_target_eval_seqs: Iterable[str] = (),
) -> str:
    if not target_gene_predictions:
        return "no_call"

    for pred in target_gene_predictions:
        if exact_sequence_match(pred.get("seq", ""), target_eval_seq):
            return "exact_recovery"

    for pred in target_gene_predictions:
        if prediction_matches_observed_positions(pred.get("seq", ""), target_eval_seq, observed_region_positions):
            return "partial_recovery"

    if cross_gene_error:
        return "cross_gene_or_paralogous_error"

    for pred in target_gene_predictions:
        pred_seq = str(pred.get("seq", "")).strip().upper()
        if sequence_in_list(pred_seq, other_truth_eval_seqs):
            return "target_missed"
        if sequence_in_list(pred_seq, known_same_gene_non_target_eval_seqs):
            return "target_missed"

    matched_count, conflicting_count, _ = covered_defining_match_counts(
        target_gene_predictions,
        target_eval_seq,
        covered_target_defining_positions(defining_positions, covered_defining_positions),
    )
    covered_count = len(covered_target_defining_positions(defining_positions, covered_defining_positions))
    return "wrong_same_gene"


def sequence_identity(x: str, y: str) -> float:
    x = str(x).strip().upper()
    y = str(y).strip().upper()
    n = min(len(x), len(y))
    if n <= 0:
        return 0.0
    return sum(1 for i in range(n) if x[i] == y[i]) / n


def build_reference_catalog(pmtr_sequences: dict[str, str], index_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for allele, seq in pmtr_sequences.items():
        gene = str(allele).split("*", 1)[0]
        try:
            gene_type = gene_type_from_gene(gene)
        except ValueError:
            continue
        rows.append({"allele": allele, "gene": gene, "gene_type": gene_type, "seq": trim_sequence(seq, gene, gene_type, index_df)})
    return pd.DataFrame(rows)


def is_cross_gene_or_paralogous_error(
    predictions: list[dict],
    target_gene: str,
    catalog: pd.DataFrame,
    delta: float,
) -> bool:
    if catalog.empty or not predictions:
        return False
    for pred in predictions:
        seq = str(pred.get("seq", "")).strip().upper()
        if not seq:
            continue
        scored = catalog.copy()
        scored["identity"] = scored["seq"].map(lambda ref: sequence_identity(seq, ref))
        best_any = scored.sort_values(["identity", "allele"], ascending=[False, True]).iloc[0]
        same_gene = scored[scored["gene"].astype(str) == str(target_gene)]
        best_same_score = 0.0 if same_gene.empty else float(same_gene["identity"].max())
        if str(best_any["gene"]) != str(target_gene) and float(best_any["identity"]) > best_same_score + delta:
            return True
    return False


def sum_numeric_column(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns or df.empty:
        return 0.0
    return float(pd.to_numeric(df[column], errors="coerce").fillna(0.0).sum())


def sum_split_numeric_column(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns or df.empty:
        return 0.0
    total = 0.0
    for value in df[column].dropna().astype(str):
        for token in re.split(r"[|;,]", value):
            token = token.strip()
            if not token:
                continue
            try:
                total += float(token)
            except ValueError:
                continue
    return total


def target_gene_read_count(evidence_csv: Path, target_gene: str, min_naive: float | None = None) -> tuple[float, int]:
    if not evidence_csv.exists():
        return 0.0, 0
    df = pd.read_csv(evidence_csv)
    if "Family" not in df.columns:
        return 0.0, 0
    if min_naive is not None:
        if "NaiveDiversityIndex" not in df.columns:
            raise ValueError(f"Mutation evidence CSV missing NaiveDiversityIndex for min_naive filtering: {evidence_csv}")
        naive = pd.to_numeric(df["NaiveDiversityIndex"], errors="coerce").fillna(0)
        df = df[naive >= float(min_naive)].copy()
    sub = df[df["Family"].astype(str) == str(target_gene)]
    count = sum_numeric_column(sub, "CloneCount")
    if count <= 0:
        count = sum_split_numeric_column(sub, "CloneCountSplit")
    return count, int(len(sub))


def target_repertoire_support(repertoire_csv: Path, target_gene: str, target_allele: str) -> tuple[int, float]:
    if not repertoire_csv.exists():
        return 0, 0.0
    df = pd.read_csv(repertoire_csv)
    gene_type = gene_type_from_gene(target_gene).lower()
    gene_col = f"{gene_type}_gene"
    allele_col = f"{gene_type}_allele"
    if gene_col not in df.columns or allele_col not in df.columns:
        return 0, 0.0
    sub = df[(df[gene_col].astype(str) == str(target_gene)) & (df[allele_col].astype(str) == str(target_allele))]
    return int(len(sub)), sum_numeric_column(sub, "read_count")


def find_inference_csv(run_root: Path, expr: str, method: str, pop: str, gene: str, seed: int) -> Path:
    return run_root / "results" / "infer" / expr / method / pop / gene / f"infer_{pop}_seed{seed}.{gene}.csv"


def find_evidence_csv(run_root: Path, expr: str, pop: str, gene: str, seed: int) -> Path:
    return run_root / "results" / "mutations" / expr / pop / gene / f"{pop}_seed{seed}.{gene}_sequences.csv"


def find_repertoire_csv(run_root: Path, expr: str, pop: str, seed: int) -> Path:
    return run_root / "samples" / expr / pop / f"seed{seed}" / f"{pop}_seed{seed}_repertoire.csv"


def resolve_genotype_csv(run_root: Path, expr: str, pop: str, seed: int, manifest_path: Path) -> Path:
    if manifest_path.exists():
        return manifest_path
    local_label = run_root / "results" / "labels" / expr / pop / f"genotype_{pop}_seed{seed}.csv"
    if local_label.exists():
        return local_label
    return manifest_path


def other_truth_sequences(
    genotype_csv: Path,
    target_gene: str,
    target_allele: str,
    target_gene_type: str,
    pmtr_sequences: dict[str, str],
    index_df: pd.DataFrame,
) -> list[str]:
    df = pd.read_csv(genotype_csv)
    sub = df[df["gene"].astype(str) == str(target_gene)]
    if sub.empty:
        return []
    row = sub.iloc[0]
    out: list[str] = []
    for slot in ["A", "B"]:
        allele = str(row.get(f"allele_{slot}", "")).strip()
        seq = str(row.get(f"seq_{slot}", "")).strip()
        if allele == target_allele:
            continue
        if allele in pmtr_sequences:
            eval_seq = trim_sequence(pmtr_sequences[allele], target_gene, target_gene_type, index_df)
        else:
            eval_seq = trim_sequence(seq, target_gene, target_gene_type, index_df)
        if eval_seq:
            out.append(str(eval_seq).upper())
    return list(dict.fromkeys(out))


def known_same_gene_non_target_sequences(
    reference_catalog: pd.DataFrame,
    target_gene: str,
    target_eval_seq: str,
) -> list[str]:
    if reference_catalog.empty:
        return []
    sub = reference_catalog[reference_catalog["gene"].astype(str) == str(target_gene)]
    target = str(target_eval_seq).strip().upper()
    out = []
    for seq in sub["seq"].dropna().astype(str):
        seq_norm = seq.strip().upper()
        if seq_norm and seq_norm != target:
            out.append(seq_norm)
    return list(dict.fromkeys(out))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", default="runs")
    parser.add_argument("--expr", default="expr_leaveout_allele_specific")
    parser.add_argument("--test-pop", default="PANEL")
    parser.add_argument("--methods", default=" ".join(DEFAULT_METHODS))
    parser.add_argument("--target-manifest", default="")
    parser.add_argument("--target-metadata", default="")
    parser.add_argument("--graph-dir", default="")
    parser.add_argument("--pmtr-ref", default="ref/pmTR_TRB_V_J_cleaned.csv")
    parser.add_argument("--index", default="ref/TRB_index.csv")
    parser.add_argument("--out-dir", default="")
    parser.add_argument("--cross-gene-delta", type=float, default=0.0)
    parser.add_argument("--min-naive", type=float, default=2)
    args = parser.parse_args()

    run_root_arg = Path(args.run_root)
    run_root = run_root_arg.resolve() if run_root_arg.is_absolute() else (Path.cwd() / run_root_arg).resolve()
    manifest_base = run_root / "results" / "leaveout" / args.expr
    target_manifest = Path(args.target_manifest) if args.target_manifest else manifest_base / "target_panel_manifest.tsv"
    target_metadata_path = Path(args.target_metadata) if args.target_metadata else manifest_base / "target_metadata.csv"
    graph_dir_override = Path(args.graph_dir) if args.graph_dir else None
    pmtr_ref_path = Path(args.pmtr_ref)
    if not pmtr_ref_path.is_absolute():
        cwd_candidate = (Path.cwd() / pmtr_ref_path).resolve()
        pmtr_ref_path = cwd_candidate if cwd_candidate.exists() else run_root / args.pmtr_ref
    index_path = Path(args.index)
    if not index_path.is_absolute():
        cwd_candidate = (Path.cwd() / index_path).resolve()
        index_path = cwd_candidate if cwd_candidate.exists() else run_root / args.index
    out_dir = Path(args.out_dir) if args.out_dir else manifest_base / "summary"
    if not out_dir.is_absolute():
        out_dir = (Path.cwd() / out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = read_table_auto(target_manifest)
    metadata = pd.read_csv(target_metadata_path)
    metadata_by_target = {str(r["target_allele"]): r for _, r in metadata.iterrows()}
    methods = split_words(args.methods)
    pmtr_sequences = load_pmtr_sequences(pmtr_ref_path)
    index_df = load_trb_index(index_path)
    reference_catalog = build_reference_catalog(pmtr_sequences, index_df)

    graph_catalog_cache: dict[str, pd.DataFrame] = {}

    rows: list[dict[str, object]] = []
    for _, mrow in manifest.iterrows():
        target_allele = str(mrow["target_allele"])
        target_gene = str(mrow["target_gene"])
        meta = metadata_by_target[target_allele]
        target_gene_type = str(meta.get("target_gene_type", "")).strip() or gene_type_from_gene(target_gene)
        seed = int(mrow["seed"])
        pop = str(mrow.get("pop", args.test_pop))
        genotype_csv = resolve_genotype_csv(run_root, args.expr, pop, seed, Path(str(mrow["genotype_csv"])))
        target_eval_seq = str(meta["target_eval_seq"]).upper()
        defining = parse_int_list(meta.get("defining_positions", ""))
        evidence_csv = find_evidence_csv(run_root, args.expr, pop, target_gene_type, seed)
        repertoire_csv = find_repertoire_csv(run_root, args.expr, pop, seed)
        coverage_by_gene = evidence_coverage_by_gene(evidence_csv, min_naive=args.min_naive)
        observed_region_by_gene = evidence_coverage_by_gene(evidence_csv, min_naive=None)
        covered_positions = coverage_by_gene.get(target_gene, set())
        observed_region_positions = {
            pos for pos in observed_region_by_gene.get(target_gene, set()) if 0 <= pos < len(target_eval_seq)
        }
        covered_defining = covered_target_defining_positions(defining, covered_positions)
        stratum = coverage_stratum(defining, covered_positions)
        read_count, evidence_count = target_gene_read_count(evidence_csv, target_gene, min_naive=args.min_naive)
        repertoire_file_available = repertoire_csv.exists()
        if repertoire_file_available:
            clonotype_count, repertoire_read_count = target_repertoire_support(repertoire_csv, target_gene, target_allele)
        else:
            clonotype_count, repertoire_read_count = pd.NA, pd.NA
        other_truth = other_truth_sequences(genotype_csv, target_gene, target_allele, target_gene_type, pmtr_sequences, index_df)
        known_same_gene_non_target = known_same_gene_non_target_sequences(reference_catalog, target_gene, target_eval_seq)
        if graph_dir_override is not None:
            graph_dir = graph_dir_override
        else:
            graph_dir = run_root / "results" / "pang" / args.expr / "LeaveoutGlobal" / target_gene_type
        graph_key = str(graph_dir)
        if graph_key not in graph_catalog_cache:
            catalog_path = graph_dir / "allele_catalog.csv"
            if not catalog_path.exists():
                raise FileNotFoundError(f"Missing expected leave-target graph allele catalog: {catalog_path}")
            graph_catalog_cache[graph_key] = pd.read_csv(catalog_path)
        graph_catalog = graph_catalog_cache[graph_key]
        graph_catalog_available = not graph_catalog.empty
        sequence_present = target_sequence_present_in_catalog(graph_catalog, target_gene, target_eval_seq)

        for method in methods:
            infer_csv = find_inference_csv(run_root, args.expr, method, pop, target_gene_type, seed)
            if not infer_csv.exists():
                raise FileNotFoundError(f"Missing expected inference CSV for {method} {target_allele}: {infer_csv}")
            preds = prediction_records(infer_csv, target_gene_type, index_df)
            target_preds = [p for p in preds if str(p.get("gene", "")) == target_gene]
            cross_gene = is_cross_gene_or_paralogous_error(
                target_preds,
                target_gene,
                reference_catalog,
                args.cross_gene_delta,
            )
            matched_count, conflicting_count, single_site_compatible = covered_defining_match_counts(
                target_preds,
                target_eval_seq,
                covered_defining,
            )
            status = classify_recovery_status(
                target_gene_predictions=target_preds,
                target_eval_seq=target_eval_seq,
                other_truth_eval_seqs=other_truth,
                covered_defining_positions=set(covered_defining),
                defining_positions=defining,
                observed_region_positions=set(observed_region_positions),
                cross_gene_error=cross_gene,
                known_same_gene_non_target_eval_seqs=known_same_gene_non_target,
            )
            rows.append(
                {
                    "method": method,
                    "target_allele": target_allele,
                    "target_gene": target_gene,
                    "target_gene_type": target_gene_type,
                    "cross_gene_assessment_scope": "effective target-gene-assigned predictions only; other same-chain gene calls do not trigger target-level cross-gene error",
                    "cross_gene_delta": args.cross_gene_delta,
                    "seed": seed,
                    "replicate": mrow.get("replicate", ""),
                    "pop": pop,
                    "status": status,
                    "exact_recovery": int(status == "exact_recovery"),
                    "partial_recovery": int(status == "partial_recovery"),
                    "target_missed": int(status == "target_missed"),
                    "wrong_same_gene": int(status == "wrong_same_gene"),
                    "cross_gene_or_paralogous_error": int(status == "cross_gene_or_paralogous_error"),
                    "no_call": int(status == "no_call"),
                    "evaluable_in_mature_region": bool(meta["evaluable_in_mature_region"]),
                    "graph_catalog_available": bool(graph_catalog_available),
                    "target_sequence_present_in_prior_graph": bool(sequence_present),
                    "n_defining_variants_in_evaluation_region": int(meta["n_defining_variants_in_evaluation_region"]),
                    "defining_positions": meta.get("defining_positions", ""),
                    "first_defining_position": meta.get("first_defining_position", ""),
                    "last_defining_position": meta.get("last_defining_position", ""),
                    "variant_position_class": meta.get("variant_position_class", ""),
                    "pmtr_level": meta.get("pmtr_level", ""),
                    "pmtr_is_new": meta.get("pmtr_is_new", ""),
                    "pmtr_is_default": meta.get("pmtr_is_default", ""),
                    "pmtr_is_new_noCDR3": meta.get("pmtr_is_new_noCDR3", ""),
                    "pmtr_is_default_noCDR3": meta.get("pmtr_is_default_noCDR3", ""),
                    "pmtr_mutations": meta.get("pmtr_mutations", ""),
                    "max_population_count": meta.get("max_population_count", ""),
                    "max_population": meta.get("max_population", ""),
                    "coverage_stratum": stratum,
                    "covered_defining_variant_count": len(covered_defining),
                    "matched_covered_defining_variant_count": matched_count,
                    "conflicting_covered_defining_variant_count": conflicting_count,
                    "single_site_compatible": bool(single_site_compatible),
                    "evidence_summary_min_naive_filter": args.min_naive,
                    "partial_status_definition": (
                        "observed-region compatible partial status: the complete inferred sequence is not an exact "
                        "target match, but at least one same-gene prediction is concordant with the target across all "
                        "MiXCR-derived observed positions for the target gene before NaiveDiversityIndex filtering."
                    ),
                    "total_defining_variant_count": len(defining),
                    "target_gene_read_count": read_count,
                    "target_gene_evidence_item_count": evidence_count,
                    "repertoire_file_available": bool(repertoire_file_available),
                    "target_allele_clonotype_count": clonotype_count,
                    "target_allele_repertoire_read_count": repertoire_read_count,
                    "infer_csv": str(infer_csv),
                    "evidence_csv": str(evidence_csv),
                    "genotype_csv": str(genotype_csv),
                }
            )

    per_sample = pd.DataFrame(rows)
    for numeric_col in ["target_allele_clonotype_count", "target_allele_repertoire_read_count"]:
        if numeric_col in per_sample.columns:
            per_sample[numeric_col] = pd.to_numeric(per_sample[numeric_col], errors="coerce")
    observed_methods = set(per_sample["method"].dropna().astype(str)) if not per_sample.empty else set()
    expected_methods = set(methods)
    if observed_methods != expected_methods:
        raise RuntimeError(f"Unexpected method set in summary: observed={sorted(observed_methods)}, expected={sorted(expected_methods)}")
    expected_rows = len(manifest) * len(methods)
    if len(per_sample) != expected_rows:
        raise RuntimeError(f"Unexpected summary row count: observed={len(per_sample)}, expected={expected_rows}")
    per_sample_path = out_dir / "per_target_method_status.csv"
    per_sample.to_csv(per_sample_path, index=False)

    status_order = [
        "exact_recovery",
        "partial_recovery",
        "target_missed",
        "wrong_same_gene",
        "cross_gene_or_paralogous_error",
        "no_call",
    ]
    group_cols = [
        "method",
        "target_allele",
        "target_gene",
        "target_gene_type",
        "evaluable_in_mature_region",
        "graph_catalog_available",
        "target_sequence_present_in_prior_graph",
        "coverage_stratum",
    ]
    summary = per_sample.groupby(group_cols, dropna=False).agg(
        n_replicates=("seed", "count"),
        exact_recovery_count=("exact_recovery", "sum"),
        partial_recovery_count=("partial_recovery", "sum"),
        target_missed_count=("target_missed", "sum"),
        wrong_same_gene_count=("wrong_same_gene", "sum"),
        cross_gene_or_paralogous_error_count=("cross_gene_or_paralogous_error", "sum"),
        no_call_count=("no_call", "sum"),
        mean_target_gene_read_count=("target_gene_read_count", "mean"),
        mean_target_gene_evidence_item_count=("target_gene_evidence_item_count", "mean"),
        mean_target_allele_clonotype_count=("target_allele_clonotype_count", "mean"),
        mean_covered_defining_variant_count=("covered_defining_variant_count", "mean"),
        mean_matched_covered_defining_variant_count=("matched_covered_defining_variant_count", "mean"),
        mean_conflicting_covered_defining_variant_count=("conflicting_covered_defining_variant_count", "mean"),
        single_site_compatible_count=("single_site_compatible", "sum"),
        total_defining_variant_count=("total_defining_variant_count", "max"),
        pmtr_level=("pmtr_level", "first"),
        pmtr_is_new=("pmtr_is_new", "first"),
        pmtr_is_default=("pmtr_is_default", "first"),
        pmtr_is_new_noCDR3=("pmtr_is_new_noCDR3", "first"),
        pmtr_is_default_noCDR3=("pmtr_is_default_noCDR3", "first"),
        pmtr_mutations=("pmtr_mutations", "first"),
        max_population_count=("max_population_count", "first"),
        max_population=("max_population", "first"),
        evidence_summary_min_naive_filter=("evidence_summary_min_naive_filter", "first"),
        cross_gene_delta=("cross_gene_delta", "first"),
        partial_status_definition=("partial_status_definition", "first"),
    ).reset_index()
    for status in status_order:
        count_col = f"{status}_count"
        rate_col = f"{status}_rate"
        summary[rate_col] = summary[count_col] / summary["n_replicates"]
    summary["single_site_compatible_rate"] = summary["single_site_compatible_count"] / summary["n_replicates"]
    summary_path = out_dir / "target_recovery_summary.csv"
    summary.to_csv(summary_path, index=False)

    overall = per_sample.groupby(["method", "status"]).size().reset_index(name="n")
    overall["rate"] = overall.groupby("method")["n"].transform(lambda x: x / x.sum())
    overall["evidence_summary_min_naive_filter"] = args.min_naive
    overall["cross_gene_delta"] = args.cross_gene_delta
    overall["partial_status_definition"] = (
        "observed-region compatible partial status; same-gene non-exact prediction matches the target "
        "across all MiXCR-derived observed positions within the evaluated region"
    )
    overall_path = out_dir / "overall_status_summary.csv"
    overall.to_csv(overall_path, index=False)

    print(f"Wrote per-sample target status: {per_sample_path}")
    print(f"Wrote target recovery summary: {summary_path}")
    print(f"Wrote overall status summary: {overall_path}")


if __name__ == "__main__":
    main()
