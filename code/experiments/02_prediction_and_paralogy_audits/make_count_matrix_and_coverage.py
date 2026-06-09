#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create source count records and covered/imputed defining-variant strata.

Run after the normal simu pipeline has produced results/infer and results/labels.
The matching status produced here is gene-aware and is used for source records,
coverage strata, and diagnostics. Final overall S2/S5 sequence-recovery tables
are rebuilt separately by rebuild_sequence_only_primary_tables.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[3]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / Path(__file__).resolve().parent.name
DEFAULT_RESULTS_ROOT = PACKAGE_ROOT / "experiments" / "01_in_silico_trbv_benchmarks" / "results"
DEFAULT_REF_ROOT = PACKAGE_ROOT / "data" / "ref"
sys.path.append(str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))
from common_pantcr_io import (  # noqa: E402
    compatible_match_score,
    coverage_stratum,
    defining_positions,
    evidence_coverage_by_gene,
    evidence_observed_bases_by_gene,
    evidence_support_counts,
    evidence_support_stratum,
    find_evidence_file,
    find_label_file,
    find_validation_fold,
    load_default_refs,
    load_pmtr_sequences,
    load_trb_index,
    prediction_records,
    sample_id_from_filename,
    sequence_length_difference,
    split_words,
    trim_sequence,
    truth_records,
)


EXPECTED_SEED_RANGES = {
    "expr_ScenarioA": (0, 49),
    "expr_ScenarioB": (0, 49),
    "expr_FullLength": (0, 49),
    "expr_ScenarioC": (50, 99),
}


def resolve_ref_path(ref_root: Path, value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    parts = path.parts
    if parts and parts[0] == "ref":
        path = Path(*parts[1:])
    return ref_root / path


def match_truth_and_predictions(truth, preds, gene_type):
    """Exact-first matching, then longest prefix/suffix-compatible partial."""
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

    for pi, p in enumerate(preds):
        pred_status = "matched" if pi in used_preds else "fp"
        pred_rows.append((pi, pred_status))

    truth_rows = [truth_rows_by_index[ti] for ti in range(len(truth))]
    return truth_rows, pred_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", default=str(DEFAULT_RESULTS_ROOT), help="Generated results root, usually experiments/01.../results.")
    ap.add_argument("--ref-root", default=str(DEFAULT_REF_ROOT), help="Reference directory, usually data/ref.")
    ap.add_argument("--expr", default="expr_ScenarioA")
    ap.add_argument("--methods", default="MiXCR FindAlleles Bayes BayesNoPrior")
    ap.add_argument("--pops", default="AFR EUR AMR EAS SAS")
    ap.add_argument("--genes", default="V J")
    ap.add_argument("--out-dir", default=str(EXPERIMENT_DIR / "generated" / "count_matrix_and_coverage_strata"))
    ap.add_argument("--pmtr", default="pmTR_TRB_V_J_cleaned.csv")
    ap.add_argument("--index", default="TRB_index.csv")
    ap.add_argument("--default-ref", default="IMGT_TRB_default.csv")
    ap.add_argument(
        "--min-naive",
        type=float,
        default=1,
        help="NaiveDiversityIndex threshold for mutation evidence rows; expr_ScenarioA/expr_ScenarioB manuscript simulation recovery used min_naive=1.",
    )
    ap.add_argument(
        "--expected-seeds-per-pop",
        type=int,
        default=50,
        help="Require exactly this many unique seed-level inference files for every method/population/gene combination.",
    )
    ap.add_argument(
        "--expected-seed-start",
        type=int,
        default=None,
        help="Optional exact expected seed range start. If omitted, known manuscript expr defaults are used.",
    )
    ap.add_argument(
        "--expected-seed-end",
        type=int,
        default=None,
        help="Optional exact expected seed range end, inclusive. If omitted, known manuscript expr defaults are used.",
    )
    args = ap.parse_args()

    results_root = Path(args.results_root)
    ref_root = Path(args.ref_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    methods = split_words(args.methods, ["MiXCR", "FindAlleles", "Bayes", "BayesNoPrior"])
    pops = split_words(args.pops, ["AFR", "EUR", "AMR", "EAS", "SAS"])
    genes = split_words(args.genes, ["V", "J"])

    index_df = load_trb_index(resolve_ref_path(ref_root, args.index))
    pmtr_map = load_pmtr_sequences(resolve_ref_path(ref_root, args.pmtr))
    default_refs_raw = load_default_refs(resolve_ref_path(ref_root, args.default_ref))

    truth_out = []
    pred_out = []
    input_manifest = []

    for method in methods:
        for pop in pops:
            for gene_type in genes:
                infer_dir = results_root / "infer" / args.expr / method / pop / gene_type
                if not infer_dir.exists():
                    input_manifest.append({
                        "method": method,
                        "population": pop,
                        "gene_type": gene_type,
                        "infer_dir": str(infer_dir),
                        "infer_file": "",
                        "label_file": "",
                        "evidence_file": "",
                        "seed": "",
                        "status": "missing_infer_dir",
                    })
                    continue
                infer_files = sorted(infer_dir.glob("*.csv"))
                if not infer_files:
                    input_manifest.append({
                        "method": method,
                        "population": pop,
                        "gene_type": gene_type,
                        "infer_dir": str(infer_dir),
                        "infer_file": "",
                        "label_file": "",
                        "evidence_file": "",
                        "seed": "",
                        "status": "missing_infer_files",
                    })
                    continue
                for infer_csv in infer_files:
                    _, seed = sample_id_from_filename(infer_csv)
                    if seed is None:
                        input_manifest.append({
                            "method": method,
                            "population": pop,
                            "gene_type": gene_type,
                            "infer_dir": str(infer_dir),
                            "infer_file": str(infer_csv),
                            "label_file": "",
                            "evidence_file": "",
                            "seed": "",
                            "status": "unparseable_seed_from_infer_filename",
                        })
                        continue
                    label_csv = find_label_file(results_root, args.expr, pop, seed)
                    if label_csv is None:
                        input_manifest.append({
                            "method": method,
                            "population": pop,
                            "gene_type": gene_type,
                            "infer_dir": str(infer_dir),
                            "infer_file": str(infer_csv),
                            "label_file": "",
                            "evidence_file": "",
                            "seed": seed,
                            "status": "missing_label_file",
                        })
                        continue
                    evidence_csv = find_evidence_file(results_root, args.expr, pop, gene_type, seed)
                    if evidence_csv is None:
                        input_manifest.append({
                            "method": method,
                            "population": pop,
                            "gene_type": gene_type,
                            "infer_dir": str(infer_dir),
                            "infer_file": str(infer_csv),
                            "label_file": str(label_csv),
                            "evidence_file": "",
                            "seed": seed,
                            "status": "missing_evidence_file",
                        })
                        continue
                    input_manifest.append({
                        "method": method,
                        "population": pop,
                        "gene_type": gene_type,
                        "infer_dir": str(infer_dir),
                        "infer_file": str(infer_csv),
                        "label_file": str(label_csv),
                        "evidence_file": str(evidence_csv),
                        "seed": seed,
                        "status": "ok",
                    })
                    coverage = evidence_coverage_by_gene(evidence_csv, min_naive=args.min_naive)
                    observed_bases = evidence_observed_bases_by_gene(evidence_csv, min_naive=args.min_naive)
                    fold = find_validation_fold(results_root, args.expr, gene_type, evidence_csv.name)

                    truth = truth_records(label_csv, gene_type, index_df, pmtr_map)
                    preds = prediction_records(infer_csv, gene_type, index_df)
                    matched_truth, matched_preds = match_truth_and_predictions(truth, preds, gene_type)

                    for ti, pi, status, compatible_candidate_count, compatible_match_len, compatible_best_tie_count in matched_truth:
                        t = truth[ti]
                        gene = t["gene"]
                        ref = trim_sequence(default_refs_raw.get(gene, ""), gene, gene_type, index_df)
                        def_pos = defining_positions(t["seq"], ref) if ref else []
                        length_diff = sequence_length_difference(t["seq"], ref) if ref else 0
                        covered = coverage.get(gene, set())
                        stratum = coverage_stratum(def_pos, covered, length_difference=length_diff)
                        support_counts = evidence_support_counts(t["seq"], def_pos, observed_bases.get(gene, {}))
                        support_stratum = evidence_support_stratum(
                            def_pos,
                            covered,
                            support_counts,
                            length_difference=length_diff,
                        )
                        pred_seq = preds[pi]["seq"] if pi is not None else ""
                        truth_out.append({
                            "expr": args.expr,
                            "min_naive_filter": args.min_naive,
                            "evidence_source_level": "PanTCR mutation-detection table after min_naive filtering",
                            "truth_unit": "unique allele sequence per gene per sample in evaluation scope",
                            "partial_status_definition": "prefix/suffix sequence-compatible partial, not conservative defining-site recovery",
                            "method": method,
                            "population": pop,
                            "seed": seed,
                            "fold": fold,
                            "gene_type": gene_type,
                            "gene": gene,
                            "truth_allele": t["allele"],
                            "truth_seq": t["seq"],
                            "truth_seq_full_catalog": t.get("seq_full_catalog", t["seq"]),
                            "truth_sequence_scope": t.get("truth_sequence_scope", "eval_trimmed"),
                            "matched_pred_seq": pred_seq,
                            "compatible_partial_candidate_count": compatible_candidate_count,
                            "compatible_partial_match_length": compatible_match_len,
                            "compatible_partial_best_tie_count": compatible_best_tie_count,
                            "status": status,
                            "coverage_stratum": stratum,
                            "coverage_stratum_definition": "stratified by allele-defining variant positions relative to the default reference",
                            "evidence_support_stratum": support_stratum,
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
                            "evidence_file": str(evidence_csv) if evidence_csv else "",
                            "infer_file": str(infer_csv),
                            "label_file": str(label_csv),
                        })

                    for pi, status in matched_preds:
                        p = preds[pi]
                        pred_out.append({
                            "expr": args.expr,
                            "min_naive_filter": args.min_naive,
                            "method": method,
                            "population": pop,
                            "seed": seed,
                            "fold": fold,
                            "gene_type": gene_type,
                            "gene": p["gene"],
                            "pred_allele": p["allele"],
                            "pred_seq": p["seq"],
                            "status": status,
                            "infer_file": str(infer_csv),
                        })

    manifest_df = pd.DataFrame(input_manifest)
    manifest_df.to_csv(out_dir / "input_completeness_manifest.csv", index=False)
    bad_inputs = manifest_df[~manifest_df["status"].eq("ok")] if not manifest_df.empty else pd.DataFrame()
    if not bad_inputs.empty:
        raise RuntimeError(
            "Missing or malformed inputs for count/coverage analysis; see "
            f"{out_dir / 'input_completeness_manifest.csv'}"
        )
    if args.expected_seeds_per_pop:
        if args.expected_seed_start is not None or args.expected_seed_end is not None:
            if args.expected_seed_start is None or args.expected_seed_end is None:
                raise ValueError("--expected-seed-start and --expected-seed-end must be provided together")
            expected_seed_start, expected_seed_end = int(args.expected_seed_start), int(args.expected_seed_end)
        elif args.expr in EXPECTED_SEED_RANGES:
            expected_seed_start, expected_seed_end = EXPECTED_SEED_RANGES[args.expr]
        else:
            expected_seed_start, expected_seed_end = None, None
        seed_counts = (
            manifest_df.groupby(["method", "population", "gene_type"])
            .agg(
                n_rows=("seed", "size"),
                n_unique_seeds=("seed", "nunique"),
            )
            .reset_index()
        )
        bad_seed_counts = seed_counts[
            (seed_counts["n_rows"] != int(args.expected_seeds_per_pop))
            | (seed_counts["n_unique_seeds"] != int(args.expected_seeds_per_pop))
        ]
        duplicate_seeds = (
            manifest_df.groupby(["method", "population", "gene_type", "seed"])
            .size()
            .reset_index(name="n")
            .query("n > 1")
        )
        seed_set_rows = []
        bad_seed_sets = pd.DataFrame()
        if expected_seed_start is not None and expected_seed_end is not None:
            expected_seeds = {f"seed{i}" for i in range(expected_seed_start, expected_seed_end + 1)}
            for keys, sub in manifest_df.groupby(["method", "population", "gene_type"], dropna=False):
                actual = set(sub["seed"].dropna().astype(str))
                missing = sorted(expected_seeds - actual, key=lambda x: int(x.replace("seed", "")))
                unexpected = sorted(actual - expected_seeds, key=lambda x: int(x.replace("seed", "")) if x.replace("seed", "").isdigit() else x)
                seed_set_rows.append({
                    "method": keys[0],
                    "population": keys[1],
                    "gene_type": keys[2],
                    "expected_seed_start": expected_seed_start,
                    "expected_seed_end": expected_seed_end,
                    "expected_seed_count": len(expected_seeds),
                    "actual_seed_count": len(actual),
                    "missing_expected_seeds": ";".join(missing),
                    "unexpected_seeds": ";".join(unexpected),
                    "seed_set_matches_expected": not missing and not unexpected,
                })
            bad_seed_sets = pd.DataFrame(seed_set_rows)
            bad_seed_sets = bad_seed_sets[~bad_seed_sets["seed_set_matches_expected"].astype(bool)]
            pd.DataFrame(seed_set_rows).to_csv(out_dir / "input_expected_seed_set_check.csv", index=False)
        if not bad_seed_counts.empty or not duplicate_seeds.empty or not bad_seed_sets.empty:
            seed_counts.to_csv(out_dir / "input_seed_cardinality_check.csv", index=False)
            duplicate_seeds.to_csv(out_dir / "input_duplicate_seed_check.csv", index=False)
            raise RuntimeError(
                "Unexpected seed cardinality for count/coverage analysis; see "
                f"{out_dir / 'input_seed_cardinality_check.csv'}"
            )

    truth_df = pd.DataFrame(truth_out)
    pred_df = pd.DataFrame(pred_out)
    truth_df.to_csv(out_dir / "per_truth_call_status.csv", index=False)
    pred_df.to_csv(out_dir / "per_prediction_status.csv", index=False)

    if not truth_df.empty:
        counts = (
            truth_df.groupby(["method", "population", "gene_type", "coverage_stratum", "status"])
            .size()
            .reset_index(name="n")
        )
        counts.insert(0, "expr", args.expr)
        counts.insert(1, "min_naive_filter", args.min_naive)
        counts.to_csv(out_dir / "summary_truth_counts_by_coverage.csv", index=False)

        evidence_counts = (
            truth_df.groupby(["method", "population", "gene_type", "evidence_support_stratum", "status"])
            .size()
            .reset_index(name="n")
        )
        evidence_counts.insert(0, "expr", args.expr)
        evidence_counts.insert(1, "min_naive_filter", args.min_naive)
        evidence_counts.to_csv(out_dir / "summary_truth_counts_by_evidence_support.csv", index=False)

        pred_counts = pred_df.groupby(["method", "population", "gene_type", "status"]).size().reset_index(name="n")
        pred_counts.insert(0, "expr", args.expr)
        pred_counts.insert(1, "min_naive_filter", args.min_naive)
        pred_counts.to_csv(out_dir / "summary_prediction_counts.csv", index=False)

        metric_rows = []
        for keys, sub in truth_df.groupby(["method", "population", "gene_type", "coverage_stratum"]):
            method, pop, gene_type, stratum = keys
            exact_tp = int((sub["status"] == "exact_tp").sum())
            partial = int((sub["status"] == "partial_recovery").sum())
            fn = int((sub["status"] == "fn").sum())
            no_call = int((sub["status"] == "no_call").sum())
            pred_sub = pred_df[
                (pred_df["method"] == method)
                & (pred_df["population"] == pop)
                & (pred_df["gene_type"] == gene_type)
            ]
            fp = int((pred_sub["status"] == "fp").sum())
            precision = exact_tp / (exact_tp + fp) if (exact_tp + fp) else 0.0
            recall = exact_tp / (exact_tp + partial + fn + no_call) if (exact_tp + partial + fn + no_call) else 0.0
            fdr = fp / (exact_tp + fp) if (exact_tp + fp) else 0.0
            metric_rows.append({
                "method": method,
                "min_naive_filter": args.min_naive,
                "population": pop,
                "gene_type": gene_type,
                "coverage_stratum": stratum,
                "exact_tp": exact_tp,
                "partial_recovery": partial,
                "fn": fn,
                "no_call": no_call,
                "fp_global_for_method_pop_gene": fp,
                "precision_exact": precision,
                "recall_exact": recall,
                "fdr_exact": fdr,
                "metric_note": "Exact-only metric within this coverage stratum; FP denominator is global for the same method/population/gene.",
            })
        pd.DataFrame(metric_rows).to_csv(out_dir / "summary_metrics_by_coverage.csv", index=False)

        evidence_metric_rows = []
        for keys, sub in truth_df.groupby(["method", "population", "gene_type", "evidence_support_stratum"]):
            method, pop, gene_type, stratum = keys
            exact_tp = int((sub["status"] == "exact_tp").sum())
            partial = int((sub["status"] == "partial_recovery").sum())
            fn = int((sub["status"] == "fn").sum())
            no_call = int((sub["status"] == "no_call").sum())
            recovered = exact_tp + partial
            total = exact_tp + partial + fn + no_call
            evidence_metric_rows.append({
                "method": method,
                "min_naive_filter": args.min_naive,
                "population": pop,
                "gene_type": gene_type,
                "evidence_support_stratum": stratum,
                "exact_tp": exact_tp,
                "partial_recovery": partial,
                "fn": fn,
                "no_call": no_call,
                "recovered_exact_or_partial": recovered,
                "n_truth_alleles": total,
                "recovered_exact_or_partial_rate": recovered / total if total else 0.0,
                "partial_status_definition": "prefix/suffix sequence-compatible partial",
            })
        pd.DataFrame(evidence_metric_rows).to_csv(out_dir / "summary_metrics_by_evidence_support.csv", index=False)

    print(f"Wrote supplemental count/coverage outputs to {out_dir}")


if __name__ == "__main__":
    main()
