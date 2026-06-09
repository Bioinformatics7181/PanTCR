#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fine evidence-support analysis for manuscript pseudo-bulk scRNA results."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Dict, List, Tuple

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


METHOD_INFERS = {
    "MiXCR": ("MiXCR-default", "infer_MiXCR.csv"),
    "FindAlleles": ("FindAlleles", "infer_findalleles.csv"),
    "PanTCR.semi": ("PanTCR.semi", "infer_PanTCR.semi.V.csv"),
}


def file_sha256(path: Path) -> str:
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


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


def read_metric_summary(summary_path: Path, allow_missing_zero: bool = False) -> Dict[str, object]:
    if not summary_path.exists():
        if allow_missing_zero:
            return {
                "precision": 0.0,
                "recall": 0.0,
                "f1": 0.0,
                "metric_file_exists": False,
                "metric_required_keys_present": False,
                "metric_missing_reason": "missing summary file; verified empty inference table with no parsed prediction records",
            }
        raise FileNotFoundError(f"Missing metric summary file: {summary_path}")
    df = pd.read_csv(summary_path)
    if not {"Metric", "Value"}.issubset(df.columns):
        raise ValueError(f"Metric summary missing Metric/Value columns: {summary_path}")
    vals = dict(zip(df["Metric"].astype(str), df["Value"]))
    required = {"Precision", "Recall", "F1-score"}
    missing = sorted(required - set(vals))
    if missing:
        raise ValueError(f"Metric summary missing required keys {missing}: {summary_path}")
    return {
        "precision": float(vals.get("Precision", 0.0)),
        "recall": float(vals.get("Recall", 0.0)),
        "f1": float(vals.get("F1-score", 0.0)),
        "metric_file_exists": True,
        "metric_required_keys_present": True,
        "metric_missing_reason": "",
    }


def metric_path(dataset_dir: Path, method: str) -> Path:
    if method == "MiXCR":
        return dataset_dir / "MiXCR-default" / "eval_MiXCR.V_summary.csv"
    if method == "FindAlleles":
        return dataset_dir / "FindAlleles" / "eval_findalleles.V_summary.csv"
    return dataset_dir / "PanTCR.semi" / "eval_PanTCR.semi.V_summary.csv"


def infer_path(dataset_dir: Path, method: str) -> Path:
    subdir, name = METHOD_INFERS[method]
    return dataset_dir / subdir / name


def dataset_id_from_dir(dataset_dir: Path) -> str:
    return dataset_dir.name.split("_", 1)[1]


def validate_sc_mapping(base: Path, mapping: pd.DataFrame) -> tuple[dict[str, Path], pd.DataFrame]:
    expected_sc_ids = {f"SC-{i}" for i in range(1, 9)}
    observed_mapping_ids = set(mapping["SC_ID"].astype(str))
    if observed_mapping_ids != expected_sc_ids:
        raise RuntimeError(
            "Manuscript pseudo-bulk mapping must contain exactly SC-1 through SC-8 after filtering; "
            f"observed={sorted(observed_mapping_ids)}"
        )
    dir_rows = []
    seen: dict[str, list[Path]] = {}
    for p in base.iterdir():
        if p.is_dir() and p.name.startswith("SC-"):
            sc_id = p.name.split("_", 1)[0]
            seen.setdefault(sc_id, []).append(p)
            dir_rows.append({
                "SC_ID": sc_id,
                "parsed_dir_DatasetID": dataset_id_from_dir(p),
                "dataset_dir": str(p),
                "dir_name": p.name,
            })
    duplicate = {k: v for k, v in seen.items() if len(v) > 1}
    if duplicate:
        msg = "; ".join(f"{k}: {', '.join(str(x) for x in v)}" for k, v in duplicate.items())
        raise RuntimeError(f"Duplicate SC directory prefixes found: {msg}")
    manuscript_ids = set(mapping["SC_ID"].astype(str))
    found_ids = set(seen)
    missing = sorted(manuscript_ids - found_ids)
    if missing:
        raise FileNotFoundError(f"Missing manuscript SC result directories: {', '.join(missing)}")
    sc_to_dir = {k: v[0] for k, v in seen.items()}
    manifest = pd.DataFrame(dir_rows)
    manifest = manifest.merge(mapping[["SC_ID", "DatasetID"]].astype(str), on="SC_ID", how="left")
    manifest["included_in_manuscript_SC1_to_SC8"] = manifest["SC_ID"].isin(manuscript_ids)
    return sc_to_dir, manifest


def build_ref_manifest(ref_root: Path, pmtr_source_used: Path) -> pd.DataFrame:
    names = ["TRB_index.csv", "IMGT_TRB_default.csv", "pmTR_TRB_V_J_cleaned.csv"]
    rows = []
    for name in names:
        path = ref_root / name
        rows.append({
            "ref_file": name,
            "path": str(path),
            "exists": path.exists(),
            "sha256": file_sha256(path),
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "used_for_pmtr_lookup": path.resolve() == pmtr_source_used.resolve() if path.exists() and pmtr_source_used.exists() else False,
        })
    return pd.DataFrame(rows)


def analyze(base: Path, ref_root: Path, min_naive: float) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapping = pd.read_csv(base / "scbulk_manuscript_mode_mapping.csv")
    mapping = mapping[mapping["SC_in_manuscript"].astype(str).str.lower() == "yes"].copy()
    sc_to_dir, mapping_manifest = validate_sc_mapping(base, mapping)

    index_df = load_trb_index(ref_root / "TRB_index.csv")
    pmtr_source = ref_root / "pmTR_TRB_V_J_cleaned.csv"
    pmtr_map = load_pmtr_sequences(pmtr_source)
    if not pmtr_map:
        raise FileNotFoundError(f"No usable cleaned pmTR reference found under {ref_root}")
    ref_manifest = build_ref_manifest(ref_root, pmtr_source)
    default_refs = load_default_refs(ref_root / "IMGT_TRB_default.csv")

    truth_rows = []
    pred_rows = []
    metric_rows = []

    for _, map_row in mapping.iterrows():
        sc_id = str(map_row["SC_ID"])
        dataset_id = str(map_row["DatasetID"])
        dataset_dir = sc_to_dir[sc_id]
        dir_dataset_id = dataset_id_from_dir(dataset_dir)
        if dataset_id != dir_dataset_id:
            raise RuntimeError(f"DatasetID mismatch for {sc_id}: mapping={dataset_id}, directory={dir_dataset_id}")
        evidence_csv = dataset_dir / f"{dataset_id}.V_sequences.csv"
        label_csv = dataset_dir / "genotype.csv"
        if not evidence_csv.exists():
            raise FileNotFoundError(f"Missing pseudo-bulk mutation evidence for {sc_id}/{dataset_id}: {evidence_csv}")
        if not label_csv.exists():
            raise FileNotFoundError(f"Missing pseudo-bulk genotype label for {sc_id}/{dataset_id}: {label_csv}")
        coverage = evidence_coverage_by_gene(evidence_csv, min_naive=min_naive)
        observed = evidence_observed_bases_by_gene(evidence_csv, min_naive=min_naive)
        truth = truth_records(label_csv, "V", index_df, pmtr_map, prefer_label_sequence=True)

        for method in METHOD_INFERS:
            infer_csv = infer_path(dataset_dir, method)
            if not infer_csv.exists():
                raise FileNotFoundError(f"Missing pseudo-bulk inference CSV for {sc_id}/{dataset_id}/{method}: {infer_csv}")
            infer_table_rows = validate_infer_schema(infer_csv)
            preds = prediction_records(infer_csv, "V", index_df)
            metrics = read_metric_summary(
                metric_path(dataset_dir, method),
                allow_missing_zero=(infer_table_rows == 0 and len(preds) == 0),
            )
            metric_rows.append({
                "SC_ID": sc_id,
                "DatasetID": dataset_id,
                "method": method,
                **metrics,
                "infer_file": str(infer_csv),
                "infer_table_rows": infer_table_rows,
                "n_prediction_records": len(preds),
                "metric_file": str(metric_path(dataset_dir, method)),
            })
            matched_truth, matched_preds = match_truth_and_predictions(truth, preds, "V")
            for ti, pi, status, compatible_candidate_count, compatible_match_len, compatible_best_tie_count in matched_truth:
                t = truth[ti]
                gene = t["gene"]
                ref = trim_sequence(default_refs.get(gene, ""), gene, "V", index_df)
                if not ref:
                    raise ValueError(f"Missing default reference sequence for {sc_id}/{dataset_id}/V/{gene}")
                def_pos = defining_positions(t["seq"], ref)
                length_diff = sequence_length_difference(t["seq"], ref)
                covered = coverage.get(gene, set())
                support_counts = evidence_support_counts(t["seq"], def_pos, observed.get(gene, {}))
                truth_rows.append({
                    "SC_ID": sc_id,
                    "DatasetID": dataset_id,
                    "truth_unit": "unique allele sequence per gene per sample in evaluation scope",
                    "partial_status_definition": "prefix/suffix sequence-compatible partial, not conservative defining-site recovery",
                    "partial_matching_rule": "best compatible prediction by prefix/suffix overlap length after exact matches",
                    "method": method,
                    "gene_type": "V",
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
                    "infer_file": str(infer_path(dataset_dir, method)),
                    "infer_table_rows": infer_table_rows,
                    "label_file": str(label_csv),
                    "evidence_file": str(evidence_csv),
                })
            for pi, status in matched_preds:
                p = preds[pi]
                pred_rows.append({
                    "SC_ID": sc_id,
                    "DatasetID": dataset_id,
                    "prediction_evidence_annotation_scope": "prediction status only; false-positive prediction site-level evidence is not annotated in this experiment",
                    "method": method,
                    "gene_type": "V",
                    "gene": p["gene"],
                    "pred_allele": p["allele"],
                    "pred_seq": p["seq"],
                    "status": status,
                    "infer_file": str(infer_path(dataset_dir, method)),
                })

    return pd.DataFrame(truth_rows), pd.DataFrame(pred_rows), pd.DataFrame(metric_rows), mapping_manifest, ref_manifest


def summarize(truth_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    out = {}
    df = truth_df.copy()
    df["recovered_exact_or_partial"] = df["status"].isin(["exact_tp", "partial_recovery"])
    out["summary_by_evidence_support"] = (
        df.groupby(["method", "evidence_support_stratum"])
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
    out["summary_by_sc_and_evidence_support"] = (
        df.groupby(["SC_ID", "DatasetID", "method", "evidence_support_stratum"])
        .agg(
            n_truth_alleles=("status", "size"),
            recovered_exact_or_partial_n=("recovered_exact_or_partial", "sum"),
            recovered_exact_or_partial_rate=("recovered_exact_or_partial", "mean"),
        )
        .reset_index()
    )
    out["summary_by_coverage"] = (
        df.groupby(["method", "coverage_stratum"])
        .agg(
            n_truth_alleles=("status", "size"),
            recovered_exact_or_partial_n=("recovered_exact_or_partial", "sum"),
            recovered_exact_or_partial_rate=("recovered_exact_or_partial", "mean"),
        )
        .reset_index()
    )
    partial = df[df["status"].eq("partial_recovery")].copy()
    if partial.empty:
        out["partial_match_diagnostics"] = pd.DataFrame(
            columns=[
                "method",
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
            partial.groupby(["method"])
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


def write_report(out_dir: Path, metrics: pd.DataFrame, summaries: Dict[str, pd.DataFrame], mapping_manifest: pd.DataFrame, ref_manifest: pd.DataFrame, min_naive: float) -> None:
    lines = []
    lines.append("# V-gene pseudo-bulk scRNA truth-recovery evidence-support analysis\n")
    lines.append("## Scope\n")
    lines.append(f"This report covers the manuscript SC-1 to SC-8 pseudo-bulk datasets. It uses the final `PanTCR.semi` results and the corresponding MiXCR/FindAlleles outputs. Only V-gene truth-recovery results are analyzed because the manuscript SC table reports V-gene metrics; no J-gene conclusion is made here. Mutation-evidence rows were filtered with `NaiveDiversityIndex >= {min_naive:g}` to match the real-data inference setting.\n")
    excluded = mapping_manifest[~mapping_manifest["included_in_manuscript_SC1_to_SC8"].astype(bool)]
    if not excluded.empty:
        lines.append("Available non-manuscript pseudo-bulk directories were excluded from this report:\n")
        lines.append(markdown_table(excluded[["SC_ID", "dir_name"]]))
    lines.append("Dataset-directory mapping used for SC-1 to SC-8:\n")
    lines.append(markdown_table(mapping_manifest[mapping_manifest["included_in_manuscript_SC1_to_SC8"].astype(bool)][["SC_ID", "DatasetID", "parsed_dir_DatasetID", "dir_name"]]))
    lines.append("Reference files used for mature-region trimming, defining-site annotation, and pmTR allele lookup:\n")
    lines.append(markdown_table(ref_manifest[["ref_file", "exists", "used_for_pmtr_lookup", "sha256", "size_bytes"]]))
    pmtr_used = ref_manifest[ref_manifest["used_for_pmtr_lookup"].astype(bool)]
    if not pmtr_used.empty:
        lines.append(f"pmTR allele lookup used the canonical cleaned catalog `{pmtr_used.iloc[0]['ref_file']}`.\n")
    lines.append("\n## Local metric extraction\n")
    show = metrics.copy()
    show["precision"] = show["precision"].map(lambda x: f"{x:.4f}")
    show["recall"] = show["recall"].map(lambda x: f"{x:.4f}")
    show["f1"] = show["f1"].map(lambda x: f"{x:.4f}")
    lines.append(markdown_table(show[["SC_ID", "DatasetID", "method", "precision", "recall", "f1", "metric_file_exists", "metric_required_keys_present", "infer_table_rows", "n_prediction_records", "metric_missing_reason"]]))
    lines.append("These values are read from generated evaluation-summary CSV files. Metric files are required to contain Precision/Recall/F1-score keys. If `metric_file_exists=False`, this is accepted only for a verified empty inference output (`infer_table_rows=0` and no parsed prediction records), and the metric is explicitly recorded as a zero-call result rather than silently imputed.\n")
    lines.append("\n## PanTCR truth-allele evidence-supported recovery strata\n")
    ev = summaries["summary_by_evidence_support"]
    pan = ev[ev["method"] == "PanTCR.semi"].copy()
    pan["exact_or_sequence_compatible_partial_rate"] = pan["recovered_exact_or_partial_rate"].map(lambda x: f"{x:.3f}")
    lines.append(markdown_table(pan[["evidence_support_stratum", "n_truth_alleles", "exact_recovery_n", "partial_recovery_n", "recovered_exact_or_partial_n", "exact_or_sequence_compatible_partial_rate", "mixed_covered_defining_variant_count", "conflicting_covered_defining_variant_count", "no_call_n", "fn_n"]]))
    prior = pan[pan["evidence_support_stratum"] == "prior_only_imputation_dependent"]
    if not prior.empty:
        n_truth = int(prior["n_truth_alleles"].sum())
        n_recovered = int(prior["recovered_exact_or_partial_n"].sum())
        lines.append(f"\nGraph-assisted context count: `{n_recovered}/{n_truth}` PanTCR truth-allele recovery contexts in `prior_only_imputation_dependent` are graph-assisted and should not be described as directly observed.\n")
    lines.append("\n## Sequence-compatible partial-match diagnostics\n")
    partial_diag = summaries.get("partial_match_diagnostics", pd.DataFrame()).copy()
    if not partial_diag.empty:
        partial_diag["mean_compatible_partial_candidate_count"] = partial_diag["mean_compatible_partial_candidate_count"].map(lambda x: f"{float(x):.3f}")
    lines.append(markdown_table(partial_diag))
    lines.append("These diagnostics audit the sequence-compatible partial category. Candidate counts, match lengths, and tie counts are reported so partial recovery is not interpreted as exact full-sequence recovery.\n")
    lines.append("\n## Interpretation\n")
    lines.append("- The pseudo-bulk analysis stratifies truth alleles and their PanTCR recovery status by mutation-evidence support when defining sites are physically covered. False-positive prediction sites are not site-annotated in this table. SC-9/SC-10 are outside the manuscript SC-1 to SC-8 scope and are not included here. Per-SC strata are written to `summary_by_sc_and_evidence_support.csv` for dataset-specific checks.\n")
    lines.append("- As in the simulation analysis, `fully_evidence_supported` means all truth-defining variant sites have compatible observed bases in the filtered pseudo-bulk mutation table. If a covered site has both target and conflicting observed bases, it is conservatively counted as mixed/conflicting rather than cleanly supported. `prior_only_imputation_dependent` means the defining sites are not physically covered and should be described as graph-assisted rather than directly observed.\n")
    lines.append("- The primary use is to clarify the boundary between filtered-mutation-evidence-supported allele evidence and graph-guided completion in real, sparse pseudo-bulk transcriptomic data.\n")
    lines.append("\n## Output files\n")
    for name in [
        "per_truth_call_status.csv",
        "per_prediction_status.csv",
        "metric_check.csv",
        "summary_by_evidence_support.csv",
        "summary_by_sc_and_evidence_support.csv",
        "summary_by_coverage.csv",
        "partial_match_diagnostics.csv",
        "scbulk_dataset_mapping_manifest.csv",
        "reference_file_manifest.csv",
    ]:
        lines.append(f"- `{name}`")
    lines.append("  Note: `per_prediction_status.csv` records prediction status only; false-positive prediction site-level evidence is not annotated.")
    (out_dir / "scbulk_fine_evidence_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    default_base = (
        PACKAGE_ROOT
        / "experiments"
        / "01_in_silico_trbv_benchmarks"
        / "benchmark_pipelines"
        / "03_pseudo_bulk_rnaseq"
        / "generated"
        / "per_dataset_results"
    )
    ap.add_argument("--base", default=str(default_base))
    ap.add_argument("--ref-root", default=str(PACKAGE_ROOT / "data" / "ref"))
    ap.add_argument("--out-dir", default=str(EXPERIMENT_DIR / "generated" / "evidence_analysis"))
    ap.add_argument("--min-naive", type=float, default=2)
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    truth_df, pred_df, metrics, mapping_manifest, ref_manifest = analyze(Path(args.base), Path(args.ref_root), min_naive=args.min_naive)
    summaries = summarize(truth_df)
    truth_df.to_csv(out_dir / "per_truth_call_status.csv", index=False, encoding="utf-8-sig")
    pred_df.to_csv(out_dir / "per_prediction_status.csv", index=False, encoding="utf-8-sig")
    metrics.to_csv(out_dir / "metric_check.csv", index=False, encoding="utf-8-sig")
    mapping_manifest.to_csv(out_dir / "scbulk_dataset_mapping_manifest.csv", index=False, encoding="utf-8-sig")
    ref_manifest.to_csv(out_dir / "reference_file_manifest.csv", index=False, encoding="utf-8-sig")
    for name, df in summaries.items():
        df.to_csv(out_dir / f"{name}.csv", index=False, encoding="utf-8-sig")
    write_report(out_dir, metrics, summaries, mapping_manifest, ref_manifest, min_naive=args.min_naive)
    print(f"Wrote scBulk evidence analysis to {out_dir}")
    print(f"truth rows={len(truth_df)} prediction rows={len(pred_df)} metric rows={len(metrics)}")


if __name__ == "__main__":
    main()
