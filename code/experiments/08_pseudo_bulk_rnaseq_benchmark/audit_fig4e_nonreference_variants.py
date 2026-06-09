#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build a Fig. 4E non-reference TRBV allele recovery/evidence audit.

The audit is truth-centric for pseudo-bulk SC-1..SC-8 under PanTCR.semi.
Non-reference means default-discordant within the evaluated mature V-region.
Variant support is classified using the retained PanTCR evidence rows
(NaiveDiversityIndex >= 2), not raw unfiltered MiXCR-observed ranges.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import pandas as pd

CODE_EXPERIMENT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = CODE_EXPERIMENT_DIR.parents[2]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / CODE_EXPERIMENT_DIR.name

sys.path.append(str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))
from common_pantcr_io import (  # noqa: E402
    defining_positions,
    evidence_observed_bases_by_gene,
    load_default_refs,
    load_trb_index,
    sequence_length_difference,
    trim_sequence,
)


def display_path(path: Path | str) -> str:
    p = Path(path)
    if not p.is_absolute():
        return str(p).replace("\\", "/")
    try:
        return str(p.resolve().relative_to(PACKAGE_ROOT)).replace("\\", "/")
    except ValueError:
        return os.path.relpath(p.resolve(), PACKAGE_ROOT).replace("\\", "/")


# These are the high-shift-similarity rows documented during the Fig. 4E
# non-reference audit. They are retained in the strict sheet, but excluded from
# the revised Fig. 4E-facing sheet.
LIKELY_SHIFT_ARTIFACTS = {
    ("SC-2", "TRBV11-1", "TRBV11-1*01"),
    ("SC-2", "TRBV6-9", "TRBV6-9*01"),
    ("SC-2", "TRBV7-3", "TRBV7-3*01"),
    ("SC-3", "TRBV7-4", "TRBV7-4*01"),
    ("SC-4", "TRBV7-4", "TRBV7-4*01"),
    ("SC-6", "TRBV7-4", "TRBV7-4*01"),
    ("SC-8", "TRBV5-3", "TRBV5-3*01"),
}


def token(pos: int, ref_base: str, truth_base: str) -> str:
    return f"S{ref_base}{pos}{truth_base}"


def join_tokens(values: Iterable[str]) -> str:
    vals = [v for v in values if v]
    return "; ".join(vals)


def classify_tokens(
    gene: str,
    truth_seq: str,
    ref_seq: str,
    observed_by_pos: Dict[int, Set[str]],
) -> dict:
    def_pos = defining_positions(truth_seq, ref_seq)
    all_tokens: List[str] = []
    clean_supported: List[str] = []
    mixed_supported: List[str] = []
    conflicting: List[str] = []
    graph_imputed: List[str] = []
    observed_bases_repr: List[str] = []

    for pos in def_pos:
        ref_base = ref_seq[pos]
        truth_base = truth_seq[pos]
        tok = token(pos, ref_base, truth_base)
        all_tokens.append(tok)
        bases = {str(b).upper() for b in observed_by_pos.get(pos, set()) if str(b).strip()}
        if bases:
            observed_bases_repr.append(f"{tok}:{'/'.join(sorted(bases))}")
        if not bases:
            graph_imputed.append(tok)
        elif bases == {truth_base}:
            clean_supported.append(tok)
        elif truth_base in bases:
            mixed_supported.append(tok)
        else:
            conflicting.append(tok)

    return {
        "Default-relative changes": join_tokens(all_tokens),
        "Evidence-supported changes": join_tokens(clean_supported),
        "Mixed observed changes": join_tokens(mixed_supported),
        "Conflicting observed changes": join_tokens(conflicting),
        "Graph-imputed changes": join_tokens(graph_imputed),
        "Observed retained bases at defining sites": join_tokens(observed_bases_repr),
        "No. of default-relative changes": len(all_tokens),
        "No. of evidence-supported changes": len(clean_supported),
        "No. of mixed observed changes": len(mixed_supported),
        "No. of conflicting observed changes": len(conflicting),
        "No. of graph-imputed changes": len(graph_imputed),
        "Support split": support_split_label(clean_supported, mixed_supported, conflicting, graph_imputed),
    }


def support_split_label(clean: List[str], mixed: List[str], conflicting: List[str], imputed: List[str]) -> str:
    if conflicting:
        return "observed conflict"
    if mixed:
        return "mixed observed support"
    if clean and not imputed:
        return "all defining changes evidence-supported"
    if clean and imputed:
        return "evidence-supported plus graph-imputed"
    if imputed and not clean:
        return "graph-imputed defining changes only"
    return "no default-relative changes"


def status_label(status: str) -> str:
    if status == "exact_tp":
        return "Exact recovery"
    if status == "partial_recovery":
        return "Partial recovery"
    if status == "no_call":
        return "No call"
    if status == "fn":
        return "Not recovered"
    return status


def build_audit(
    truth_csv: Path,
    semi_truth_csv: Path,
    ref_root: Path,
    out_dir: Path,
    min_naive: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    index_df = load_trb_index(ref_root / "TRB_index.csv")
    default_refs_raw = load_default_refs(ref_root / "IMGT_TRB_default.csv")
    default_df = pd.read_csv(ref_root / "IMGT_TRB_default.csv")
    default_alleles = {
        str(row["Gene"]).strip(): str(row.get("Selected Allele", "")).strip()
        for _, row in default_df.iterrows()
        if pd.notna(row.get("Gene"))
    }
    default_refs = {
        gene: trim_sequence(seq, gene, "V", index_df)
        for gene, seq in default_refs_raw.items()
    }

    df = pd.read_csv(truth_csv)
    df = df[(df["method"] == "PanTCR.semi") & (df["gene_type"] == "V")].copy()
    df["n_defining_positions"] = pd.to_numeric(df["n_defining_positions"], errors="coerce").fillna(0).astype(int)
    df = df[df["n_defining_positions"] > 0].copy()

    observed_cache: Dict[str, Dict[str, Dict[int, Set[str]]]] = {}
    rows: List[dict] = []
    for _, r in df.iterrows():
        sc_id = str(r["SC_ID"])
        dataset_id = str(r["DatasetID"])
        gene = str(r["gene"])
        truth_allele = str(r["truth_allele"])
        truth_seq = str(r["truth_seq"]).strip().upper()
        evidence_file = Path(str(r["evidence_file"]))
        if not evidence_file.is_absolute():
            evidence_file = Path.cwd() / evidence_file
        cache_key = str(evidence_file)
        if cache_key not in observed_cache:
            observed_cache[cache_key] = evidence_observed_bases_by_gene(evidence_file, min_naive=min_naive)
        observed_by_gene = observed_cache[cache_key]
        ref_seq = default_refs.get(gene, "")
        if not ref_seq:
            raise ValueError(f"Missing default reference for {gene}")
        length_diff = sequence_length_difference(truth_seq, ref_seq)
        variant_info = classify_tokens(gene, truth_seq, ref_seq, observed_by_gene.get(gene, {}))
        artifact_key = (sc_id, gene, truth_allele)
        likely_artifact = artifact_key in LIKELY_SHIFT_ARTIFACTS
        recovered_exact = str(r["status"]) == "exact_tp"
        row = {
            "SC ID": sc_id,
            "Dataset ID": dataset_id,
            "Gene": gene,
            "Truth allele": truth_allele,
            "Compared default allele": default_alleles.get(gene, ""),
            "Compared reference": "MiXCR default allele for the same TRBV gene",
            "PanTCR recovery status": status_label(str(r["status"])),
            "Recovered by PanTCR": "Yes" if recovered_exact else "No",
            "Evidence support stratum": str(r["evidence_support_stratum"]),
            "Included in revised Fig. 4E non-reference set": "No" if likely_artifact else "Yes",
            "Strict default-discordant row": "Yes",
            "Likely coordinate/boundary-shift artifact": "Yes" if likely_artifact else "No",
            "Length difference vs compared reference": length_diff,
            **variant_info,
            "Matched PanTCR sequence": str(r.get("matched_pred_seq", "")),
            "Truth sequence": truth_seq,
            "Evidence file": display_path(evidence_file),
            "Source truth-status file": display_path(truth_csv),
        }
        rows.append(row)

    detail = pd.DataFrame(rows)
    sort_cols = ["SC ID", "Gene", "Truth allele"]
    detail = detail.sort_values(sort_cols).reset_index(drop=True)
    fig4e_detail = detail[detail["Included in revised Fig. 4E non-reference set"] == "Yes"].copy()

    summary = summarize(fig4e_detail, label="Revised Fig. 4E included")
    strict_summary = summarize(detail, label="Strict default-discordant")
    combined_summary, unrecovered_detail = build_combined_recovery_split_summary(fig4e_detail, semi_truth_csv)
    definitions = pd.DataFrame(
        [
            {
                "Term": "Non-reference truth allele",
                "Definition": "A truth TRBV allele whose evaluated mature V-region sequence differs from the MiXCR default/reference allele sequence of the same TRBV gene.",
            },
            {
                "Term": "Revised Fig. 4E included",
                "Definition": "Strict default-discordant truth alleles after excluding documented likely coordinate/boundary-shift artifacts.",
            },
            {
                "Term": "Default-relative change token",
                "Definition": "S<reference base><0-based coordinate><truth base>, where the reference base is from the same-gene MiXCR default allele in the evaluated mature V-region.",
            },
            {
                "Term": "Evidence-supported changes",
                "Definition": "Default-relative changes whose truth base is cleanly present in retained PanTCR sample-level evidence after NaiveDiversityIndex filtering.",
            },
            {
                "Term": "Mixed observed changes",
                "Definition": "Default-relative changes where the retained evidence contains the truth base and at least one conflicting base.",
            },
            {
                "Term": "Conflicting observed changes",
                "Definition": "Default-relative changes where retained evidence exists at the coordinate but does not contain the truth base.",
            },
            {
                "Term": "Graph-imputed changes",
                "Definition": "Default-relative changes without retained PanTCR evidence at the coordinate; if recovered, the truth base at that coordinate was completed by graph-guided inference/output rather than direct retained evidence.",
            },
            {
                "Term": "No. of Non-default Truth Alleles",
                "Definition": "Truth-centric count of unique TRBV allele sequences per gene per sample whose evaluated mature V-region sequence differs from the same-gene MiXCR default allele. This audit count is not intended to replace aggregate TP/FN denominators in the main performance tables.",
            },
            {
                "Term": "Observed-defining-site Supported",
                "Definition": "A PanTCR-exact-recovered non-default truth allele for which every default-relative defining change has the truth base present in retained sample-level evidence. This category includes both clean and mixed observed evidence, because the defining variant is physically represented in retained evidence.",
            },
            {
                "Term": "Observed + Graph-completed",
                "Definition": "A PanTCR-exact-recovered non-default truth allele for which at least one default-relative defining change is supported by retained sample-level evidence and at least one defining change is completed without retained evidence.",
            },
            {
                "Term": "Graph-completed Only",
                "Definition": "A PanTCR-exact-recovered non-default truth allele for which the default-relative defining changes do not have retained sample-level truth-base support and are completed by graph-guided inference/output.",
            },
        ]
    )
    return detail, fig4e_detail, pd.concat([summary, strict_summary], ignore_index=True), combined_summary, unrecovered_detail, definitions


def recovered_evidence_split(row: pd.Series) -> str:
    n = int(pd.to_numeric(row.get("n_defining_positions", 0), errors="coerce") or 0)
    matched = int(pd.to_numeric(row.get("matched_covered_defining_variant_count", 0), errors="coerce") or 0)
    if n <= 0:
        return "Reference-like"
    if matched >= n:
        return "Observed-defining-site supported"
    if matched == 0:
        return "Graph-completed only"
    return "Observed + graph-completed"


def summarize_recovery_split(df: pd.DataFrame, benchmark: str, dataset_order: List[str]) -> pd.DataFrame:
    rows = []
    for dataset in dataset_order:
        sub = df[df["Dataset"] == dataset].copy()
        if sub.empty:
            continue
        total = int(len(sub))
        recovered = sub[sub["Recovered by PanTCR"] == "Yes"].copy()
        if "Recovery split" not in recovered.columns and not recovered.empty:
            recovered["Recovery split"] = recovered.apply(recovered_evidence_split, axis=1)
        counts = recovered["Recovery split"].value_counts().to_dict() if not recovered.empty else {}
        exact = int(len(recovered))
        unrecovered = sub[sub["Recovered by PanTCR"] != "Yes"].copy()
        unrecovered_counts = unrecovered["Unrecovered observation category"].value_counts().to_dict() if not unrecovered.empty and "Unrecovered observation category" in unrecovered.columns else {}
        rows.append(
            {
                "Benchmark": benchmark,
                "Dataset": dataset,
                "No. of Non-default Truth Alleles": total,
                "PanTCR Exact Recovered": exact,
                "PanTCR Recovery Rate": exact / total if total else 0,
                "Observed-defining-site Supported": int(counts.get("Observed-defining-site supported", 0)),
                "Observed + Graph-completed": int(counts.get("Observed + graph-completed", 0)),
                "Graph-completed Only": int(counts.get("Graph-completed only", 0)),
                "Not Recovered": total - exact,
                "Not Recovered: Truth-base Observed": int(unrecovered_counts.get("Truth-base observed at defining site(s)", 0)),
                "Not Recovered: Covered Without Truth-base": int(unrecovered_counts.get("Covered defining site(s), truth base not retained", 0)),
                "Not Recovered: No Retained Defining-site Coverage": int(unrecovered_counts.get("No retained defining-site coverage", 0)),
            }
        )
    return pd.DataFrame(rows)


def unrecovered_observation_category(n: int, covered: int, matched: int) -> str:
    if matched > 0:
        return "Truth-base observed at defining site(s)"
    if covered > 0:
        return "Covered defining site(s), truth base not retained"
    return "No retained defining-site coverage"


def build_combined_recovery_split_summary(fig4e_detail: pd.DataFrame, semi_truth_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    parts = []
    unrecovered_parts = []

    if semi_truth_csv.exists():
        semi = pd.read_csv(semi_truth_csv)
        semi = semi[
            (semi["method"] == "Bayes")
            & (semi["gene_type"] == "V")
            & (pd.to_numeric(semi["n_defining_positions"], errors="coerce").fillna(0) > 0)
        ].copy()
        for col in [
            "n_defining_positions",
            "covered_defining_variant_count",
            "matched_covered_defining_variant_count",
        ]:
            semi[col] = pd.to_numeric(semi[col], errors="coerce").fillna(0).astype(int)
        semi["Dataset"] = semi["dataset"].astype(str)
        semi["Recovered by PanTCR"] = semi["status"].map(lambda x: "Yes" if str(x) == "exact_tp" else "No")
        semi["Recovery split"] = semi.apply(recovered_evidence_split, axis=1)
        semi["Unrecovered observation category"] = semi.apply(
            lambda r: unrecovered_observation_category(
                int(r["n_defining_positions"]),
                int(r["covered_defining_variant_count"]),
                int(r["matched_covered_defining_variant_count"]),
            ),
            axis=1,
        )
        parts.append(summarize_recovery_split(semi, "Semi-synthetic AIRR-seq", ["AIRR-1", "AIRR-2"]))
        semi_unrec = semi[semi["Recovered by PanTCR"] != "Yes"].copy()
        unrecovered_parts.append(
            semi_unrec.assign(
                Benchmark="Semi-synthetic AIRR-seq",
                Dataset=semi_unrec["Dataset"],
                Sample=semi_unrec["sample_id"].astype(str),
                Gene=semi_unrec["gene"].astype(str),
                **{"Truth allele": semi_unrec["truth_allele"].astype(str)},
                **{"PanTCR status": semi_unrec["status"].map(status_label)},
                **{"No. of default-relative changes": semi_unrec["n_defining_positions"].astype(int)},
                **{"No. of retained covered defining sites": semi_unrec["covered_defining_variant_count"].astype(int)},
                **{"No. of retained truth-base defining sites": semi_unrec["matched_covered_defining_variant_count"].astype(int)},
                **{"Unrecovered observation category": semi_unrec["Unrecovered observation category"]},
                **{"Evidence support stratum": semi_unrec["evidence_support_stratum"].astype(str)},
            )[
                [
                    "Benchmark",
                    "Dataset",
                    "Sample",
                    "Gene",
                    "Truth allele",
                    "PanTCR status",
                    "No. of default-relative changes",
                    "No. of retained covered defining sites",
                    "No. of retained truth-base defining sites",
                    "Unrecovered observation category",
                    "Evidence support stratum",
                ]
            ]
        )

    pseudo = fig4e_detail.copy()
    pseudo["Dataset"] = pseudo["SC ID"].astype(str)
    pseudo["Recovery split"] = ""
    pseudo.loc[
        (pseudo["Recovered by PanTCR"] == "Yes") & (pseudo["No. of evidence-supported changes"].astype(int) + pseudo["No. of mixed observed changes"].astype(int) >= pseudo["No. of default-relative changes"].astype(int)),
        "Recovery split",
    ] = "Observed-defining-site supported"
    pseudo.loc[
        (pseudo["Recovered by PanTCR"] == "Yes") & (pseudo["No. of graph-imputed changes"].astype(int) >= pseudo["No. of default-relative changes"].astype(int)),
        "Recovery split",
    ] = "Graph-completed only"
    pseudo.loc[
        (pseudo["Recovered by PanTCR"] == "Yes") & (pseudo["Recovery split"] == ""),
        "Recovery split",
    ] = "Observed + graph-completed"
    pseudo["Unrecovered observation category"] = pseudo.apply(
        lambda r: unrecovered_observation_category(
            int(r["No. of default-relative changes"]),
            int(r["No. of default-relative changes"]) - int(r["No. of graph-imputed changes"]),
            int(r["No. of evidence-supported changes"]) + int(r["No. of mixed observed changes"]),
        ),
        axis=1,
    )
    pseudo_summary = summarize_recovery_split(
        pseudo,
        "Pseudo-bulk RNA-seq",
        [f"SC-{i}" for i in range(1, 9)],
    )
    if not pseudo_summary.empty:
        overall = {
            "Benchmark": "Pseudo-bulk RNA-seq",
            "Dataset": "Overall",
            "No. of Non-default Truth Alleles": int(pseudo_summary["No. of Non-default Truth Alleles"].sum()),
            "PanTCR Exact Recovered": int(pseudo_summary["PanTCR Exact Recovered"].sum()),
            "Observed-defining-site Supported": int(pseudo_summary["Observed-defining-site Supported"].sum()),
            "Observed + Graph-completed": int(pseudo_summary["Observed + Graph-completed"].sum()),
            "Graph-completed Only": int(pseudo_summary["Graph-completed Only"].sum()),
            "Not Recovered": int(pseudo_summary["Not Recovered"].sum()),
            "Not Recovered: Truth-base Observed": int(pseudo_summary["Not Recovered: Truth-base Observed"].sum()),
            "Not Recovered: Covered Without Truth-base": int(pseudo_summary["Not Recovered: Covered Without Truth-base"].sum()),
            "Not Recovered: No Retained Defining-site Coverage": int(pseudo_summary["Not Recovered: No Retained Defining-site Coverage"].sum()),
        }
        overall["PanTCR Recovery Rate"] = overall["PanTCR Exact Recovered"] / overall["No. of Non-default Truth Alleles"]
        pseudo_summary = pd.concat([pseudo_summary, pd.DataFrame([overall])], ignore_index=True)
    parts.append(pseudo_summary)
    pseudo_unrec = pseudo[pseudo["Recovered by PanTCR"] != "Yes"].copy()
    unrecovered_parts.append(
        pseudo_unrec.assign(
            Benchmark="Pseudo-bulk RNA-seq",
            Dataset=pseudo_unrec["SC ID"].astype(str),
            Sample=pseudo_unrec["Dataset ID"].astype(str),
            Gene=pseudo_unrec["Gene"].astype(str),
            **{"Truth allele": pseudo_unrec["Truth allele"].astype(str)},
            **{"PanTCR status": pseudo_unrec["PanTCR recovery status"].astype(str)},
            **{"No. of retained covered defining sites": pseudo_unrec["No. of default-relative changes"].astype(int) - pseudo_unrec["No. of graph-imputed changes"].astype(int)},
            **{"No. of retained truth-base defining sites": pseudo_unrec["No. of evidence-supported changes"].astype(int) + pseudo_unrec["No. of mixed observed changes"].astype(int)},
            **{"Unrecovered observation category": pseudo_unrec["Unrecovered observation category"]},
            **{"Evidence support stratum": pseudo_unrec["Evidence support stratum"].astype(str)},
        )[
            [
                "Benchmark",
                "Dataset",
                "Sample",
                "Gene",
                "Truth allele",
                "PanTCR status",
                "No. of default-relative changes",
                "No. of retained covered defining sites",
                "No. of retained truth-base defining sites",
                "Unrecovered observation category",
                "Evidence support stratum",
            ]
        ]
    )

    if not parts:
        return pd.DataFrame(), pd.DataFrame()
    summary = pd.concat(parts, ignore_index=True)
    unrecovered_detail = pd.concat(unrecovered_parts, ignore_index=True) if unrecovered_parts else pd.DataFrame()
    return summary, unrecovered_detail


def summarize(detail: pd.DataFrame, label: str) -> pd.DataFrame:
    if detail.empty:
        return pd.DataFrame()
    g = (
        detail.groupby(["SC ID"], dropna=False)
        .agg(
            **{
                "No. of non-reference truth alleles": ("Truth allele", "size"),
                "No. recovered by PanTCR": ("Recovered by PanTCR", lambda s: int((s == "Yes").sum())),
                "No. unrecovered by PanTCR": ("Recovered by PanTCR", lambda s: int((s != "Yes").sum())),
                "No. of default-relative changes": ("No. of default-relative changes", "sum"),
                "No. of evidence-supported changes": ("No. of evidence-supported changes", "sum"),
                "No. of mixed observed changes": ("No. of mixed observed changes", "sum"),
                "No. of conflicting observed changes": ("No. of conflicting observed changes", "sum"),
                "No. of graph-imputed changes": ("No. of graph-imputed changes", "sum"),
            }
        )
        .reset_index()
    )
    total = pd.DataFrame(
        [
            {
                "SC ID": "Overall",
                "No. of non-reference truth alleles": int(len(detail)),
                "No. recovered by PanTCR": int((detail["Recovered by PanTCR"] == "Yes").sum()),
                "No. unrecovered by PanTCR": int((detail["Recovered by PanTCR"] != "Yes").sum()),
                "No. of default-relative changes": int(detail["No. of default-relative changes"].sum()),
                "No. of evidence-supported changes": int(detail["No. of evidence-supported changes"].sum()),
                "No. of mixed observed changes": int(detail["No. of mixed observed changes"].sum()),
                "No. of conflicting observed changes": int(detail["No. of conflicting observed changes"].sum()),
                "No. of graph-imputed changes": int(detail["No. of graph-imputed changes"].sum()),
            }
        ]
    )
    out = pd.concat([g, total], ignore_index=True)
    out.insert(0, "Set", label)
    out["PanTCR recovery rate"] = out["No. recovered by PanTCR"] / out["No. of non-reference truth alleles"]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--truth-csv",
        default=str(EXPERIMENT_DIR / "generated" / "evidence_analysis" / "per_truth_call_status.csv"),
    )
    ap.add_argument(
        "--semi-truth-csv",
        default=str(
            PACKAGE_ROOT
            / "experiments"
            / "07_semi_synthetic_airr_benchmark"
            / "generated"
            / "evidence_analysis"
            / "per_truth_call_status.csv"
        ),
    )
    ap.add_argument("--ref-root", default=str(PACKAGE_ROOT / "data" / "ref"))
    ap.add_argument(
        "--out-dir",
        default=str(EXPERIMENT_DIR / "generated" / "fig4e_nonreference_variant_audit"),
    )
    ap.add_argument("--min-naive", type=float, default=2)
    args = ap.parse_args()

    detail, fig4e_detail, summary, combined_summary, unrecovered_detail, definitions = build_audit(
        Path(args.truth_csv),
        Path(args.semi_truth_csv),
        Path(args.ref_root),
        Path(args.out_dir),
        min_naive=args.min_naive,
    )
    out_dir = Path(args.out_dir)
    detail.to_csv(out_dir / "strict_default_discordant_allele_audit.csv", index=False, encoding="utf-8-sig")
    fig4e_detail.to_csv(out_dir / "fig4e_nonreference_allele_audit.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(out_dir / "fig4e_nonreference_summary.csv", index=False, encoding="utf-8-sig")
    combined_summary.to_csv(out_dir / "figure4_nondefault_recovery_split_summary.csv", index=False, encoding="utf-8-sig")
    unrecovered_detail.to_csv(out_dir / "figure4_unrecovered_nondefault_alleles.csv", index=False, encoding="utf-8-sig")
    definitions.to_csv(out_dir / "definitions.csv", index=False, encoding="utf-8-sig")
    recovered = fig4e_detail[fig4e_detail["Recovered by PanTCR"] == "Yes"].copy()
    table_cols = [
        "SC ID",
        "Gene",
        "Truth allele",
        "Compared default allele",
        "Default-relative changes",
        "Evidence-supported changes",
        "Mixed observed changes",
        "Graph-imputed changes",
        "Matched PanTCR sequence",
    ]
    table_rows = recovered[table_cols].copy()
    table_rows = table_rows.rename(
        columns={
            "SC ID": "Dataset",
            "Default-relative changes": "Change Relative to Default Allele",
            "Evidence-supported changes": "Evidence-supported Defining Site",
            "Mixed observed changes": "Mixed Observed Defining Site",
            "Graph-imputed changes": "Graph-imputed Defining Site",
            "Matched PanTCR sequence": "Inferred Sequence",
        }
    )
    for col in [
        "Evidence-supported Defining Site",
        "Mixed Observed Defining Site",
        "Graph-imputed Defining Site",
    ]:
        table_rows[col] = table_rows[col].replace("", "None")
    table_rows.to_csv(out_dir / "fig4e_recovered_nonreference_alleles_for_table.csv", index=False, encoding="utf-8-sig")
    recovered.to_csv(out_dir / "fig4e_recovered_nonreference_alleles.csv", index=False, encoding="utf-8-sig")
    with (out_dir / "workbook_data.json").open("w", encoding="utf-8") as fh:
        json.dump(
            {
                "combined_summary": combined_summary.fillna("").to_dict(orient="records"),
                "unrecovered_detail": unrecovered_detail.fillna("").to_dict(orient="records"),
                "summary": summary.fillna("").to_dict(orient="records"),
                "table_recovered": table_rows.fillna("").to_dict(orient="records"),
                "recovered": recovered.fillna("").to_dict(orient="records"),
                "fig4e_detail": fig4e_detail.fillna("").to_dict(orient="records"),
                "strict_detail": detail.fillna("").to_dict(orient="records"),
                "definitions": definitions.fillna("").to_dict(orient="records"),
            },
            fh,
            ensure_ascii=False,
        )
    with (out_dir / "manifest.json").open("w", encoding="utf-8") as fh:
        json.dump(
            {
                "truth_csv": display_path(args.truth_csv),
                "ref_root": display_path(args.ref_root),
                "min_naive": args.min_naive,
                "strict_default_discordant_rows": int(len(detail)),
                "revised_fig4e_rows": int(len(fig4e_detail)),
                "strict_recovered": int((detail["Recovered by PanTCR"] == "Yes").sum()),
                "revised_fig4e_recovered": int((fig4e_detail["Recovered by PanTCR"] == "Yes").sum()),
            },
            fh,
            ensure_ascii=False,
            indent=2,
        )
    print(f"Wrote Fig. 4E non-reference variant audit to {out_dir}")
    print(
        "strict_rows={}; revised_fig4e_rows={}; revised_recovered={}".format(
            len(detail),
            len(fig4e_detail),
            int((fig4e_detail["Recovered by PanTCR"] == "Yes").sum()),
        )
    )


if __name__ == "__main__":
    main()
