#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare population-agnostic leave-allele-out simulation inputs.

This script writes generated experiment files: genotypes, manifests, and target
metadata. It reads the shared PanTCR reference files without modifying them.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

CODE_EXPERIMENT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = CODE_EXPERIMENT_DIR.parents[2]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / CODE_EXPERIMENT_DIR.name
sys.path.insert(0, str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))

from common_pantcr_io import gene_type_from_gene, load_default_refs, load_pmtr_sequences, load_trb_index, trim_sequence


DEFAULT_TARGETS = [
    "TRBV7-7*03",
    "TRBV6-6*02",
    "TRBV12-2*02",
    "TRBV11-2*02",
    "TRBV24-1*02",
    "TRBV30*02",
    "TRBV5-1*02",
    "TRBV5-4*05",
    "TRBV5-5*02",
    "TRBV6-1*03",
    "TRBV6-6*08",
    "TRBV7-9*09",
    "TRBV7-9*10",
    "TRBV12-2*03",
    "TRBV12-3*05",
    "TRBV12-5*02",
]

POP_COLUMNS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def parse_gene_name(allele: str) -> str:
    return str(allele).split("*", 1)[0].strip()


def gene_type_or_empty(gene: str) -> str:
    try:
        return gene_type_from_gene(gene)
    except ValueError:
        return ""


def parse_target_alleles(target_alleles: str | None, target_file: Path | None = None) -> list[str]:
    tokens: list[str] = []
    if target_alleles:
        tokens.extend(x for x in re.split(r"[\s,]+", target_alleles.strip()) if x)
    if target_file is not None:
        with target_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue
                tokens.extend(x for x in re.split(r"[\s,]+", line) if x)
    if not tokens:
        tokens = list(DEFAULT_TARGETS)

    seen: set[str] = set()
    out: list[str] = []
    for token in tokens:
        if token not in seen:
            out.append(token)
            seen.add(token)
    return out


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().isin({"TRUE", "T", "1", "YES", "Y"})


def add_gene_and_weights(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "gene" not in out.columns:
        out["gene"] = out["allele"].map(parse_gene_name)
    if "pooled_count" not in out.columns:
        counts = []
        for col in POP_COLUMNS:
            if col in out.columns:
                counts.append(pd.to_numeric(out[col], errors="coerce").fillna(0.0))
        if counts:
            out["pooled_count"] = sum(counts)
        else:
            out["pooled_count"] = 1.0
    out["sequence"] = out["sequence"].astype(str).str.upper()
    return out


def allele_eval_sequence(
    row: pd.Series,
    pmtr_sequences: dict[str, str],
    index_df: pd.DataFrame,
    gene_type: str = "V",
) -> str:
    allele = str(row["allele"]).strip()
    gene = str(row["gene"]).strip()
    seq = pmtr_sequences.get(allele)
    if seq:
        return trim_sequence(seq, gene, gene_type, index_df)
    return trim_sequence(str(row["sequence"]), gene, gene_type, index_df)


def sequence_position_class(positions: Iterable[int]) -> str:
    pos = sorted(int(x) for x in positions)
    if not pos:
        return "none"
    labels: list[str] = []
    if any(p < 75 for p in pos):
        labels.append("early")
    if any(75 <= p < 180 for p in pos):
        labels.append("middle")
    if any(p >= 180 for p in pos):
        labels.append("late")
    return "/".join(labels)


def defining_positions(seq: str, ref: str) -> list[int]:
    n = min(len(seq), len(ref))
    return [i for i in range(n) if seq[i] != ref[i]]


def weighted_choice_indices(group: pd.DataFrame, rng: np.random.Generator, mode: str) -> np.ndarray:
    n = len(group)
    if n == 0:
        raise ValueError("Cannot sample from an empty allele group.")
    if mode == "uniform":
        probs = np.ones(n, dtype=float) / n
    else:
        weights = pd.to_numeric(group["pooled_count"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        total = float(weights.sum())
        probs = weights / total if total > 0 else np.ones(n, dtype=float) / n
    return rng.choice(np.arange(n), size=2, replace=True, p=probs)


def choose_reference_background_row(group: pd.DataFrame) -> pd.Series:
    for col in ["is_default_noCDR3", "is_default"]:
        if col in group.columns:
            default_rows = group[bool_series(group[col])]
            if not default_rows.empty:
                return default_rows.sort_values(["pooled_count", "allele"], ascending=[False, True]).iloc[0]
    return group.sort_values(["pooled_count", "allele"], ascending=[False, True]).iloc[0]


def sample_background_genotype(
    allele_pool: pd.DataFrame,
    rng: np.random.Generator,
    sampling_mode: str,
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for gene, group in allele_pool.groupby("gene", sort=True):
        group = group.reset_index(drop=True)
        idx_a, idx_b = weighted_choice_indices(group, rng, sampling_mode)
        row_a = group.iloc[int(idx_a)]
        row_b = group.iloc[int(idx_b)]
        rows.append(
            {
                "gene": gene,
                "allele_A": row_a["allele"],
                "seq_A": row_a["sequence"],
                "allele_B": row_b["allele"],
                "seq_B": row_b["sequence"],
            }
        )
    return pd.DataFrame(rows)


def force_target_allele(
    background: pd.DataFrame,
    target_row: pd.Series,
    background_row: pd.Series,
) -> pd.DataFrame:
    out = background.copy()
    gene = str(target_row["gene"])
    mask = out["gene"].astype(str) == gene
    if not mask.any():
        raise ValueError(f"Target gene is missing from background genotype: {gene}")
    out.loc[mask, "allele_A"] = str(background_row["allele"])
    out.loc[mask, "seq_A"] = str(background_row["sequence"])
    out.loc[mask, "allele_B"] = str(target_row["allele"])
    out.loc[mask, "seq_B"] = str(target_row["sequence"])
    return out


def resolve_path(path_value: str | Path, base: Path | None = None) -> Path:
    path = Path(path_value)
    if not path.is_absolute() and base is not None:
        path = base / path
    return path.resolve()


def write_manifest(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def build_target_metadata(
    clean_df: pd.DataFrame,
    targets: list[str],
    pmtr_sequences: dict[str, str],
    default_refs: dict[str, str],
    index_df: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, pd.Series]]:
    rows: list[dict[str, object]] = []
    target_rows: dict[str, pd.Series] = {}
    for target in targets:
        hits = clean_df[clean_df["allele"].astype(str) == target]
        if hits.empty:
            raise ValueError(f"Target allele not found in pmTR table: {target}")
        row = hits.iloc[0].copy()
        gene = str(row["gene"])
        gene_type = gene_type_from_gene(gene)
        target_eval_seq = allele_eval_sequence(row, pmtr_sequences, index_df, gene_type)
        default_full = default_refs.get(gene, "")
        default_eval_seq = trim_sequence(default_full, gene, gene_type, index_df) if default_full else ""
        def_pos = defining_positions(target_eval_seq, default_eval_seq)
        pop_counts = {
            pop: float(pd.to_numeric(pd.Series([row.get(pop, 0)]), errors="coerce").fillna(0.0).iloc[0])
            for pop in POP_COLUMNS
        }
        max_pop = max(pop_counts, key=pop_counts.get)
        target_rows[target] = row
        rows.append(
            {
                "target_allele": target,
                "target_gene": gene,
                "target_gene_type": gene_type,
                "target_eval_seq": target_eval_seq,
                "target_eval_len": len(target_eval_seq),
                "default_eval_seq": default_eval_seq,
                "default_eval_len": len(default_eval_seq),
                "n_defining_variants_in_evaluation_region": len(def_pos),
                "defining_positions": ";".join(str(x) for x in def_pos),
                "first_defining_position": min(def_pos) if def_pos else "",
                "last_defining_position": max(def_pos) if def_pos else "",
                "variant_position_class": sequence_position_class(def_pos),
                "evaluable_in_mature_region": bool(len(def_pos) >= 1),
                "pmtr_level": row.get("level", ""),
                "pmtr_is_new": row.get("is_new", ""),
                "pmtr_is_default": row.get("is_default", ""),
                "pmtr_is_new_noCDR3": row.get("is_new_noCDR3", ""),
                "pmtr_is_default_noCDR3": row.get("is_default_noCDR3", ""),
                "pmtr_mutations": row.get("mutations", ""),
                "max_population_count": pop_counts[max_pop],
                "max_population": max_pop,
            }
        )
    return pd.DataFrame(rows), target_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--simu-root", default="simu")
    parser.add_argument("--run-root", default="runs")
    parser.add_argument("--expr", default="expr_leaveout_allele_specific")
    parser.add_argument("--target-alleles", default=" ".join(DEFAULT_TARGETS))
    parser.add_argument("--target-file", default="")
    parser.add_argument("--n-train", type=int, default=40)
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--seed-start-train", type=int, default=0)
    parser.add_argument("--seed-start-test", type=int, default=10000)
    parser.add_argument("--train-pop", default="TRAIN")
    parser.add_argument("--test-pop", default="PANEL")
    parser.add_argument("--sampling-mode", choices=["pooled", "uniform"], default="pooled")
    parser.add_argument("--random-seed", type=int, default=20260521)
    parser.add_argument("--allow-unevaluable", action="store_true")
    parser.add_argument("--pmtr", default="ref/pmTR_TRB_V_J_cleaned.csv")
    parser.add_argument("--pmtr-ref", default="ref/pmTR_TRB_V_J_cleaned.csv")
    parser.add_argument("--default-ref", default="ref/IMGT_TRB_default.csv")
    parser.add_argument("--index", default="ref/TRB_index.csv")
    args = parser.parse_args()

    simu_root = resolve_path(args.simu_root, EXPERIMENT_DIR)
    run_root = resolve_path(args.run_root, EXPERIMENT_DIR)
    pmtr_path = resolve_path(args.pmtr, simu_root)
    pmtr_ref_path = resolve_path(args.pmtr_ref, simu_root)
    default_ref_path = resolve_path(args.default_ref, simu_root)
    index_path = resolve_path(args.index, simu_root)
    target_file = resolve_path(args.target_file, EXPERIMENT_DIR) if args.target_file else None

    for path in [pmtr_path, pmtr_ref_path, default_ref_path, index_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    targets = parse_target_alleles(args.target_alleles, target_file)
    rng = np.random.default_rng(args.random_seed)

    clean_df = add_gene_and_weights(pd.read_csv(pmtr_path))
    index_df = load_trb_index(index_path)
    pmtr_sequences = load_pmtr_sequences(pmtr_ref_path)
    default_refs = load_default_refs(default_ref_path)

    target_meta, target_rows = build_target_metadata(clean_df, targets, pmtr_sequences, default_refs, index_df)
    if not args.allow_unevaluable:
        bad = target_meta[~target_meta["evaluable_in_mature_region"]]
        if not bad.empty:
            raise ValueError(
                "Some targets have no defining variants in the evaluation region: "
                + ", ".join(bad["target_allele"].astype(str).tolist())
            )

    target_eval_by_gene = {
        (str(r["target_gene"]), str(r["target_eval_seq"])) for _, r in target_meta.iterrows()
    }
    clean_df["gene_type"] = clean_df["gene"].map(gene_type_or_empty)
    clean_df["eval_seq"] = clean_df.apply(
        lambda r: allele_eval_sequence(r, pmtr_sequences, index_df, str(r["gene_type"]))
        if str(r["gene_type"]) in {"V", "J"}
        else str(r["sequence"]).upper(),
        axis=1,
    )
    clean_df["is_target_name"] = clean_df["allele"].astype(str).isin(targets)
    clean_df["is_target_eval_sequence"] = clean_df.apply(
        lambda r: (str(r["gene"]), str(r["eval_seq"])) in target_eval_by_gene,
        axis=1,
    )
    allele_pool = clean_df[~clean_df["is_target_name"] & ~clean_df["is_target_eval_sequence"]].copy()

    missing_genes = sorted(set(clean_df["gene"]) - set(allele_pool["gene"]))
    if missing_genes:
        raise ValueError("Target exclusion removed all alleles for genes: " + ", ".join(missing_genes))

    output_root = run_root
    manifest_dir = output_root / "results" / "leaveout" / args.expr
    samples_root = output_root / "samples" / args.expr
    manifest_dir.mkdir(parents=True, exist_ok=True)

    target_meta_path = manifest_dir / "target_metadata.csv"
    target_meta.to_csv(target_meta_path, index=False)

    excluded_rows = clean_df[clean_df["is_target_name"] | clean_df["is_target_eval_sequence"]][
        ["allele", "gene", "gene_type", "eval_seq", "is_target_name", "is_target_eval_sequence"]
    ].copy()
    excluded_rows.to_csv(manifest_dir / "excluded_target_sequences.csv", index=False)

    simulation_rows: list[dict[str, object]] = []
    training_rows: list[dict[str, object]] = []
    target_rows_out: list[dict[str, object]] = []

    for i in range(args.n_train):
        seed = args.seed_start_train + i
        pop = args.train_pop
        sample_dir = samples_root / pop / f"seed{seed}"
        sample_dir.mkdir(parents=True, exist_ok=True)
        genotype_csv = sample_dir / f"genotype_{pop}_seed{seed}.csv"
        genotype = sample_background_genotype(allele_pool, rng, args.sampling_mode)
        genotype.to_csv(genotype_csv, index=False)
        row = {
            "role": "TRAIN",
            "target_allele": "NA",
            "target_gene": "NA",
            "seed": seed,
            "replicate": "NA",
            "pop": pop,
            "genotype_csv": str(genotype_csv.resolve()),
            "sample_dir": str(sample_dir.resolve()),
            "sample_prefix": f"{pop}_seed{seed}",
        }
        simulation_rows.append(row)
        training_rows.append(row)

    background_by_gene = {
        gene: choose_reference_background_row(group.reset_index(drop=True))
        for gene, group in allele_pool.groupby("gene", sort=True)
    }

    for target_index, target in enumerate(targets):
        target_row = target_rows[target]
        target_gene = str(target_row["gene"])
        background_row = background_by_gene[target_gene]
        for rep in range(args.replicates):
            seed = args.seed_start_test + target_index * args.replicates + rep
            pop = args.test_pop
            sample_dir = samples_root / pop / f"seed{seed}"
            sample_dir.mkdir(parents=True, exist_ok=True)
            genotype_csv = sample_dir / f"genotype_{pop}_seed{seed}.csv"
            background = sample_background_genotype(allele_pool, rng, args.sampling_mode)
            genotype = force_target_allele(background, target_row, background_row)
            genotype.to_csv(genotype_csv, index=False)
            row = {
                "role": "TEST",
                "target_allele": target,
                "target_gene": target_gene,
                "seed": seed,
                "replicate": rep,
                "pop": pop,
                "genotype_csv": str(genotype_csv.resolve()),
                "sample_dir": str(sample_dir.resolve()),
                "sample_prefix": f"{pop}_seed{seed}",
            }
            simulation_rows.append(row)
            target_rows_out.append(row)

    write_manifest(manifest_dir / "simulation_manifest.tsv", simulation_rows)
    write_manifest(manifest_dir / "training_manifest.tsv", training_rows)
    write_manifest(manifest_dir / "target_panel_manifest.tsv", target_rows_out)

    config = {
        "expr": args.expr,
        "run_root": str(output_root.resolve()),
        "simu_root": str(simu_root),
        "n_train": args.n_train,
        "replicates": args.replicates,
        "sampling_mode": args.sampling_mode,
        "random_seed": args.random_seed,
        "train_pop": args.train_pop,
        "test_pop": args.test_pop,
        "target_alleles": " ".join(targets),
    }
    pd.DataFrame([config]).to_csv(manifest_dir / "experiment_config.csv", index=False)

    print(f"Wrote target metadata: {target_meta_path}")
    print(f"Wrote simulation manifest: {manifest_dir / 'simulation_manifest.tsv'}")
    print(f"Training samples: {len(training_rows)}")
    print(f"Forced target test samples: {len(target_rows_out)}")


if __name__ == "__main__":
    main()
