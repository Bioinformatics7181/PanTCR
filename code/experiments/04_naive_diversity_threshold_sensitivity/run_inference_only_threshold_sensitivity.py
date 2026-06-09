#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rerun inference-only min_naive sensitivity on fixed existing fold graphs.

This differs from a full graph-rebuild sensitivity analysis. Here the fold
graph is kept fixed and only the inference `--min_naive` changes.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

PACKAGE_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SIMU_ROOT = PACKAGE_ROOT / "experiments" / "01_in_silico_trbv_benchmarks"
sys.path.append(str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))
from common_pantcr_io import sample_id_from_filename, split_words  # noqa: E402


def run(cmd: list[str], dry_run: bool = False) -> None:
    print(" ".join(str(x) for x in cmd))
    if not dry_run:
        subprocess.run(cmd, check=True)


def append_manifest(manifest_path: Path, row: dict) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(manifest_path, mode="a", header=not manifest_path.exists(), index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--simu-root", default=str(DEFAULT_SIMU_ROOT))
    ap.add_argument("--results-root", default="", help="Directory containing staged results/validation and results/pang")
    ap.add_argument("--scripts-root", default="", help="Directory containing code/pantcr/infer_genotype_bayes.py")
    ap.add_argument("--expr", default="expr_ScenarioA")
    ap.add_argument("--genes", default="V J")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--thresholds", default="0 1 2 3")
    ap.add_argument("--penalty-K", dest="penalty_K", default="0", help="Synthetic-run K=2 penalty passed as --penalty_K.")
    ap.add_argument("--pi-min", dest="pi_min", default="0.1", help="Synthetic-run minimum mixture proportion passed as --pi_min.")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = Path(args.simu_root)
    results = Path(args.results_root) if args.results_root else root / "results"
    scripts_root = Path(args.scripts_root) if args.scripts_root else PACKAGE_ROOT / "code" / "pantcr"
    infer_script = scripts_root / "infer_genotype_bayes.py"
    if not infer_script.exists():
        raise FileNotFoundError(f"Missing inference script: {infer_script}")

    validation = results / "validation" / args.expr
    fixed_graph_root = results / "pang" / args.expr
    if not validation.exists():
        raise FileNotFoundError(f"Missing validation input: {validation}")
    if not fixed_graph_root.exists():
        raise FileNotFoundError(f"Missing fixed graph input: {fixed_graph_root}")

    genes = split_words(args.genes, ["V", "J"])
    thresholds = [int(x) for x in split_words(args.thresholds, ["0", "1", "2", "3"])]
    manifest_path = results / "supplemental_inference_only_threshold" / args.expr / "inference_only_threshold_manifest.csv"
    if manifest_path.exists() and not args.dry_run:
        manifest_path.unlink()

    for threshold in thresholds:
        for gene in genes:
            for fold in range(args.folds):
                val_dir = validation / gene / f"fold_{fold}"
                graph_dir = fixed_graph_root / gene / f"fold_{fold}"
                if not val_dir.exists():
                    raise FileNotFoundError(f"Missing validation fold: {val_dir}")
                if not graph_dir.exists():
                    raise FileNotFoundError(f"Missing fixed graph fold: {graph_dir}")
                samples = sorted(val_dir.glob("*.csv"))
                append_manifest(
                    manifest_path,
                    {
                        "expr": args.expr,
                        "gene_type": gene,
                        "fold": fold,
                        "inference_min_naive": threshold,
                        "fixed_graph_dir": str(graph_dir),
                        "n_validation_samples": len(samples),
                    },
                )
                if not samples:
                    raise RuntimeError(f"No validation samples found: {val_dir}")

                for sample_csv in samples:
                    pop, seed = sample_id_from_filename(sample_csv)
                    if pop is None or seed is None:
                        continue
                    graph_method = f"Bayes_fixedGraph_inferMinNaive{threshold}"
                    graph_out_dir = results / "infer" / args.expr / graph_method / pop / gene
                    graph_out_dir.mkdir(parents=True, exist_ok=True)
                    graph_out = graph_out_dir / f"infer_{pop}_{seed}.{gene}.csv"
                    if not (args.resume and graph_out.exists()):
                        run(
                            [
                                args.python,
                                str(infer_script),
                                "--sample_csv",
                                str(sample_csv),
                                "--population_id",
                                pop,
                                "--min_naive",
                                str(threshold),
                                "--penalty_K",
                                str(args.penalty_K),
                                "--pi_min",
                                str(args.pi_min),
                                "--pangenome_dir",
                                str(graph_dir),
                                "--out",
                                str(graph_out),
                            ],
                            dry_run=args.dry_run,
                        )

                    noprior_method = f"BayesNoPrior_inferMinNaive{threshold}"
                    noprior_out_dir = results / "infer" / args.expr / noprior_method / pop / gene
                    noprior_out_dir.mkdir(parents=True, exist_ok=True)
                    noprior_out = noprior_out_dir / f"infer_{pop}_{seed}.{gene}.csv"
                    if not (args.resume and noprior_out.exists()):
                        run(
                            [
                                args.python,
                                str(infer_script),
                                "--sample_csv",
                                str(sample_csv),
                                "--population_id",
                                pop,
                                "--min_naive",
                                str(threshold),
                                "--penalty_K",
                                str(args.penalty_K),
                                "--pi_min",
                                str(args.pi_min),
                                "--out",
                                str(noprior_out),
                            ],
                            dry_run=args.dry_run,
                        )


if __name__ == "__main__":
    main()
