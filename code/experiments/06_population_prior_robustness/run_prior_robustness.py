#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run matched/mismatched/global/no-prior PanTCR inference on existing folds."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SIMU_ROOT = PACKAGE_ROOT / "experiments" / "01_in_silico_trbv_benchmarks"
sys.path.append(str(PACKAGE_ROOT / "code" / "experiments" / "00_benchmark_utils"))
from common_pantcr_io import sample_id_from_filename, split_words  # noqa: E402


DEFAULT_MISMATCH = {
    "AFR": "EUR",
    "EUR": "AFR",
    "EAS": "EUR",
    "AMR": "EAS",
    "SAS": "EUR",
}


def parse_map(text: str) -> dict:
    if not text:
        return dict(DEFAULT_MISMATCH)
    out = {}
    for item in text.replace(",", " ").split():
        if ":" not in item:
            continue
        a, b = item.split(":", 1)
        out[a] = b
    return out or dict(DEFAULT_MISMATCH)


def run(cmd, dry_run=False):
    print(" ".join(str(x) for x in cmd))
    if not dry_run:
        subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simu-root", default=str(DEFAULT_SIMU_ROOT))
    ap.add_argument("--results-root", default="", help="Optional directory containing results/validation and results/pang")
    ap.add_argument("--scripts-root", default="", help="Optional directory containing code/pantcr/infer_genotype_bayes.py")
    ap.add_argument("--expr", default="expr_ScenarioA")
    ap.add_argument("--genes", default="V J")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--min-naive", type=int, default=1)
    ap.add_argument("--penalty-K", dest="penalty_K", default="0", help="Synthetic-run K=2 penalty passed as --penalty_K.")
    ap.add_argument("--pi-min", dest="pi_min", default="0.1", help="Synthetic-run minimum mixture proportion passed as --pi_min.")
    ap.add_argument("--modes", default="matched mismatch global no_prior")
    ap.add_argument("--mismatch-map", default="")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = Path(args.simu_root)
    results = Path(args.results_root) if args.results_root else root / "results"
    scripts_root = Path(args.scripts_root) if args.scripts_root else PACKAGE_ROOT / "code" / "pantcr"
    infer_script = scripts_root / "infer_genotype_bayes.py"
    if not infer_script.exists():
        raise FileNotFoundError(f"Missing inference script: {infer_script}")
    mismatch_map = parse_map(args.mismatch_map)
    genes = split_words(args.genes, ["V", "J"])
    modes = split_words(args.modes, ["matched", "mismatch", "global", "no_prior"])

    for gene in genes:
        for fold in range(args.folds):
            val_dir = results / "validation" / args.expr / gene / f"fold_{fold}"
            graph_dir = results / "pang" / args.expr / gene / f"fold_{fold}"
            if not val_dir.exists():
                continue
            for sample_csv in sorted(val_dir.glob("*.csv")):
                pop, seed = sample_id_from_filename(sample_csv)
                if pop is None or seed is None:
                    continue
                for mode in modes:
                    if mode == "matched":
                        method_name = "Bayes_priorMatched"
                        pop_arg = pop
                        use_graph = True
                    elif mode == "mismatch":
                        method_name = "Bayes_priorMismatch"
                        pop_arg = mismatch_map.get(pop, "EUR")
                        use_graph = True
                    elif mode == "global":
                        method_name = "Bayes_priorGlobal"
                        pop_arg = ""
                        use_graph = True
                    elif mode == "no_prior":
                        method_name = "BayesNoPrior_no_prior"
                        pop_arg = ""
                        use_graph = False
                    else:
                        continue

                    if use_graph and not graph_dir.exists():
                        raise FileNotFoundError(
                            f"Missing graph directory for {args.expr} {gene} fold_{fold}: {graph_dir}"
                        )

                    out_dir = results / "infer" / args.expr / method_name / pop / gene
                    out_dir.mkdir(parents=True, exist_ok=True)
                    out_file = out_dir / f"infer_{pop}_{seed}.{gene}.csv"

                    cmd = [
                        args.python,
                        str(infer_script),
                        "--sample_csv",
                        str(sample_csv),
                        "--min_naive",
                        str(args.min_naive),
                        "--penalty_K",
                        str(args.penalty_K),
                        "--pi_min",
                        str(args.pi_min),
                        "--out",
                        str(out_file),
                    ]
                    if pop_arg:
                        cmd += ["--population_id", pop_arg]
                    if use_graph:
                        cmd += ["--pangenome_dir", str(graph_dir)]
                    run(cmd, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
