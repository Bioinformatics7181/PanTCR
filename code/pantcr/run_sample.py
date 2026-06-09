#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run PanTCR on one MiXCR-processed TRB clone table.

This convenience wrapper is intended for the GitHub-facing application example.
It executes the three user-facing PanTCR steps for V and J genes:

1. coordinate-projected observation extraction;
2. graph-guided Bayesian allele inference;
3. post-inference retained-evidence/imputation audit.

It does not run MiXCR itself. Users should first generate a MiXCR clone table
with the RNA-seq workflow described in the README.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd: list[str]) -> None:
    print("Running:", " ".join(str(part) for part in cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    package_root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser(description="Run PanTCR inference and evidence audit on one MiXCR TRB clone table.")
    ap.add_argument("--clones", "--clone-table", dest="clones", required=True, type=Path, help="MiXCR clones_TRB.tsv file.")
    ap.add_argument("--sample-id", required=True, help="Sample identifier used in output filenames.")
    ap.add_argument("--out-dir", required=True, type=Path, help="Output directory.")
    ap.add_argument("--graph-root", type=Path, default=package_root / "data" / "pretrained_graphs" / "trb_airr1",
                    help="Graph root containing V/ and J/ subdirectories.")
    ap.add_argument("--reference-dir", type=Path, default=package_root / "data" / "ref",
                    help="Reference directory containing IMGT_TRBV_pro.tsv, IMGT_TRBJ_pro.tsv, and IMGT_TRB_default.csv.")
    ap.add_argument("--population-id", default="", help="Optional population label for population-specific graph weights.")
    ap.add_argument("--min-naive", type=int, default=2, help="Naive-diversity threshold for retained evidence during inference/audit.")
    ap.add_argument("--penalty-K", type=float, default=3.0, help="Extra penalty for K=2 versus K=1.")
    ap.add_argument("--pi-min", type=float, default=0.1, help="Minimum minor component proportion for K=2 calls.")
    ap.add_argument("--max-candidates", type=int, default=40, help="Maximum candidate allele paths per gene.")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    scripts = Path(__file__).resolve().parent

    for gene_type in ("V", "J"):
        ref_name = "IMGT_TRBV_pro.tsv" if gene_type == "V" else "IMGT_TRBJ_pro.tsv"
        prefix = args.out_dir / f"{args.sample_id}.{gene_type}"
        observation_csv = prefix.with_name(prefix.name + "_sequences.csv")
        infer_csv = args.out_dir / f"infer_{args.sample_id}.{gene_type}.csv"
        audit_prefix = args.out_dir / f"audit_{args.sample_id}.{gene_type}"

        run_cmd([
            sys.executable,
            str(scripts / "collect_mutation.py"),
            "--input", str(args.clones),
            "--gene", gene_type,
            "--ref", str(args.reference_dir / ref_name),
            "--prefix", str(prefix),
        ])

        run_cmd([
            sys.executable,
            str(scripts / "infer_genotype_bayes.py"),
            "--sample_csv", str(observation_csv),
            "--out", str(infer_csv),
            "--pangenome_dir", str(args.graph_root / gene_type),
            "--population_id", args.population_id,
            "--min_naive", str(args.min_naive),
            "--penalty_K", str(args.penalty_K),
            "--pi_min", str(args.pi_min),
            "--max_candidates", str(args.max_candidates),
        ])

        run_cmd([
            sys.executable,
            str(scripts / "audit_inference_evidence.py"),
            "--infer-csv", str(infer_csv),
            "--sample-csv", str(observation_csv),
            "--reference-csv", str(args.reference_dir / "IMGT_TRB_default.csv"),
            "--out-prefix", str(audit_prefix),
            "--min-naive", str(args.min_naive),
        ])

    print(f"PanTCR sample run completed: {args.out_dir}")


if __name__ == "__main__":
    main()
