#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild the final leave-target-allele recovery source CSVs.

This package-facing script matches the reported S12/S13 definition:
exact recovery, observed-region compatible partial recovery, and not recovered
or incompatible. By default it reads the generated summary from the
leave-allele-out experiment run directory and writes source CSVs to the
experiment `generated/summary` folder unless `--out-dir` is supplied.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import csv

import pandas as pd


CODE_EXPERIMENT_DIR = Path(__file__).resolve().parent
PACKAGE_ROOT = CODE_EXPERIMENT_DIR.parents[2]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / CODE_EXPERIMENT_DIR.name
DEFAULT_SUMMARY_DIR = (
    EXPERIMENT_DIR
    / "runs"
    / "results"
    / "leaveout"
    / "expr_leaveout_allele_specific"
    / "summary"
)
DEFAULT_METADATA_CSV = (
    EXPERIMENT_DIR
    / "runs"
    / "results"
    / "leaveout"
    / "expr_leaveout_allele_specific"
    / "target_metadata.csv"
)
DEFAULT_OUT_DIR = EXPERIMENT_DIR / "generated" / "summary"

METHOD_ORDER = ["MiXCR", "FindAlleles", "BayesNoPrior", "PanTCRLeaveout"]
METHOD_LABEL = {
    "MiXCR": "MiXCR-default",
    "FindAlleles": "FindAlleles",
    "BayesNoPrior": "PanTCR-NP",
    "PanTCRLeaveout": "PanTCR",
}


def read_final_per_sample(summary_dir: Path) -> pd.DataFrame:
    path = summary_dir / "per_target_method_status.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing final leave-target per-sample status: {path}")
    df = pd.read_csv(path)
    required = {"method", "target_allele", "target_gene", "status"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    return df


def read_metadata(metadata_csv: Path) -> pd.DataFrame:
    if not metadata_csv.exists():
        raise FileNotFoundError(f"Missing target metadata: {metadata_csv}")
    meta = pd.read_csv(metadata_csv)
    required = {"target_allele", "target_gene", "pmtr_mutations"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"{metadata_csv} missing required columns: {sorted(missing)}")
    return meta


def default_relative_changes(target_seq: object, default_seq: object) -> str:
    target = str(target_seq).strip().upper()
    default = str(default_seq).strip().upper()
    changes = []
    for pos, (target_base, default_base) in enumerate(zip(target, default)):
        if target_base != default_base:
            changes.append(f"S{default_base}{pos}{target_base}")
    if len(target) != len(default):
        changes.append(f"LEN{len(default)}>{len(target)}")
    return ";".join(changes) if changes else "Default"


def count_status(frame: pd.DataFrame, status: str) -> int:
    return int(frame["status"].astype(str).eq(status).sum())


def not_recovered_count(frame: pd.DataFrame) -> int:
    known = frame["status"].astype(str).isin(["exact_recovery", "partial_recovery"])
    return int((~known).sum())


def rate(value: int, total: int) -> str:
    return f"{value / total:.3f}" if total else "0.000"


def build_overall(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method in METHOD_ORDER:
        sub = per[per["method"].astype(str).eq(method)]
        total = int(len(sub))
        exact = count_status(sub, "exact_recovery")
        partial = count_status(sub, "partial_recovery")
        not_recovered = not_recovered_count(sub)
        rows.append(
            {
                "Method": METHOD_LABEL[method],
                "NO. of Target Test Cases": total,
                "Exact Recovery": exact,
                "Observed-region Compatible Partial Recovery": partial,
                "Not Recovered or Incompatible": not_recovered,
                "Exact Rate": rate(exact, total),
                "Exact or Partial Rate": rate(exact + partial, total),
            }
        )
    return pd.DataFrame(rows)


def build_by_target(per: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    if {"target_eval_seq", "default_eval_seq"}.issubset(meta.columns):
        meta["default_relative_change"] = meta.apply(
            lambda row: default_relative_changes(row["target_eval_seq"], row["default_eval_seq"]),
            axis=1,
        )
    else:
        meta["default_relative_change"] = meta["pmtr_mutations"].astype(str)
    meta = meta.set_index("target_allele", drop=False)
    targets = sorted(per["target_allele"].dropna().astype(str).unique())
    rows = []
    for target in targets:
        if target not in set(per["target_allele"].astype(str)):
            continue
        for method in METHOD_ORDER:
            sub = per[per["target_allele"].astype(str).eq(target) & per["method"].astype(str).eq(method)]
            total = int(len(sub))
            exact = count_status(sub, "exact_recovery")
            partial = count_status(sub, "partial_recovery")
            not_recovered = not_recovered_count(sub)
            mrow = meta.loc[target]
            rows.append(
                {
                    "Target allele": target,
                    "Target gene": str(mrow["target_gene"]),
                    "Default-relative Change": str(mrow["default_relative_change"]),
                    "NO. of Target Test Cases": total,
                    "Method": METHOD_LABEL[method],
                    "Exact Recovery": exact,
                    "Observed-region Compatible Partial Recovery": partial,
                    "Not Recovered or Incompatible": not_recovered,
                    "Exact Rate": rate(exact, total),
                    "Exact or Partial Rate": rate(exact + partial, total),
                }
            )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary-dir", type=Path, default=DEFAULT_SUMMARY_DIR)
    ap.add_argument("--metadata-csv", type=Path, default=DEFAULT_METADATA_CSV)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    per = read_final_per_sample(args.summary_dir)
    meta = read_metadata(args.metadata_csv)
    overall = build_overall(per)
    by_target = build_by_target(per, meta)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    overall_path = args.out_dir / "leaveout_allele_recovery_overall.csv"
    by_target_path = args.out_dir / "leaveout_allele_recovery_by_target.csv"
    overall.to_csv(overall_path, index=False, encoding="utf-8")
    by_target.to_csv(by_target_path, index=False, encoding="utf-8")
    print(overall_path)
    print(by_target_path)


if __name__ == "__main__":
    main()
