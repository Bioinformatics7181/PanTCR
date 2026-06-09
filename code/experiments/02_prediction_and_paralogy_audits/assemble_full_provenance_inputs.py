#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assemble the normalized provenance input tree used by experiment 02.

The script does not calculate recovery metrics. It copies the generated
per-truth and per-prediction records from upstream experiment outputs into the
directory layout expected by the S5/S6 observed-region and sequence-only
primary-table builders.
"""

from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[3]
EXPERIMENT_DIR = PACKAGE_ROOT / "experiments" / "02_prediction_and_paralogy_audits"

DEFAULT_COUNT_ROOT = EXPERIMENT_DIR / "generated" / "count_matrix_and_coverage_strata"
DEFAULT_SEMI_EVIDENCE_DIR = (
    PACKAGE_ROOT / "experiments" / "07_semi_synthetic_airr_benchmark" / "generated" / "evidence_analysis"
)
DEFAULT_SCBULK_EVIDENCE_DIR = (
    PACKAGE_ROOT / "experiments" / "08_pseudo_bulk_rnaseq_benchmark" / "generated" / "evidence_analysis"
)
DEFAULT_OUT_ROOT = EXPERIMENT_DIR / "generated" / "full_provenance_inputs"

REQUIRED_FILES = ("per_truth_call_status.csv", "per_prediction_status.csv")
MANAGED_RELATIVE_DIRS = (
    Path("01_count_matrix_and_coverage_strata"),
    Path("12_semi_simu_evidence_analysis"),
    Path("13_scbulk_real_evidence_analysis"),
)

SCENARIO_DIRS = {
    "expr_ScenarioA": "scenario_a",
    "expr_ScenarioB": "scenario_b",
    "expr_ScenarioC": "scenario_c",
    "expr_FullLength": "full_length",
}


@dataclass(frozen=True)
class CopySpec:
    label: str
    source_dir: Path
    dest_dir: Path


def require_files(source_dir: Path, label: str) -> list[Path]:
    missing = [name for name in REQUIRED_FILES if not (source_dir / name).is_file()]
    if missing:
        missing_text = ", ".join(missing)
        raise FileNotFoundError(f"{label} is missing required file(s): {missing_text} in {source_dir}")
    return [source_dir / name for name in REQUIRED_FILES]


def copy_required_files(spec: CopySpec) -> list[dict[str, str]]:
    source_files = require_files(spec.source_dir, spec.label)
    spec.dest_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for source in source_files:
        dest = spec.dest_dir / source.name
        shutil.copy2(source, dest)
        rows.append(
            {
                "label": spec.label,
                "source_file": str(source),
                "dest_file": str(dest),
                "bytes": str(dest.stat().st_size),
            }
        )
    return rows


def build_copy_specs(
    count_root: Path,
    semi_evidence_dir: Path,
    scbulk_evidence_dir: Path,
    out_root: Path,
) -> list[CopySpec]:
    specs: list[CopySpec] = []
    for expr, scenario_dir in SCENARIO_DIRS.items():
        specs.append(
            CopySpec(
                label=f"in_silico:{expr}",
                source_dir=count_root / expr,
                dest_dir=out_root / "01_count_matrix_and_coverage_strata" / scenario_dir,
            )
        )
    specs.append(
        CopySpec(
            label="semi_synthetic_airr",
            source_dir=semi_evidence_dir,
            dest_dir=out_root / "12_semi_simu_evidence_analysis" / "results",
        )
    )
    specs.append(
        CopySpec(
            label="pseudo_bulk_rnaseq",
            source_dir=scbulk_evidence_dir,
            dest_dir=out_root / "13_scbulk_real_evidence_analysis" / "results",
        )
    )
    return specs


def assemble_inputs(
    count_root: Path,
    semi_evidence_dir: Path,
    scbulk_evidence_dir: Path,
    out_root: Path,
    clean: bool = True,
) -> pd.DataFrame:
    if clean:
        for rel_dir in MANAGED_RELATIVE_DIRS:
            target = out_root / rel_dir
            if target.exists():
                shutil.rmtree(target)
        manifest_path = out_root / "input_manifest.csv"
        if manifest_path.exists():
            manifest_path.unlink()

    rows: list[dict[str, str]] = []
    for spec in build_copy_specs(count_root, semi_evidence_dir, scbulk_evidence_dir, out_root):
        rows.extend(copy_required_files(spec))
    manifest = pd.DataFrame(rows)
    out_root.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out_root / "input_manifest.csv", index=False)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count-root", default=str(DEFAULT_COUNT_ROOT))
    parser.add_argument("--semi-evidence-dir", default=str(DEFAULT_SEMI_EVIDENCE_DIR))
    parser.add_argument("--scbulk-evidence-dir", default=str(DEFAULT_SCBULK_EVIDENCE_DIR))
    parser.add_argument("--out-root", default=str(DEFAULT_OUT_ROOT))
    parser.add_argument("--no-clean", action="store_true", help="Do not remove previously assembled managed subdirectories before copying.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = assemble_inputs(
        count_root=Path(args.count_root),
        semi_evidence_dir=Path(args.semi_evidence_dir),
        scbulk_evidence_dir=Path(args.scbulk_evidence_dir),
        out_root=Path(args.out_root),
        clean=not args.no_clean,
    )
    print(f"Wrote normalized full-provenance inputs to {Path(args.out_root)}")
    print(f"Copied {len(manifest)} files")


if __name__ == "__main__":
    main()
