#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared helpers for PanTCR supplementary analyses."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd


def split_words(value: str | Sequence[str] | None, default: Sequence[str]) -> List[str]:
    if value is None:
        return list(default)
    if isinstance(value, (list, tuple)):
        return [str(x) for x in value]
    value = str(value).replace(",", " ").strip()
    return value.split() if value else list(default)


def sample_id_from_filename(path: Path) -> Tuple[Optional[str], Optional[str]]:
    """Return (population, seed string like seed0) from common PanTCR filenames."""
    m = re.search(r"([A-Za-z]+)_(seed\d+)", path.name)
    if not m:
        return None, None
    return m.group(1), m.group(2)


def gene_type_from_gene(gene: object) -> str:
    text = str(gene)
    if "TRBV" in text:
        return "V"
    if "TRBJ" in text:
        return "J"
    raise ValueError(f"Unsupported TRB gene name: {gene!r}")


def seed_number(seed_string: str) -> str:
    m = re.search(r"(\d+)", seed_string)
    return m.group(1) if m else seed_string


def read_table_auto(path: Path, nrows=None) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False)


def require_columns(df: pd.DataFrame, required: Set[str], path: Path, context: str) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{context} missing required columns {missing}: {path}")


def load_trb_index(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for col in ["FR1Begin", "CDR3Begin", "FR4Begin"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_default_refs(path: Path) -> Dict[str, str]:
    df = pd.read_csv(path)
    out = {}
    for _, row in df.iterrows():
        gene = str(row.get("Gene", "")).strip()
        seq = str(row.get("Sequence", "")).strip().upper()
        if gene and seq:
            out[gene] = seq
    return out


def load_pmtr_sequences(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.exists():
        return {}
    df = pd.read_csv(path)
    if "allele" not in df.columns:
        raise ValueError(f"pmTR reference file missing required column 'allele': {path}")
    if "sequence_trimmed" in df.columns:
        seq_col = "sequence_trimmed"
    elif "sequence" in df.columns:
        seq_col = "sequence"
    else:
        raise ValueError(f"pmTR reference file missing required sequence column 'sequence_trimmed' or 'sequence': {path}")
    allele_col = "allele"
    return {
        str(r[allele_col]).strip(): str(r[seq_col]).strip().upper()
        for _, r in df.iterrows()
        if pd.notna(r.get(allele_col)) and pd.notna(r.get(seq_col))
    }


def trim_sequence(seq: str, gene: str, gene_type: str, index_df: pd.DataFrame) -> str:
    seq = str(seq).strip().upper()
    if not seq or seq == "NAN":
        return ""
    row = index_df[index_df["baseGene"].astype(str) == str(gene)]
    if row.empty:
        return seq
    r = row.iloc[0]
    if gene_type == "V":
        if pd.notna(r.get("FR1Begin")) and pd.notna(r.get("CDR3Begin")):
            cut = int(r["CDR3Begin"]) - int(r["FR1Begin"])
            return seq[:cut] if cut > 0 else seq
        return seq
    if gene_type == "J":
        start = int(r["FR4Begin"]) if pd.notna(r.get("FR4Begin")) else 0
        return seq[start:] if len(seq) > start else seq
    return seq


def truth_records(
    label_csv: Path,
    gene_type: str,
    index_df: pd.DataFrame,
    pmtr_map: Dict[str, str],
    prefer_label_sequence: bool = False,
) -> List[dict]:
    df = pd.read_csv(label_csv)
    require_columns(df, {"gene", "allele_A", "allele_B", "seq_A", "seq_B"}, label_csv, "Genotype label CSV")
    filt = f"TRB{gene_type}"
    df = df[df["gene"].astype(str).str.contains(filt, na=False)].copy()
    rows = []
    for _, r in df.iterrows():
        gene = str(r["gene"])
        for slot in ["A", "B"]:
            allele = str(r.get(f"allele_{slot}", "")).strip()
            raw_seq = str(r.get(f"seq_{slot}", "")).strip()
            # Synthetic labels may store allele names with short or incomplete
            # sequence fields, so the default uses the catalog sequence when it
            # is available. Biological proxy-truth benchmarks should instead
            # keep the sample-level label sequence as truth because a same-name
            # catalog allele can differ from the observed/proxy truth sequence.
            if prefer_label_sequence:
                full_seq = str(raw_seq or pmtr_map.get(allele, "")).strip().upper()
            else:
                full_seq = str(pmtr_map.get(allele) or raw_seq).strip().upper()
            seq = trim_sequence(full_seq, gene, gene_type, index_df)
            if not seq and allele and allele.upper() not in {"NAN", "NONE", "NA"}:
                raise ValueError(f"Unresolved nonempty allele {allele!r} without usable sequence in {label_csv}")
            if seq:
                rows.append({
                    "gene": gene,
                    "slot": slot,
                    "allele": allele,
                    "seq": seq,
                    "seq_full_catalog": full_seq,
                    "truth_sequence_scope": "eval_trimmed",
                })
    # Deduplicate homozygous truth alleles by gene+seq.
    seen = set()
    unique = []
    for row in rows:
        key = (row["gene"], row["seq"])
        if key not in seen:
            unique.append(row)
            seen.add(key)
    return unique


def prediction_records(infer_csv: Path, gene_type: str, index_df: pd.DataFrame) -> List[dict]:
    if not infer_csv.exists():
        return []
    df = pd.read_csv(infer_csv)
    require_columns(df, {"gene"}, infer_csv, "Inference CSV")
    seq_cols = [col for col in ["seq_A", "seq_B"] if col in df.columns]
    if not seq_cols:
        raise ValueError(f"Inference CSV has no sequence columns seq_A/seq_B: {infer_csv}")
    filt = f"TRB{gene_type}"
    df = df[df["gene"].astype(str).str.contains(filt, na=False)].copy()
    rows = []
    for _, r in df.iterrows():
        gene = str(r["gene"])
        for slot in ["A", "B"]:
            seq_col = f"seq_{slot}"
            allele_col = f"allele_{slot}"
            if seq_col not in df.columns or pd.isna(r.get(seq_col)):
                continue
            seq = trim_sequence(str(r.get(seq_col, "")), gene, gene_type, index_df)
            if seq:
                rows.append({
                    "gene": gene,
                    "slot": slot,
                    "allele": str(r.get(allele_col, "")).strip(),
                    "seq": seq,
                })
    # Deduplicate homozygous predictions by gene+seq.
    seen = set()
    unique = []
    for row in rows:
        key = (row["gene"], row["seq"])
        if key not in seen:
            unique.append(row)
            seen.add(key)
    return unique


def compatible_match(x: str, y: str, gene_type: str) -> bool:
    return compatible_match_score(x, y, gene_type) > 0


def compatible_match_score(x: str, y: str, gene_type: str) -> int:
    x = str(x).strip().upper()
    y = str(y).strip().upper()
    n = min(len(x), len(y))
    if n <= 0:
        return 0
    if gene_type == "J":
        return n if x[-n:] == y[-n:] else 0
    return n if x[:n] == y[:n] else 0


def parse_ranges_part(part: str) -> List[Tuple[int, int]]:
    if part is None or pd.isna(part):
        return []
    out = []
    for seg in str(part).split(";"):
        seg = seg.strip()
        if not seg:
            continue
        toks = seg.split(",")
        if len(toks) != 2:
            raise ValueError(f"Malformed ObservedRanges segment {seg!r}; expected 'start,end'")
        try:
            a, b = int(toks[0]), int(toks[1])
        except ValueError:
            raise ValueError(f"Malformed ObservedRanges segment {seg!r}; start/end must be integers") from None
        if b < a:
            a, b = b, a
        out.append((a, b))
    return out


def covered_positions_from_observed(observed_ranges: Iterable[str]) -> Set[int]:
    covered: Set[int] = set()
    for obs in observed_ranges:
        for piece in str(obs).split("|"):
            for a, b in parse_ranges_part(piece):
                covered.update(range(a, b))
    return covered


def filter_evidence_rows(df: pd.DataFrame, min_naive: Optional[float] = None) -> pd.DataFrame:
    """Apply the same NaiveDiversityIndex row filter used by PanTCR inference."""
    if min_naive is None:
        return df
    if "NaiveDiversityIndex" not in df.columns:
        raise ValueError("Mutation evidence CSV missing required column 'NaiveDiversityIndex' for min_naive filtering")
    out = df.copy()
    out["NaiveDiversityIndex"] = pd.to_numeric(out["NaiveDiversityIndex"], errors="coerce").fillna(0)
    return out[out["NaiveDiversityIndex"] >= float(min_naive)].copy()


def evidence_coverage_by_gene(evidence_csv: Optional[Path], min_naive: Optional[float] = None) -> Dict[str, Set[int]]:
    if evidence_csv is None or not evidence_csv.exists():
        return {}
    df = pd.read_csv(evidence_csv)
    require_columns(df, {"Family", "ObservedRanges"}, evidence_csv, "Mutation evidence CSV")
    df = filter_evidence_rows(df, min_naive=min_naive)
    out: Dict[str, Set[int]] = {}
    for gene, sub in df.groupby("Family"):
        out[str(gene)] = covered_positions_from_observed(sub["ObservedRanges"].dropna().astype(str))
    return out


def evidence_observed_bases_by_gene(evidence_csv: Optional[Path], min_naive: Optional[float] = None) -> Dict[str, Dict[int, Set[str]]]:
    """Map gene -> position -> observed bases from mutation-detection rows."""
    if evidence_csv is None or not evidence_csv.exists():
        return {}
    df = pd.read_csv(evidence_csv)
    require_columns(df, {"Family", "ObservedRanges", "Sequence"}, evidence_csv, "Mutation evidence CSV")
    df = filter_evidence_rows(df, min_naive=min_naive)
    out: Dict[str, Dict[int, Set[str]]] = {}
    for row in df.itertuples(index=False):
        gene = str(getattr(row, "Family", ""))
        seq = str(getattr(row, "Sequence", "")).strip().upper()
        obs = getattr(row, "ObservedRanges", "")
        if not gene or not seq:
            continue
        gene_map = out.setdefault(gene, {})
        for piece in str(obs).split("|"):
            for a, b in parse_ranges_part(piece):
                for pos in range(a, b):
                    if 0 <= pos < len(seq):
                        gene_map.setdefault(pos, set()).add(seq[pos])
    return out


def defining_positions(seq: str, ref: str) -> List[int]:
    seq = str(seq).strip().upper()
    ref = str(ref).strip().upper()
    n = min(len(seq), len(ref))
    return [i for i in range(n) if seq[i] != ref[i]]


def sequence_length_difference(seq: str, ref: str) -> int:
    return len(str(seq).strip().upper()) - len(str(ref).strip().upper())


def coverage_stratum(def_pos: List[int], covered: Set[int], length_difference: int = 0) -> str:
    if length_difference != 0 and not def_pos:
        return "length_difference_unresolved"
    if not def_pos:
        return "reference_like_no_defining_variant"
    if not covered:
        return "no_physical_coverage"
    n_cov = sum(1 for p in def_pos if p in covered)
    if n_cov == len(def_pos):
        return "fully_covered"
    if n_cov == 0:
        return "imputation_dependent"
    return "partially_covered"


def evidence_support_counts(seq: str, def_pos: List[int], observed_bases: Dict[int, Set[str]]) -> Dict[str, int]:
    seq = str(seq).strip().upper()
    covered = 0
    matched = 0
    clean_matched = 0
    conflicting = 0
    mixed = 0
    for pos in def_pos:
        bases = {str(b).upper() for b in observed_bases.get(pos, set()) if str(b).strip()}
        if not bases:
            continue
        covered += 1
        truth_base = seq[pos] if 0 <= pos < len(seq) else ""
        has_match = bool(truth_base and truth_base in bases)
        has_conflict = bool(not truth_base or any(base != truth_base for base in bases))
        if has_match:
            matched += 1
            if not has_conflict:
                clean_matched += 1
        if has_conflict:
            conflicting += 1
        if has_match and has_conflict:
            mixed += 1
    return {
        "covered_defining_variant_count": covered,
        "matched_covered_defining_variant_count": matched,
        "clean_matched_covered_defining_variant_count": clean_matched,
        "conflicting_covered_defining_variant_count": conflicting,
        "mixed_covered_defining_variant_count": mixed,
    }


def evidence_support_stratum(
    def_pos: List[int],
    covered_positions: Set[int],
    counts: Dict[str, int],
    length_difference: int = 0,
) -> str:
    if length_difference != 0 and not def_pos:
        return "length_difference_unresolved"
    if not def_pos:
        return "reference_like_no_defining_variant"
    covered = int(counts.get("covered_defining_variant_count", 0))
    matched = int(counts.get("matched_covered_defining_variant_count", 0))
    conflicting = int(counts.get("conflicting_covered_defining_variant_count", 0))
    if not covered_positions:
        return "no_physical_coverage"
    if covered == 0:
        return "prior_only_imputation_dependent"
    if conflicting > 0:
        return "conflicting_covered_evidence"
    if matched == len(def_pos):
        return "fully_evidence_supported"
    if matched == covered:
        return "partially_evidence_supported"
    return "covered_but_not_target_supporting"


def find_label_file(results_root: Path, expr: str, pop: str, seed: str) -> Optional[Path]:
    label_dir = results_root / "labels" / expr / pop
    if not label_dir.exists():
        return None
    matches = sorted(label_dir.glob(f"genotype_*{seed}.csv"))
    if not matches:
        matches = sorted(label_dir.glob(f"genotype_*seed{seed_number(seed)}.csv"))
    return matches[0] if matches else None


def find_evidence_file(results_root: Path, expr: str, pop: str, gene: str, seed: str) -> Optional[Path]:
    d = results_root / "mutations" / expr / pop / gene
    if not d.exists():
        return None
    patterns = [
        f"{pop}_{seed}.{gene}_sequences.csv",
        f"*{pop}*{seed}*.{gene}_sequences.csv",
        f"*seed{seed_number(seed)}*.{gene}_sequences.csv",
    ]
    for pat in patterns:
        hits = sorted(d.glob(pat))
        if hits:
            return hits[0]
    return None


def find_validation_fold(results_root: Path, expr: str, gene: str, evidence_name: str) -> Optional[int]:
    base = results_root / "validation" / expr / gene
    if not base.exists():
        return None
    for fold_dir in sorted(base.glob("fold_*")):
        if (fold_dir / evidence_name).exists():
            m = re.search(r"fold_(\d+)", fold_dir.name)
            return int(m.group(1)) if m else None
    return None
