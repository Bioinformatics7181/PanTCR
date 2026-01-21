#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
simulate_bias_advanced.py

Generates a molecular pool FASTA file based on a Repertoire CSV, implementing 
"continuous 3' bias / 5' truncation + post-J C-region padding (prefix from C-start)".
Also outputs 'simulation_log.csv' for traceability, detailing truncation, 
padding, and length parameters for each molecule.

Core Features:
- Hierarchical Sampling (Optional): Support for sample-level parameters + molecule-level continuous sampling.
- Molecule-level States: Two states (Intact vs. Degraded). Truncation length 'd' follows a Truncated Normal distribution (continuous, rather than discrete steps).
- Length Control: Target length 'L' follows a Truncated Normal distribution (narrow); the final length is strictly equal to 'L_final'.
- C-Region Padding: Prefixes are taken strictly from the start of the C-gene sequence to respect the authenticity of the J->C junction. Supports 'c_need=0' (no padding).
- Meta Logging: Comprehensive logs that can pass through ground-truth columns from the original CSV.

Input CSV Requirements:
  Must contain at least: 'clone_id', 'full_sequence', 'read_count'.

Usage Example:
  python simulate_bias_advanced.py -i rep.csv -o molecules.fasta --log simulation_log.csv \
    --leader_trim_len 55 --p_degraded 0.7 --p_degraded_sd 0.08 --sample_col sample_id \
    --len_mean 340 --len_sd 6 --len_min 330 --len_max 360 \
    --cut_degraded_mean 70 --cut_degraded_sd 20 --cut_max 140 \
    --c_max 80 --min_keep_len 200 --seed 42
"""

import argparse
import csv
import hashlib
import random
import sys
from typing import Dict, Any, Tuple, List, Optional

import pandas as pd


# ==========================================
# TRBC
# ==========================================
TRBC1_1_SEQ = "AGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTTCCCTGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACGGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATGGTCAAGAGAAAGGATTTCTGA"
TRBC2_1_SEQ = "AGGACCTGAAAAACGTGTTCCCACCCAAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTACCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCTGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTCACCTCCGAGTCTTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCTTGCTAGGGAAGGCCACCTTGTATGCCGTGCTGGTCAGTGCCCTCGTGCTGATGGCCATGGTCAAGAGAAAGGATTCCAGAGGCTAG"
C_SEQS: List[str] = [TRBC1_1_SEQ, TRBC2_1_SEQ]


# ===================== Tools =====================
def stable_hash_int(s: str) -> int:
    h = hashlib.md5(s.encode("utf-8")).hexdigest()
    return int(h[:8], 16)


def sample_truncnorm(mean: float, sd: float, low: float, high: float, rng: random.Random, max_tries: int = 5000) -> float:
    if sd <= 0:
        return float(min(max(mean, low), high))
    for _ in range(max_tries):
        x = rng.gauss(mean, sd)
        if low <= x <= high:
            return x
    x = rng.gauss(mean, sd)
    return float(min(max(x, low), high))


def clamp_int(x: int, low: int, high: int) -> int:
    return max(low, min(high, x))


def clamp_float(x: float, low: float, high: float) -> float:
    return max(low, min(high, x))


# ===================== Main =====================
def generate_one_molecule(
    base_seq: str,
    rng: random.Random,
    params: Dict[str, Any],
    sample_params: Dict[str, Any],
) -> Tuple[str, Dict[str, Any]]:
    """
    Return:
      final_seq: final sequence
      meta: logs for checking
    """
    base_len = len(base_seq)

    if params["len_min"] < params["min_keep_len"]:
        raise ValueError(f"len_min({params['len_min']}) must be >= min_keep_len({params['min_keep_len']})")

    p_deg = sample_params["p_degraded_sample"]

    u_state = rng.random()
    state = "degraded" if u_state < p_deg else "intact"

    L_drawn = int(round(sample_truncnorm(
        params["len_mean"], params["len_sd"], params["len_min"], params["len_max"], rng
    )))
    L = clamp_int(L_drawn, params["len_min"], params["len_max"])

    if state == "intact":
        cut_drawn = int(round(sample_truncnorm(
            params["cut_intact_mean"], params["cut_intact_sd"],
            params["cut_min"], params["cut_max_intact"], rng
        )))
    else:
        cut_drawn = int(round(sample_truncnorm(
            sample_params["cut_degraded_mean_sample"], params["cut_degraded_sd"],
            params["cut_min"], params["cut_max"], rng
        )))

    max_cut_allowed = max(0, base_len - params["min_keep_len"])
    d = clamp_int(cut_drawn, 0, max_cut_allowed)

    cut_total = d
    seq_cut = base_seq[cut_total:]
    kept_len_preC = len(seq_cut)

    extra_cut = 0
    clamped = False
    clamp_reason: Optional[str] = None
    resampled = 0

    if kept_len_preC > L:
        extra_cut = kept_len_preC - L
        cut_total = clamp_int(cut_total + extra_cut, 0, max_cut_allowed)
        seq_cut = base_seq[cut_total:]
        kept_len_preC = len(seq_cut)
        if kept_len_preC > L:
            clamped = True
            clamp_reason = "cannot_cut_enough_to_reach_L"
            L = kept_len_preC

    c_need = L - kept_len_preC

    while c_need > params["c_max"] and resampled < params["max_resample"]:
        resampled += 1

        L_drawn2 = int(round(sample_truncnorm(
            params["len_mean"], params["len_sd"], params["len_min"], params["len_max"], rng
        )))
        L2 = clamp_int(L_drawn2, params["len_min"], params["len_max"])

        if state == "intact":
            cut_drawn2 = int(round(sample_truncnorm(
                params["cut_intact_mean"], params["cut_intact_sd"],
                params["cut_min"], params["cut_max_intact"], rng
            )))
        else:
            cut_drawn2 = int(round(sample_truncnorm(
                sample_params["cut_degraded_mean_sample"], params["cut_degraded_sd"],
                params["cut_min"], params["cut_max"], rng
            )))
        d2 = clamp_int(cut_drawn2, 0, max_cut_allowed)

        cut_total2 = d2
        seq_cut2 = base_seq[cut_total2:]
        kept2 = len(seq_cut2)

        extra_cut2 = 0
        if kept2 > L2:
            extra_cut2 = kept2 - L2
            cut_total2 = clamp_int(cut_total2 + extra_cut2, 0, max_cut_allowed)
            seq_cut2 = base_seq[cut_total2:]
            kept2 = len(seq_cut2)
            if kept2 > L2:
                L2 = kept2

        c_need2 = L2 - kept2

        if c_need2 <= params["c_max"]:
            L = L2
            L_drawn = L_drawn2
            cut_drawn = cut_drawn2
            cut_total = cut_total2
            seq_cut = seq_cut2
            kept_len_preC = kept2
            extra_cut = extra_cut2
            c_need = c_need2
            break

    if c_need > params["c_max"]:
        clamped = True
        clamp_reason = "c_need_exceeds_c_max"
        L = kept_len_preC + params["c_max"]
        c_need = params["c_max"]

    c_idx = -1
    c_fill = ""
    if c_need > 0:
        c_idx = rng.randrange(len(C_SEQS))
        c_src = C_SEQS[c_idx]
        if len(c_src) < c_need:
            c_fill = (c_src * (c_need // len(c_src) + 1))[:c_need]
        else:
            c_fill = c_src[:c_need]

    final_seq = seq_cut + c_fill
    L_final = len(final_seq)

    meta = {
        "state": state,
        "u_state": u_state,
        "p_degraded_sample": p_deg,
        "cut_degraded_mean_sample": sample_params["cut_degraded_mean_sample"],

        "base_len": base_len,
        "L_drawn": L_drawn,
        "L_final": L_final,

        "cut_drawn": cut_drawn,
        "cut_total": cut_total,
        "extra_cut": extra_cut,

        "kept_len_preC": kept_len_preC,
        "c_need": c_need,
        "c_idx": (c_idx + 1) if c_idx >= 0 else 0,

        "resampled": resampled,
        "clamped": int(clamped),
        "clamp_reason": clamp_reason or "",
    }
    return final_seq, meta


def build_sample_params(
    sample_key: str,
    global_rng: random.Random,
    params: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Draw sample-level parameters for each sample to introduce global variance 
    between samples while maintaining molecule-level heterogeneity within each sample.
    
    - p_degraded_sample: Slight fluctuation around the global p_degraded.
    - cut_degraded_mean_sample: Slight fluctuation around the global cut_degraded_mean.
    """
    sub_seed = params["seed"] + stable_hash_int(sample_key)
    rng = random.Random(sub_seed)

    p_deg = sample_truncnorm(
        params["p_degraded"], params["p_degraded_sd"], 0.0, 1.0, rng
    )
    p_deg = clamp_float(p_deg, 0.0, 1.0)

    cut_deg_mean = sample_truncnorm(
        params["cut_degraded_mean"], params["cut_degraded_mean_sd"],
        params["cut_min"], params["cut_max"], rng
    )
    cut_deg_mean = clamp_float(cut_deg_mean, params["cut_min"], params["cut_max"])

    return {
        "sample_key": sample_key,
        "p_degraded_sample": float(p_deg),
        "cut_degraded_mean_sample": float(cut_deg_mean),
        "sample_seed": sub_seed,
    }


# ===================== Main =====================
def main() -> None:
    ap = argparse.ArgumentParser(description="CSV->FASTA: 5' continuous truncation + post-J C-region padding (prefix), with comprehensive meta logging.")

    ap.add_argument("-i", "--input_csv", required=True, help="Input repertoire.csv")
    ap.add_argument("-o", "--out_fasta", required=True, help="Output molecules.fasta")
    ap.add_argument("--log", dest="out_meta_csv", required=True, help="Output simulation_log.csv")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")

    ap.add_argument("--sample_col", type=str, default="", help="Column name for Sample ID in CSV (e.g., 'sample_id'). If empty, treated as a single sample.")

    ap.add_argument("--leader_trim_len", type=int, default=0, help="Initial 5' leader trimming length (set to 0 if full_sequence does not contain a leader)")

    # Length profile parameters
    ap.add_argument("--len_mean", type=float, default=340.0, help="Mean target length")
    ap.add_argument("--len_sd", type=float, default=6.0, help="Standard deviation of target length")
    ap.add_argument("--len_min", type=int, default=330, help="Minimum allowed target length")
    ap.add_argument("--len_max", type=int, default=360, help="Maximum allowed target length")

    ap.add_argument("--c_max", type=int, default=80, help="Maximum length allowed for C-region padding")

    ap.add_argument("--min_keep_len", type=int, default=200, help="Minimum sequence length to keep after truncation")

    # Degraded state parameters
    ap.add_argument("--p_degraded", type=float, default=0.7, help="Mean global proportion of degraded molecules")
    ap.add_argument("--p_degraded_sd", type=float, default=0.0, help="Inter-sample fluctuation (SD) of p_degraded. 0 means no sample-level variance.")

    # Truncation/Cut parameters
    ap.add_argument("--cut_min", type=int, default=0, help="Absolute minimum truncation length")
    ap.add_argument("--cut_max", type=int, default=140, help="Absolute maximum truncation length")

    # Intact molecules truncation profile
    ap.add_argument("--cut_max_intact", type=int, default=30, help="Max truncation length for 'intact' molecules")
    ap.add_argument("--cut_intact_mean", type=float, default=5.0, help="Mean truncation for 'intact' molecules")
    ap.add_argument("--cut_intact_sd", type=float, default=5.0, help="SD of truncation for 'intact' molecules")

    # Degraded molecules truncation profile
    ap.add_argument("--cut_degraded_mean", type=float, default=70.0, help="Mean truncation for 'degraded' molecules")
    ap.add_argument("--cut_degraded_sd", type=float, default=20.0, help="SD of truncation for 'degraded' molecules")
    ap.add_argument("--cut_degraded_mean_sd", type=float, default=0.0, help="Inter-sample fluctuation (SD) of cut_degraded_mean. 0 means no sample-level variance.")

    ap.add_argument("--max_resample", type=int, default=5, help="Maximum number of resampling attempts when required C-padding (c_need) exceeds c_max")

    args = ap.parse_args()
    global_rng = random.Random(args.seed)

    df = pd.read_csv(args.input_csv)
    required = {"clone_id", "full_sequence", "read_count"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing mandatory columns in input CSV file: {sorted(missing)}")

    sample_col = args.sample_col.strip()
    if sample_col and sample_col not in df.columns:
        raise ValueError(f"--sample_col={sample_col} but no this column in the csv file")

    passthrough_cols = [c for c in [
        "sample_id", "v_gene", "j_gene", "v_allele", "j_allele", "cdr3", "frequency", "method"
    ] if c in df.columns]

    if sample_col and sample_col not in passthrough_cols:
        passthrough_cols.insert(0, sample_col)

    params = {
        "seed": int(args.seed),
        "sample_col": sample_col,
        "leader_trim_len": int(args.leader_trim_len),

        "len_mean": float(args.len_mean),
        "len_sd": float(args.len_sd),
        "len_min": int(args.len_min),
        "len_max": int(args.len_max),

        "c_max": int(args.c_max),
        "min_keep_len": int(args.min_keep_len),

        "p_degraded": float(args.p_degraded),
        "p_degraded_sd": float(args.p_degraded_sd),

        "cut_min": int(args.cut_min),
        "cut_max": int(args.cut_max),

        "cut_max_intact": int(args.cut_max_intact),
        "cut_intact_mean": float(args.cut_intact_mean),
        "cut_intact_sd": float(args.cut_intact_sd),

        "cut_degraded_mean": float(args.cut_degraded_mean),
        "cut_degraded_sd": float(args.cut_degraded_sd),
        "cut_degraded_mean_sd": float(args.cut_degraded_mean_sd),

        "max_resample": int(args.max_resample),
    }

    if sample_col:
        sample_keys = df[sample_col].astype(str).fillna("NA").unique().tolist()
    else:
        sample_keys = ["__SINGLE_SAMPLE__"]

    sample_param_map: Dict[str, Dict[str, Any]] = {}
    for sk in sample_keys:
        sample_param_map[sk] = build_sample_params(sk, global_rng, params)

    total_mols = int(df["read_count"].sum())
    print(f"[INFO] clones={len(df)}, total_molecules={total_mols}, seed={args.seed}")
    if sample_col:
        print(f"[INFO] sample_col={sample_col}, n_samples={len(sample_keys)}")
        pvals = [sample_param_map[sk]["p_degraded_sample"] for sk in sample_keys]
        print(f"[INFO] p_degraded_sample range: {min(pvals):.3f} .. {max(pvals):.3f}")

    # meta
    meta_cols = [
        "clone_id", "mol_index",
        "raw_len", "leader_trim_len",
        "base_len",

        *passthrough_cols,

        "sample_key",
        "sample_seed",
        "p_degraded_sample",
        "cut_degraded_mean_sample",

        "state",
        "u_state",
        "L_drawn", "L_final",
        "cut_drawn",
        "cut_total",
        "cut_total_raw",  # cut_total + leader_trim_len
        "extra_cut",
        "kept_len_preC",
        "c_need",
        "c_idx",
        "resampled",
        "clamped",
        "clamp_reason",
    ]

    severe_count = 0

    with open(args.out_fasta, "w") as f_fa, open(args.out_meta_csv, "w", newline="") as f_log:
        w = csv.DictWriter(f_log, fieldnames=meta_cols)
        w.writeheader()

        for _, row in df.iterrows():
            clone_id = str(row["clone_id"])
            raw_seq = str(row["full_sequence"]).strip().upper()
            raw_len = len(raw_seq)
            n = int(row["read_count"])
            if n <= 0:
                continue

            # leader trim
            if args.leader_trim_len > 0 and raw_len > args.leader_trim_len:
                base_seq = raw_seq[args.leader_trim_len:]
            else:
                base_seq = raw_seq
            base_len = len(base_seq)

            # sample key
            if sample_col:
                sk = str(row[sample_col]) if pd.notna(row[sample_col]) else "NA"
            else:
                sk = "__SINGLE_SAMPLE__"
            sp = sample_param_map[sk]

            for i in range(n):
                final_seq, meta = generate_one_molecule(base_seq, global_rng, params, sp)
                if meta["state"] == "degraded":
                    severe_count += 1

                header = (
                    f">{clone_id}|mol={i}|st={meta['state']}"
                    f"|L={meta['L_final']}|cut={meta['cut_total']}|C={meta['c_need']}"
                )
                f_fa.write(header + "\n")
                f_fa.write(final_seq + "\n")

                out_row: Dict[str, Any] = {
                    "clone_id": clone_id,
                    "mol_index": i,
                    "raw_len": raw_len,
                    "leader_trim_len": args.leader_trim_len,
                    "base_len": base_len,
                    "sample_key": sp["sample_key"],
                    "sample_seed": sp["sample_seed"],
                    "cut_total_raw": meta["cut_total"] + args.leader_trim_len,
                }

                for c in passthrough_cols:
                    if c in row:
                        out_row[c] = row[c]
                    else:
                        out_row[c] = ""

                out_row.update(meta)
                w.writerow(out_row)

    print("[DONE] Simulation Complete.")
    print(f"  Total Molecules: {total_mols}")
    print(f"  Degraded Ratio:  {severe_count/total_mols:.2%} (global mean target: {args.p_degraded})")
    print(f"  Output FASTA:    {args.out_fasta}")
    print(f"  Log CSV:         {args.out_meta_csv}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
