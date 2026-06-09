#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
degrade_fastq_pairs_subsample.py

Function:
- Subsample paired-end FASTQ and simulate degradation by truncation.
- Default: Truncate R1 (5' end). Optional: Truncate R2 (default 3' end).
- Output: Degraded R1/R2 FASTQ + Processing log CSV.
"""

import argparse
import csv
import gzip
import random
import sys
from typing import Tuple, Iterator

def open_text_maybe_gz(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def read_fastq(handle) -> Iterator[Tuple[str, str, str, str]]:
    """Yield (header, seq, plus, qual) without trailing newlines."""
    while True:
        h = handle.readline()
        if not h:
            break
        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        if not q:
            raise ValueError("Incomplete FASTQ format: missing lines at the end")
        yield h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")

def norm_read_id(header: str) -> str:
    if header.startswith("@"):
        return header[1:].split()[0]
    return header.split()[0]

def sample_truncnorm(mean: float, sd: float, low: float, high: float, rng: random.Random, max_tries: int = 5000) -> float:
    """Simple rejection sampling for truncated normal distribution."""
    if sd <= 0:
        return float(min(max(mean, low), high))
    for _ in range(max_tries):
        x = rng.gauss(mean, sd)
        if low <= x <= high:
            return x
    x = rng.gauss(mean, sd)
    return float(min(max(x, low), high))

def trim_seq_qual(seq: str, qual: str, trim: int, direction: str) -> Tuple[str, str]:
    if trim <= 0:
        return seq, qual
    if trim >= len(seq):
        return "", ""
    if direction == "5":
        return seq[trim:], qual[trim:]
    if direction == "3":
        return seq[:-trim], qual[:-trim]
    raise ValueError(f"Unknown direction={direction} (Use 5 or 3)")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-1", "--r1", required=True, help="R1 FASTQ(.gz)")
    ap.add_argument("-2", "--r2", required=True, help="R2 FASTQ(.gz)")
    ap.add_argument("-o1", "--out_r1", required=True, help="Output R1 FASTQ(.gz)")
    ap.add_argument("-o2", "--out_r2", required=True, help="Output R2 FASTQ(.gz)")
    ap.add_argument("--log", required=True, help="Output log CSV")

    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--keep_prob", type=float, default=1.0, help="Probability to keep a read-pair (0-1).")
    ap.add_argument("--p_degraded", type=float, default=0.7, help="Probability of degradation among kept pairs.")

    # Truncation distribution parameters
    ap.add_argument("--cut_intact_mean", type=float, default=0.0)
    ap.add_argument("--cut_intact_sd", type=float, default=0.0)
    ap.add_argument("--cut_max_intact", type=int, default=0)

    ap.add_argument("--cut_degraded_mean", type=float, default=70.0)
    ap.add_argument("--cut_degraded_sd", type=float, default=12.0)
    ap.add_argument("--cut_max", type=int, default=200)
    ap.add_argument("--cut_min", type=int, default=0)

    # Trim settings
    ap.add_argument("--trim_r1", action="store_true", help="Enable R1 truncation")
    ap.add_argument("--r1_dir", choices=["5", "3"], default="5")
    ap.add_argument("--r1_fixed_trim", type=int, default=0)

    ap.add_argument("--trim_r2", action="store_true", help="Enable R2 truncation")
    ap.add_argument("--r2_dir", choices=["5", "3"], default="3")
    ap.add_argument("--r2_fixed_trim", type=int, default=0)

    ap.add_argument("--min_read_len", type=int, default=80)
    ap.add_argument("--max_pairs", type=int, default=0)

    args = ap.parse_args()

    if not (0.0 <= args.keep_prob <= 1.0):
        raise ValueError("--keep_prob must be in [0,1]")
    if not (0.0 <= args.p_degraded <= 1.0):
        raise ValueError("--p_degraded must be in [0,1]")

    if not args.trim_r1 and not args.trim_r2:
        args.trim_r1 = True

    rng = random.Random(args.seed)

    fields = [
        "pair_index", "read_id", "kept_after_subsample",
        "state", "d_rand",
        "r1_trim_total", "r1_dir", "r1_len_before", "r1_len_after",
        "r2_trim_total", "r2_dir", "r2_len_before", "r2_len_after",
        "dropped", "drop_reason"
    ]

    total = 0
    kept_after_subsample = 0
    kept_final = 0
    dropped_subsample = 0
    dropped_short = 0
    degraded_n = 0
    intact_n = 0

    with open_text_maybe_gz(args.r1, "rt") as f1, open_text_maybe_gz(args.r2, "rt") as f2, \
         open_text_maybe_gz(args.out_r1, "wt") as o1, open_text_maybe_gz(args.out_r2, "wt") as o2, \
         open(args.log, "w", newline="") as flog:

        writer = csv.DictWriter(flog, fieldnames=fields)
        writer.writeheader()

        for idx, (rec1, rec2) in enumerate(zip(read_fastq(f1), read_fastq(f2))):
            if args.max_pairs and idx >= args.max_pairs:
                break
            total += 1

            h1, s1, p1, q1 = rec1
            h2, s2, p2, q2 = rec2
            rid1 = norm_read_id(h1)
            rid2 = norm_read_id(h2)
            rid = rid1 if rid1 == rid2 else rid1

            r1_before = len(s1)
            r2_before = len(s2)

            # Step 0: Subsampling
            if rng.random() >= args.keep_prob:
                dropped_subsample += 1
                writer.writerow({
                    "pair_index": idx, "read_id": rid, "kept_after_subsample": 0,
                    "state": "", "d_rand": 0,
                    "r1_trim_total": 0, "r1_dir": "", "r1_len_before": r1_before, "r1_len_after": 0,
                    "r2_trim_total": 0, "r2_dir": "", "r2_len_before": r2_before, "r2_len_after": 0,
                    "dropped": 1, "drop_reason": "subsampled_out",
                })
                continue

            kept_after_subsample += 1

            # Step 1: Decide degradation state
            state = "degraded" if (rng.random() < args.p_degraded) else "intact"
            if state == "degraded":
                degraded_n += 1
                d = int(round(sample_truncnorm(
                    args.cut_degraded_mean, args.cut_degraded_sd,
                    args.cut_min, args.cut_max, rng
                )))
            else:
                intact_n += 1
                d = int(round(sample_truncnorm(
                    args.cut_intact_mean, args.cut_intact_sd,
                    args.cut_min, args.cut_max_intact, rng
                )))

            # Step 2: Apply truncation
            r1_trim_total = 0
            r2_trim_total = 0

            if args.trim_r1:
                r1_trim_total = max(0, args.r1_fixed_trim) + max(0, d)
                s1, q1 = trim_seq_qual(s1, q1, r1_trim_total, args.r1_dir)

            if args.trim_r2:
                r2_trim_total = max(0, args.r2_fixed_trim) + max(0, d)
                s2, q2 = trim_seq_qual(s2, q2, r2_trim_total, args.r2_dir)

            r1_after = len(s1)
            r2_after = len(s2)

            # Step 3: Filter by minimum length
            if (r1_after < args.min_read_len) or (r2_after < args.min_read_len):
                dropped_short += 1
                writer.writerow({
                    "pair_index": idx, "read_id": rid, "kept_after_subsample": 1,
                    "state": state, "d_rand": d,
                    "r1_trim_total": r1_trim_total, "r1_dir": args.r1_dir if args.trim_r1 else "",
                    "r1_len_before": r1_before, "r1_len_after": r1_after,
                    "r2_trim_total": r2_trim_total, "r2_dir": args.r2_dir if args.trim_r2 else "",
                    "r2_len_before": r2_before, "r2_len_after": r2_after,
                    "dropped": 1, "drop_reason": "too_short_after_trim",
                })
                continue

            # Step 4: Write output
            kept_final += 1
            o1.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
            o2.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")

            writer.writerow({
                "pair_index": idx, "read_id": rid, "kept_after_subsample": 1,
                "state": state, "d_rand": d,
                "r1_trim_total": r1_trim_total, "r1_dir": args.r1_dir if args.trim_r1 else "",
                "r1_len_before": r1_before, "r1_len_after": r1_after,
                "r2_trim_total": r2_trim_total, "r2_dir": args.r2_dir if args.trim_r2 else "",
                "r2_len_before": r2_before, "r2_len_after": r2_after,
                "dropped": 0, "drop_reason": "",
            })

    print("DONE")
    print(f"Total pairs read:             {total}")
    print(f"Kept after subsample:         {kept_after_subsample} ({(kept_after_subsample/total if total else 0):.2%}), target={args.keep_prob}")
    print(f"  Degraded among kept:        {degraded_n} ({(degraded_n/kept_after_subsample if kept_after_subsample else 0):.2%}), target={args.p_degraded}")
    print(f"  Intact among kept:          {intact_n} ({(intact_n/kept_after_subsample if kept_after_subsample else 0):.2%})")
    print(f"Dropped by subsampling out:  {dropped_subsample} ({(dropped_subsample/total if total else 0):.2%})")
    print(f"Dropped by too-short trim:   {dropped_short} ({(dropped_short/total if total else 0):.2%})")
    print(f"Final kept pairs written:    {kept_final} ({(kept_final/total if total else 0):.2%})")
    print(f"Out R1: {args.out_r1}")
    print(f"Out R2: {args.out_r2}")
    print(f"Log:    {args.log}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)