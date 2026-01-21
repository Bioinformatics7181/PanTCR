#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_pangenome_graph_v7.py

Build a pangenome graph from per-sample allele CSVs and a metadata table.
This version preserves v6 behavior and additionally records edge exposure/support
in edge_support_long.csv. Output files are written under --out_dir.
"""

import argparse
import sys
import gc
import hashlib
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def stable_id(s: str, n: int = 16) -> str:
    return hashlib.sha1(s.encode("utf-8")).hexdigest()[:n]

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

_TOKEN_RE = re.compile(
    r"(S[ACGTN]\d+[ACGTN])"
    r"|"
    r"(I\d+[ACGTN]+?(?=S|I|D|$))"
    r"|"
    r"(D\d+[ACGTN]+?(?=S|I|D|$))"
)

def tokenize_mutation_info(mut_info: str) -> List[str]:
    if mut_info is None:
        return []
    s = str(mut_info).strip()
    if s == "" or s.lower() == "nan":
        return []
    parts = [p.strip() for p in s.split(";") if p.strip()]
    toks: List[str] = []
    if len(parts) > 1:
        for p in parts:
            toks.extend([m.group(0) for m in _TOKEN_RE.finditer(p) if m.group(0)])
        return toks
    return [m.group(0) for m in _TOKEN_RE.finditer(s) if m.group(0)]

def parse_mutations(mut_string: str) -> List[dict]:
    """Parse compact mutation annotations into a list of edit operations."""
    muts = []
    if not mut_string:
        return muts

    sub_pat = r'S(?P<sub_ref>[ATGC])(?P<sub_pos>\d+)(?P<sub_target>[ATGCN])'

    del_pat = r'D(?P<del_seq>[ATGC]+)(?P<del_pos>\d+)'

    ins_pat = r'I(?P<ins_pos>\d+)(?P<ins_seq>[ATGC]+)'

    # Combine all patterns using OR (|) operator
    # re.finditer will find matches in the order they appear in the string.
    combined_pattern = re.compile(f'({sub_pat}|{del_pat}|{ins_pat})')

    # Iterate over all matches in input order
    for m in combined_pattern.finditer(mut_string):

        # Determine the mutation type based on which named group matched
        if m.group('sub_pos'):
            raw_pos = int(m.group('sub_pos'))
            rel_pos = raw_pos

            muts.append({
                'type': 'S',
                'pos': rel_pos,
                'ref': m.group('sub_ref'),
                'target': m.group('sub_target')
            })

        elif m.group('del_pos'):
            deleted_seq = m.group('del_seq')
            raw_pos = int(m.group('del_pos'))
            rel_pos = raw_pos

            muts.append({
                'type': 'D',
                'pos': rel_pos,
                'seq': deleted_seq
            })

        elif m.group('ins_pos'):
            raw_pos = int(m.group('ins_pos'))
            rel_pos = raw_pos

            muts.append({
                'type': 'I',
                'pos': rel_pos,
                'seq': m.group('ins_seq')
            })

    # Sort by position descending for safe application
    # Do not sort by position to preserve the original S/I ordering at identical positions
    return muts

def split_preserve_empty(s: str, sep: str) -> List[str]:
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return [""]
    return str(s).split(sep)

def parse_ranges_part(part: str) -> List[Tuple[int, int]]:
    """Return list of (start,end) with end exclusive semantics."""
    if part is None:
        return []
    p = str(part).strip()
    if p == "":
        return []
    segs = []
    for seg in p.split(";"):
        seg = seg.strip()
        if not seg:
            continue
        ab = seg.split(",")
        if len(ab) != 2:
            continue
        try:
            a = int(ab[0]); b = int(ab[1])
        except:
            continue
        if b < a:
            a, b = b, a
        segs.append((a, b))
    return segs

def parse_split_numbers(s: str) -> List[Optional[float]]:
    parts = split_preserve_empty(s, "|")
    out: List[Optional[float]] = []
    for p in parts:
        p2 = str(p).strip()
        if p2 == "":
            out.append(None)
        else:
            try:
                out.append(float(p2))
            except:
                out.append(None)
    return out

def load_metadata(meta_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(meta_path)
    required = {"sample_id", "population_id"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata is missing required columns: {sorted(missing)}")
    if "filename" not in meta.columns:
        meta["filename"] = meta["sample_id"].astype(str) + ".csv"
    return meta

def read_sample_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)

def compute_cluster_scores(df: pd.DataFrame, sample_id: str, population_id: str, min_naive: int) -> pd.DataFrame:
    out = df.copy()
    out["sample_id"] = sample_id
    out["population_id"] = population_id
    out = out[out["NaiveDiversityIndex"].fillna(0).astype(float) >= float(min_naive)].copy()
    if out.empty:
        return out
    out["allele_id"] = out.apply(lambda r: stable_id(f'{r["Family"]}|{r["Sequence"]}'), axis=1)
    out["cluster_score_raw"] = out["DiversityIndex"].fillna(0).astype(float) + 2.0 * out["NaiveDiversityIndex"].fillna(0).astype(float)
    denom = out.groupby(["Family"])["cluster_score_raw"].transform("sum").replace(0, np.nan)
    out["cluster_score_norm"] = (out["cluster_score_raw"] / denom).fillna(0.0)
    return out

def target2ref(target: str, mutations: List[dict]) -> str:
    target = list(str(target))
    for m in mutations:
        p = m["pos"]
        if m["type"] == "S":
            if 0 <= p < len(target):
                target[p] = m["ref"]
        elif m["type"] == "D":
            if 0 <= p <= len(target):
                target[p:p] = list(m["seq"])
        elif m["type"] == "I":
            if 0 <= p <= len(target):
                del target[p:p + len(m["seq"])]
    return "".join(target)

def allele_string_map_from_mut(ref: str, mut_info: str) -> Dict[int, str]:
    L = len(ref)
    state = {i: ref[i] for i in range(L)}
    muts = parse_mutations(mut_info)

    # S/D first
    for m in muts:
        p = m["pos"]
        if m["type"] == "S":
            if 0 <= p < L:
                state[p] = m["target"]
        elif m["type"] == "D":
            for j in range(len(m["seq"])):
                if 0 <= p + j < L:
                    state[p + j] = "-"

    # I then (prepend)
    for m in muts:
        if m["type"] != "I":
            continue
        p = m["pos"]
        ins = m["seq"]
        if 0 <= p < L:
            state[p] = ins + state[p]
        elif p == L and L > 0:
            # best-effort for end insertion
            state[L - 1] = state[L - 1] + ins

    return state



def tuple_list_from_ref(ref: str, mutations: List[dict]) -> List[Tuple[Optional[int], str]]:
    seq = [(i, b) for i, b in enumerate(ref)]
    for m in mutations[::-1]:
        p = m["pos"]
        if m["type"] == "S":
            for j, (ri, _) in enumerate(seq):
                if ri == p:
                    seq[j] = (ri, m["target"])
                    break
        elif m["type"] == "D":
            L = len(m["seq"])
            seq = [(ri, b) for (ri, b) in seq if not (ri is not None and p <= ri < p + L)]
        elif m["type"] == "I":
            ins = [(None, b) for b in m["seq"]]
            if p >= len(ref):
                seq.extend(ins)
            else:
                insert_at = None
                for j, (ri, _) in enumerate(seq):
                    if ri == p:
                        insert_at = j
                        break
                if insert_at is None:
                    seq.extend(ins)
                else:
                    seq[insert_at:insert_at] = ins
    return seq

def write_family_gfa(family: str, ref: str, fam_alleles: pd.DataFrame, out_gfa: Path, write_paths: bool) -> None:
    lines = ["H\tVN:Z:1.0"]
    edges = set()
    Lref = len(ref)

    for i, b in enumerate(ref):
        nid = f"{family}_R_{i}"
        lines.append(f"S\t{nid}\t{b}\tFN:Z:{family}\tTP:Z:ref\tP0:i:{i}")
        if i > 0:
            edges.add((f"{family}_R_{i-1}", nid))

    ins_nodes_written = set()
    snp_nodes_written = set()
    allele_paths = {}

    for _, r in fam_alleles.iterrows():
        allele_id = r["allele_id"]
        muts = parse_mutations(r["MutationInfo"])
        tuples = tuple_list_from_ref(ref, muts)

        path_nodes = []
        i = 0
        while i < len(tuples):
            ri, b = tuples[i]
            if ri is None:
                j = i
                ins_chars = []
                while j < len(tuples) and tuples[j][0] is None:
                    ins_chars.append(tuples[j][1])
                    j += 1
                ins_seq = "".join(ins_chars)

                next_ref = None
                for k in range(j, len(tuples)):
                    if tuples[k][0] is not None:
                        next_ref = tuples[k][0]
                        break
                anchor0 = int(next_ref) if next_ref is not None else Lref

                hid = stable_id(f"{anchor0}|{ins_seq}", 10)
                nid = f"{family}_I_{anchor0}_{hid}"
                if nid not in ins_nodes_written:
                    lines.append(f"S\t{nid}\t{ins_seq}\tFN:Z:{family}\tTP:Z:ins\tA0:i:{anchor0}\tIS:Z:{ins_seq}")
                    ins_nodes_written.add(nid)
                path_nodes.append(nid)
                i = j
                continue

            pos0 = int(ri)
            if 0 <= pos0 < Lref and b == ref[pos0]:
                path_nodes.append(f"{family}_R_{pos0}")
            else:
                nid = f"{family}_S_{pos0}_{b}"
                if 0 <= pos0 < Lref and nid not in snp_nodes_written:
                    lines.append(f"S\t{nid}\t{b}\tFN:Z:{family}\tTP:Z:snp\tP0:i:{pos0}\tRB:Z:{ref[pos0]}\tBB:Z:{b}")
                    snp_nodes_written.add(nid)
                path_nodes.append(nid)
            i += 1

        for a, b in zip(path_nodes[:-1], path_nodes[1:]):
            edges.add((a, b))
        allele_paths[allele_id] = [n + "+" for n in path_nodes]

    for (a, b) in sorted(edges):
        lines.append(f"L\t{a}\t+\t{b}\t+\t0M")

    if write_paths:
        for aid, nodes in allele_paths.items():
            if nodes:
                lines.append(f"P\t{family}_{aid}\t{','.join(nodes)}\t*")

    out_gfa.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_naive", type=int, default=1)
    ap.add_argument("--write_paths", action="store_true")
    ap.add_argument("--skip_graph", action="store_true")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    ensure_dir(out_dir)
    graphs_dir = out_dir / "graphs"
    ensure_dir(graphs_dir)

    meta = load_metadata(Path(args.metadata))

    stat_csv_path = out_dir / "sample_allele_stat.csv"
    stat_columns = [
        "sample_id","population_id","source_file","allele_id","Family","ReferenceAllele",
        "CloneCount","DiversityIndex","NaiveDiversityIndex","cluster_score_raw","cluster_score_norm",
        "MutationInfo","Sequence","ObservedRanges",
        "CloneCountSplit","DiversityIndexSplit","NaiveDiversityIndexSplit"
    ]
    pd.DataFrame([], columns=stat_columns).to_csv(stat_csv_path, index=False)

    catalog_dict = {}
    fam_refs = defaultdict(list)
    allele_aggr = defaultdict(lambda: {"samples": set(), "pops": set(), "w_sum": 0.0})
    parse_warnings = []

    # Warn if CSV files exist in in_dir but are not listed in metadata (unused for graph building).
    meta_files = set(meta['filename'].astype(str))
    dir_csv_files = sorted([p.name for p in in_dir.glob('*.csv')])
    meta_basename = Path(args.metadata).name
    unlisted = [fn for fn in dir_csv_files if (fn not in meta_files) and (fn != meta_basename)]
    if len(unlisted) == 0:
        print(f"[INFO] 0 sample CSV files in {in_dir} were unused (not listed in metadata).", file=sys.stderr)
    else:
        print(f"[INFO] {len(unlisted)} sample CSV files in {in_dir} were unused (not listed in metadata). They will be written to parse_warnings.csv.", file=sys.stderr)
        for fn in unlisted:
            parse_warnings.append({'warning_type': 'unlisted_in_metadata', 'filename': fn})

    files_to_process = []
    for _, m in meta.iterrows():
        sample_id = str(m["sample_id"])
        pop_id = str(m["population_id"])
        filename = str(m["filename"])
        fpath = in_dir / filename
        if not fpath.exists():
            continue
        files_to_process.append((sample_id, pop_id, filename, fpath))

        sdf = read_sample_csv(fpath)
        sdf = compute_cluster_scores(sdf, sample_id, pop_id, args.min_naive)
        if sdf.empty:
            continue
        sdf["source_file"] = filename
        sdf[stat_columns].to_csv(stat_csv_path, mode="a", header=False, index=False)

        for r in sdf.itertuples(index=False):
            aid = getattr(r, "allele_id")
            fam = str(getattr(r, "Family"))
            if aid not in catalog_dict:
                catalog_dict[aid] = {
                    "allele_id": aid,
                    "GeneType": getattr(r, "GeneType") if hasattr(r, "GeneType") else "",
                    "Family": fam,
                    "ReferenceAllele": getattr(r, "ReferenceAllele") if hasattr(r, "ReferenceAllele") else "",
                    "Sequence": getattr(r, "Sequence"),
                    "MutationInfo": getattr(r, "MutationInfo") if getattr(r, "MutationInfo") is not None else "",
                    "Status": getattr(r, "Status") if hasattr(r, "Status") else "",
                }
            allele_aggr[aid]["samples"].add(sample_id)
            allele_aggr[aid]["pops"].add(pop_id)
            allele_aggr[aid]["w_sum"] += float(getattr(r, "cluster_score_norm"))

            fam_refs[fam].append((str(getattr(r, "Sequence")), str(getattr(r, "MutationInfo")) if getattr(r, "MutationInfo") is not None else ""))

            obs = split_preserve_empty(getattr(r, "ObservedRanges") if hasattr(r, "ObservedRanges") else "", "|")
            ccs = split_preserve_empty(getattr(r, "CloneCountSplit") if hasattr(r, "CloneCountSplit") else "", "|")
            if len(obs) != len(ccs):
                parse_warnings.append({"sample_id": sample_id, "Family": fam, "len_obs": len(obs), "len_clone": len(ccs)})

        del sdf
        gc.collect()

    pd.DataFrame(parse_warnings).to_csv(out_dir / "parse_warnings.csv", index=False)

    if not catalog_dict:
        pd.DataFrame([], columns=["note"]).to_csv(out_dir / "allele_catalog.csv", index=False)
        pd.DataFrame([], columns=["note"]).to_csv(out_dir / "node_support_long.csv", index=False)
        pd.DataFrame([], columns=["note"]).to_csv(out_dir / "edge_support_long.csv", index=False)
        return

    catalog_df = pd.DataFrame(list(catalog_dict.values()))

    stats_rows = []
    for aid in catalog_df["allele_id"].tolist():
        ag = allele_aggr.get(aid, None)
        if ag is None:
            stats_rows.append({"allele_id": aid, "n_samples": 0, "n_pops": 0, "support_weight": 0.0})
        else:
            stats_rows.append({"allele_id": aid, "n_samples": len(ag["samples"]), "n_pops": len(ag["pops"]), "support_weight": float(ag["w_sum"])})
    catalog_df = catalog_df.merge(pd.DataFrame(stats_rows), on="allele_id", how="left")
    catalog_df.to_csv(out_dir / "allele_catalog.csv", index=False)

    ref_map = {}
    ref_inc_rows = []
    for fam, items in fam_refs.items():
        refs_found = []
        for seq, mi in items[:200]:
            ref_seq = target2ref(seq, parse_mutations(mi))
            if ref_seq:
                refs_found.append(ref_seq)
        if not refs_found:
            ref_map[fam] = ""
            continue
        vc = pd.Series(refs_found).value_counts()
        ref_map[fam] = vc.index[0]
        if len(vc) > 1:
            ref_inc_rows.append({"Family": fam, "n_ref_variants": int(len(vc))})
    pd.DataFrame(ref_inc_rows).to_csv(out_dir / "ref_inconsistency.csv", index=False)

    # precompute allele state maps
    allele_state_map: Dict[str, Dict[int, str]] = {}
    for r in catalog_df.itertuples(index=False):
        fam = str(getattr(r, "Family"))
        ref = ref_map.get(fam, "")
        if ref == "":
            continue
        aid = getattr(r, "allele_id")
        mi = getattr(r, "MutationInfo") if getattr(r, "MutationInfo") is not None else ""
        allele_state_map[aid] = allele_string_map_from_mut(ref, mi)

    # aggregate: (pop,fam,pos0,allele_string) -> [w_sum, c_sum, n_samples]
    pop_acc = {}

    # aggregate: (pop,fam,pos0,from_state,to_state) -> [w_sum, c_sum, n_samples]
    pop_edge_acc = {}

    stat_all = pd.read_csv(stat_csv_path)
    for (sample_id, pop_id, fam), sdf in stat_all.groupby(["sample_id","population_id","Family"]):
        ref = ref_map.get(fam, "")
        if ref == "":
            continue
        Lref = len(ref)

        local = {}  # (pos0, allele_string) -> [w_sum, c_sum]

        local_edge = {}  # (pos0, from_state, to_state) -> [w_sum, c_sum]

        for r in sdf.itertuples(index=False):
            allele_id = getattr(r, "allele_id")
            cluster_w = float(getattr(r, "cluster_score_norm"))
            if cluster_w == 0.0:
                continue
            state_map = allele_state_map.get(allele_id, None)
            if state_map is None:
                continue

            obs_parts = split_preserve_empty(getattr(r, "ObservedRanges"), "|")
            div_parts = parse_split_numbers(getattr(r, "DiversityIndexSplit"))
            naive_parts = parse_split_numbers(getattr(r, "NaiveDiversityIndexSplit"))

            mlen = max(len(obs_parts), len(div_parts), len(naive_parts))
            obs_parts = obs_parts + [""] * (mlen - len(obs_parts))
            div_parts = div_parts + [None] * (mlen - len(div_parts))
            naive_parts = naive_parts + [None] * (mlen - len(naive_parts))

            seq_scores = []
            for i in range(mlen):
                d = 0.0 if div_parts[i] is None else float(div_parts[i])
                n = 0.0 if naive_parts[i] is None else float(naive_parts[i])
                seq_scores.append(d + 2.0 * n)
            sum_seq = sum(s for s in seq_scores if s > 0)

            for i in range(mlen):
                if sum_seq > 0:
                    frac = (seq_scores[i] / sum_seq) if seq_scores[i] > 0 else 0.0
                else:
                    frac = 1.0 / float(mlen) if mlen > 0 else 0.0
                alloc = cluster_w * frac
                if alloc == 0.0:
                    continue

                ranges = parse_ranges_part(obs_parts[i])
                if not ranges:
                    # Empty segments contribute to allocation but do not emit positions
                    continue

                for (a, b) in ranges:
                    a2 = max(0, int(a))
                    b2 = min(Lref, int(b))
                    if b2 <= a2:
                        continue
                    for pos0 in range(a2, b2):  # [start,end)
                        # exposure
                        k_exp = (pos0, "*")
                        if k_exp in local:
                            local[k_exp][0] += alloc
                            local[k_exp][1] += 1
                        else:
                            local[k_exp] = [alloc, 1]
                        # support
                        allele_str = state_map.get(pos0, ref[pos0])
                        k_sup = (pos0, allele_str)
                        if k_sup in local:
                            local[k_sup][0] += alloc
                            local[k_sup][1] += 1
                        else:
                            local[k_sup] = [alloc, 1]

                    # edge exposure/support for adjacent positions within this range
                    # note: only counts edges where both pos0 and pos0+1 are covered by the same contiguous range
                    for pos0_e in range(a2, b2 - 1):
                        from_state = state_map.get(pos0_e, ref[pos0_e])
                        to_state = state_map.get(pos0_e + 1, ref[pos0_e + 1])
                        # exposure conditioned on from_state
                        k_eexp = (pos0_e, from_state, "*")
                        if k_eexp in local_edge:
                            local_edge[k_eexp][0] += alloc
                            local_edge[k_eexp][1] += 1
                        else:
                            local_edge[k_eexp] = [alloc, 1]
                        # support for specific transition
                        k_esup = (pos0_e, from_state, to_state)
                        if k_esup in local_edge:
                            local_edge[k_esup][0] += alloc
                            local_edge[k_esup][1] += 1
                        else:
                            local_edge[k_esup] = [alloc, 1]


        for (pos0, allele_str), (w_sum, c_sum) in local.items():
            k = (str(pop_id), str(fam), int(pos0), str(allele_str))
            if k in pop_acc:
                pop_acc[k][0] += float(w_sum)
                pop_acc[k][1] += int(c_sum)
                pop_acc[k][2] += 1
            else:
                pop_acc[k] = [float(w_sum), int(c_sum), 1]

        for (pos0, from_state, to_state), (w_sum, c_sum) in local_edge.items():
            k = (str(pop_id), str(fam), int(pos0), str(from_state), str(to_state))
            if k in pop_edge_acc:
                pop_edge_acc[k][0] += float(w_sum)
                pop_edge_acc[k][1] += int(c_sum)
                pop_edge_acc[k][2] += 1
            else:
                pop_edge_acc[k] = [float(w_sum), int(c_sum), 1]

        del sdf
        del local
        del local_edge
        gc.collect()

    out_rows = []
    for (pop_id, fam, pos0, allele_str), (w_sum, c_sum, n_samp) in pop_acc.items():
        out_rows.append({
            "population_id": pop_id,
            "family": fam,
            "pos0": int(pos0),
            "allele_string": allele_str,
            "weight_sum": float(w_sum),
            "count_sum": int(c_sum),
            "n_samples": int(n_samp),
        })
    pd.DataFrame(out_rows).to_csv(out_dir / "node_support_long.csv", index=False)

    edge_rows = []
    for (pop_id, fam, pos0, from_state, to_state), (w_sum, c_sum, n_samp) in pop_edge_acc.items():
        edge_rows.append({
            "population_id": pop_id,
            "family": fam,
            "pos0": int(pos0),
            "from_state": from_state,
            "to_state": to_state,
            "weight_sum": float(w_sum),
            "count_sum": int(c_sum),
            "n_samples": int(n_samp),
        })
    pd.DataFrame(edge_rows).to_csv(out_dir / "edge_support_long.csv", index=False)

    if not args.skip_graph:
        for fam, fam_df in catalog_df.groupby("Family"):
            ref = ref_map.get(fam, "")
            if ref == "":
                continue
            out_gfa = graphs_dir / f"{fam}.gfa"
            write_family_gfa(fam, ref, fam_df, out_gfa, args.write_paths)

    print("Done.")


if __name__ == "__main__":
    main()
