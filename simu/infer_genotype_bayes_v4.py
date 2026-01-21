#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Robust but simple genotype inference (K=1 or K=2) for sample-level pangenome calling.
(Modified V3)

Logic: 
  - Uses robust Bayesian inference (V3 features: edge priors, pi estimation, candidate unions).
Output: 
  - Formatted strictly like V2 (gene, allele_A, seq_A, allele_B, seq_B...).
  - Extra V3 info (pi, notes, raw alleles) appended as trailing columns.

  --pangenome_dir   (expects node_support_long.csv and optionally edge_support_long.csv; graphs/<FAMILY>.gfa optional)
  --sample_csv
  --population_id
  --Kmax
  --min_naive
  --eps
  --penalty_K
  --out

Adds optional:
  --pi_min (default 0.1)  minimum mixture component for calling heterozygous
  --max_candidates (default 40) cap for candidate alleles per family
"""

from __future__ import annotations
import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, Iterable
from collections import defaultdict

import numpy as np
import pandas as pd


# ----------------------------
# Utilities
# ----------------------------

def logsumexp(xs: List[float]) -> float:
    m = max(xs)
    if not math.isfinite(m):
        return -math.inf
    s = sum(math.exp(x - m) for x in xs)
    return m + math.log(s) if s > 0 else -math.inf


def stable_id(s: str, n: int = 16) -> str:
    import hashlib
    h = hashlib.sha1(s.encode("utf-8")).hexdigest()
    return h[:n]


def parse_ranges_part(s: str) -> List[Tuple[int, int]]:
    """Parse a ranges string like '0,159;178,283' into [(0,159),(178,283)]."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return []
    s = str(s).strip()
    if s == "":
        return []
    out = []
    for part in s.split(";"):
        part = part.strip()
        if part == "":
            continue
        try:
            a, b = part.split(",")
            out.append((int(a), int(b)))
        except Exception:
            continue
    return out


def parse_observed_ranges(s: str) -> List[List[Tuple[int, int]]]:
    """
    ObservedRanges format in sample csv looks like:
      '0,159;178,283|0,284'  (split pieces separated by '|')
    Return list of list of ranges per split piece.
    """
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return []
    s = str(s).strip()
    if s == "":
        return []
    parts = s.split("|")
    return [parse_ranges_part(p) for p in parts]


def split_pipe_list(s: str) -> List[Optional[float]]:
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return []
    s = str(s).strip()
    if s == "":
        return []
    out: List[Optional[float]] = []
    for p in s.split("|"):
        p = p.strip()
        if p == "" or p.lower() == "nan":
            out.append(None)
        else:
            try:
                out.append(float(p))
            except Exception:
                out.append(None)
    return out


# ----------------------------
# Mutation parsing (0-based ref coords)
# ----------------------------

_MUT_TOKEN_RE = re.compile(r"^([SID])(.+)$")
MutationToken = Tuple[int, str]

def parse_mutation_info(mut_info: str) -> List[MutationToken]:
    """
    Parses mutation string into a simple LIST of (Pos, Token).
    CRITICAL: Does NOT sort. Preserves the exact order from the input string.
    """
    if mut_info is None or (isinstance(mut_info, float) and np.isnan(mut_info)):
        return []
    s = str(mut_info).strip()
    if s == "" or s.lower() == "germline":
        return []

    # Regex to extract individual tokens
    sub_pat = r'S(?P<sub_ref>[ATGC])(?P<sub_pos>\d+)(?P<sub_target>[ATGCN])'
    del_pat = r'D(?P<del_seq>[ATGC]+)(?P<del_pos>\d+)'
    ins_pat = r'I(?P<ins_pos>\d+)(?P<ins_seq>[ATGC]+)'
    combined_pattern = re.compile(f'({sub_pat}|{del_pat}|{ins_pat})')
    
    tokens: List[MutationToken] = []
    
    # Iterate in order found in string
    for m in combined_pattern.finditer(s):
        full_token = m.group(0)
        pos = -1
        if m.group('sub_pos'): pos = int(m.group('sub_pos'))
        elif m.group('del_pos'): pos = int(m.group('del_pos'))
        elif m.group('ins_pos'): pos = int(m.group('ins_pos'))
        
        if pos >= 0:
            tokens.append((pos, full_token))
            
    return tokens


def token_pos_span(tok: str) -> Tuple[int, int]:
    """
    Return (pos, span_len) in ref coordinates for conflict checking.
    - SNP token: S<ref><pos><alt> => span 1 at pos
    - INS token: I<pos><ins> => span 1 at pos (state at pos becomes ins+ref[pos])
    - DEL token: D<pos><seq> => span len(seq) starting at pos
    """
    if tok.startswith("I"):
        # I<pos><ins>
        m = re.match(r"^I(\d+)([A-Za-z]+)$", tok)
        if not m:
            return (-1, 0)
        pos = int(m.group(1))
        return (pos, 1)
    if tok.startswith("D"):
        m = re.match(r"^D([A-Za-z]+)(\d+)$", tok)
        if not m:
            return (-1, 0)
        pos = int(m.group(2))
        seq = m.group(1)
        return (pos, len(seq))
    if tok.startswith("S"):
        # SNP token: e.g. SC269T (ref base C, pos 269, alt T)
        m = re.match(r"^S([A-Za-z])(\d+)([A-Za-z])$", tok)
        if not m:
            return (-1, 0)
        pos = int(m.group(2))
        return (pos, 1)
    return (-1, 0)


def allele_map_from_tokens(ref: str, tokens: List[MutationToken]) -> Dict[int, str]:
    """
    Apply mutations from the list strictly in order.
    Uses CUMULATIVE update logic for insertions.
    """
    amap: Dict[int, str] = {}
    L = len(ref)
    
    tokens_by_pos = defaultdict(list)
    for pos, tok in tokens:
        tokens_by_pos[pos].append(tok)

    for pos in sorted(tokens_by_pos.keys()):
        if not (0 <= pos < L): continue

        pos_tokens = tokens_by_pos[pos]
        ref_seq = list(ref[pos].upper())
        for tok in reversed(pos_tokens):
            if tok.startswith("S"):
                m = re.match(r"^S([A-Za-z])(\d+)([A-Za-z])$", tok)
                if m and ref_seq:
                    alt = m.group(3).upper()
                    ref_seq[0] = alt
            elif tok.startswith("D"):
                m = re.match(r"^D([A-Za-z]+)(\d+)$", tok)
                if m:
                    seq = m.group(1).upper()
                    del ref_seq[0 :1]
            elif tok.startswith("I"):
                m = re.match(r"^I(\d+)([A-Za-z]+)$", tok)
                if m:
                    ins = m.group(2).upper()
                    ref_seq[0: 0] = list(ins)
        if ref_seq:
            amap[pos] = "".join(ref_seq)
        else:
            amap[pos] = "-"
    return amap


def tokens_to_mutinfo(tokens: List[MutationToken]) -> str:
    """Reconstruct mutation string from list."""
    if not tokens:
        return "Germline"
    return ";".join([t[1] for t in tokens])


def allele_state_at(ref: str, amap: Dict[int, str], pos: int) -> str:
    return amap.get(pos, ref[pos].upper())


# ----------------------------
# Reference retrieval
# ----------------------------

def parse_gfa_ref(gfa_path: Path, family: str) -> str:
    ref_bases = {}
    with gfa_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("S\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            seq = parts[2]
            tags = parts[3:]
            if f"FN:Z:{family}" not in tags:
                continue
            if "TP:Z:ref" not in tags:
                continue
            pos0 = None
            for t in tags:
                if t.startswith("P0:i:"):
                    try:
                        pos0 = int(t.split(":")[-1])
                    except:
                        pos0 = None
                    break
            if pos0 is None:
                continue
            if len(seq) != 1:
                continue
            ref_bases[pos0] = seq
    if not ref_bases:
        return ""
    L = max(ref_bases.keys()) + 1
    return "".join(ref_bases.get(i, "N") for i in range(L))


def parse_mutations(mut_string: str) -> List[dict]:
    """
        Parses mutations (Compact Format) and adjusts coordinates.
        Formats handled:
          - Substitution: ST6A (Ref T at 6 -> A)
          - Deletion:      DA28 (Ref A at 28 deleted)
          - Insertion:     I28C (Insert C at 28)
        """
    muts = []
    if not mut_string:
        return muts

    # 1. Substitution: S + Ref + Pos + Target (e.g., SC12T)
    #    Group Name: SUB
    sub_pat = r'S(?P<sub_ref>[ATGC])(?P<sub_pos>\d+)(?P<sub_target>[ATGCN])'

    # 2. Deletion: D + RefBase(s) + Pos (e.g., DA28)
    #    Group Name: DEL
    del_pat = r'D(?P<del_seq>[ATGC]+)(?P<del_pos>\d+)'

    # 3. Insertion: I + Pos + InsertedSeq (e.g., I28C)
    #    Group Name: INS
    ins_pat = r'I(?P<ins_pos>\d+)(?P<ins_seq>[ATGC]+)'

    # Combine all patterns using OR (|) operator
    # re.finditer will find matches in the order they appear in the string.
    combined_pattern = re.compile(f'({sub_pat}|{del_pat}|{ins_pat})')

    # Traverse all matches in the string to ensure strict input order
    for m in combined_pattern.finditer(mut_string):

        # Determine the mutation type based on which named group matched
        if m.group('sub_pos'):
            # --- Substitution Match ---
            raw_pos = int(m.group('sub_pos'))
            rel_pos = raw_pos

            muts.append({
                'type': 'S',
                'pos': rel_pos,
                'ref': m.group('sub_ref'),
                'target': m.group('sub_target')
            })

        elif m.group('del_pos'):
            # --- Deletion Match ---
            deleted_seq = m.group('del_seq')
            raw_pos = int(m.group('del_pos'))
            rel_pos = raw_pos

            muts.append({
                'type': 'D',
                'pos': rel_pos,
                'seq': deleted_seq
            })

        elif m.group('ins_pos'):
            # --- Insertion Match ---
            raw_pos = int(m.group('ins_pos'))
            rel_pos = raw_pos

            muts.append({
                'type': 'I',
                'pos': rel_pos,
                'seq': m.group('ins_seq')
            })

    # Position sorting disabled for safe application as it might tamper order
    # just output by default order
    return muts


def target2ref(target, mutations):
    # target + mutations -> ref
    target = list(target)
    for m in mutations:
        p = m['pos']
        if m['type'] == 'S':
            target[p] = m['ref']
        elif m['type'] == 'D':
            target[p:p] = list(m['seq'])
        elif m['type'] == 'I':
            del target[p : p + len(m['seq'])]

    return "".join(target)


def fallback_ref_from_sample(sample_df: pd.DataFrame) -> str:
    """
    If graphs/<family>.gfa not available, try to find a Germline row sequence;
    otherwise take the most frequent/longest sequence.
    """
    if sample_df.empty:
        return ""
    
    if "Sequence" not in sample_df.columns or "MutationInfo" not in sample_df.columns:
        return ""

    # 2. Get first row data (iloc[0])
    first_row = sample_df.iloc[0]
    
    # Get Target Sequence
    target_seq = str(first_row["Sequence"]).strip().upper()
    if not target_seq:
        return ""

    # Get Mutation String
    mut_string = str(first_row["MutationInfo"]).strip()

    # 3. If mutation is germline/nan/empty, target is reference
    if not mut_string or mut_string.lower() == 'nan' or mut_string.lower() == 'germline':
        return target_seq

    try:
        # 4. Parse mutation string
        mutations = parse_mutations(mut_string)
        # 5. Back-calculate reference (Target + Mutations -> Ref)
        ref_seq = target2ref(target_seq, mutations)
        return ref_seq.upper()
    except Exception as e:
        # Return empty on parse error to avoid crash
        return ""


# ----------------------------
# Priors from weights
# ----------------------------

class NodePriors:
    def __init__(self, fam: str, node_df: pd.DataFrame, pop_id: Optional[str], alpha: float = 1e-6):
        self.fam = fam
        self.alpha = float(alpha)
        df = node_df[node_df["family"].astype(str) == str(fam)].copy()
        if pop_id is not None and "population_id" in df.columns:
            sub = df[df["population_id"].astype(str) == str(pop_id)]
            if not sub.empty:
                df = sub
        # exposure per pos
        df["allele_string"] = df["allele_string"].astype(str)
        self.exposure: Dict[int, float] = {}
        self.support: Dict[Tuple[int, str], float] = {}
        self.states: Dict[int, Set[str]] = {}
        for r in df.itertuples(index=False):
            pos = int(r.pos0)
            a = str(r.allele_string)
            w = float(r.weight_sum)
            if a == "*":
                self.exposure[pos] = self.exposure.get(pos, 0.0) + w
            else:
                key = (pos, a)
                self.support[key] = self.support.get(key, 0.0) + w
                self.states.setdefault(pos, set()).add(a)

    def p(self, pos: int, state: str) -> float:
        expo = self.exposure.get(pos, 0.0)
        stset = self.states.get(pos, set())
        k = max(1, len(stset))
        sup = self.support.get((pos, state), 0.0)
        # allow unseen states with alpha mass
        return (sup + self.alpha) / (expo + self.alpha * k) if (expo + self.alpha * k) > 0 else 1.0 / k

    def argmax_state(self, pos: int, default: str) -> str:
        stset = self.states.get(pos, set())
        if not stset:
            return default
        best = None
        bestp = -1.0
        for s in stset:
            pp = self.p(pos, s)
            if pp > bestp:
                bestp = pp
                best = s
        return best if best is not None else default


class EdgePriors:
    def __init__(self, fam: str, edge_df: pd.DataFrame, pop_id: Optional[str], alpha: float = 1e-6):
        self.fam = fam
        self.alpha = float(alpha)
        df = edge_df[edge_df["family"].astype(str) == str(fam)].copy()
        if pop_id is not None and "population_id" in df.columns:
            sub = df[df["population_id"].astype(str) == str(pop_id)]
            if not sub.empty:
                df = sub
        df["from_state"] = df["from_state"].astype(str)
        df["to_state"] = df["to_state"].astype(str)
        self.exposure: Dict[Tuple[int, str], float] = {}
        self.support: Dict[Tuple[int, str, str], float] = {}
        self.to_states: Dict[Tuple[int, str], Set[str]] = {}
        for r in df.itertuples(index=False):
            pos = int(r.pos0)
            u = str(r.from_state)
            v = str(r.to_state)
            w = float(r.weight_sum)
            if v == "*":
                self.exposure[(pos, u)] = self.exposure.get((pos, u), 0.0) + w
            else:
                self.support[(pos, u, v)] = self.support.get((pos, u, v), 0.0) + w
                self.to_states.setdefault((pos, u), set()).add(v)

    def has_from(self, pos: int, u: str) -> bool:
        return (pos, u) in self.exposure

    def p(self, pos: int, u: str, v: str) -> float:
        expo = self.exposure.get((pos, u), 0.0)
        tset = self.to_states.get((pos, u), set())
        k = max(1, len(tset))
        sup = self.support.get((pos, u, v), 0.0)
        return (sup + self.alpha) / (expo + self.alpha * k) if (expo + self.alpha * k) > 0 else 1.0 / k


def allele_log_prior_markov(ref: str, amap: Dict[int, str], nodep: Optional[NodePriors], edgep: Optional[EdgePriors]) -> float:
    if nodep is None:
        return 0.0
    L = len(ref)
    if L == 0:
        return -math.inf
    lp = math.log(max(nodep.p(0, allele_state_at(ref, amap, 0)), 1e-300))
    for pos in range(0, L - 1):
        u = allele_state_at(ref, amap, pos)
        v = allele_state_at(ref, amap, pos + 1)
        if edgep is not None and edgep.has_from(pos, u):
            pp = edgep.p(pos, u, v)
        else:
            pp = nodep.p(pos + 1, v)
        lp += math.log(max(pp, 1e-300))
    return lp


# ----------------------------
# Missing-region completion (conditional MAP/Viterbi under Markov prior)
# ----------------------------

def _allowed_states_at_pos(ref: str,
                           base_amap: Dict[int, str],
                           covered: Set[int],
                           pos: int,
                           nodep: Optional[NodePriors]) -> List[str]:
    """
    Build candidate state list at a position for Viterbi filling.
    - If pos is covered: hard-constrain to the allele state implied by base_amap/ref.
    - If pos is uncovered: allow node prior states at this pos, plus ref/base state as fallback.
    """
    ref_state = ref[pos].upper()
    base_state = allele_state_at(ref, base_amap, pos)

    if pos in covered:
        return [base_state]

    stset: Set[str] = set()
    if nodep is not None:
        stset |= set(nodep.states.get(pos, set()))
    # Always allow ref and base state (important when MiXCR "fills" outside ObservedRanges)
    stset.add(ref_state)
    stset.add(base_state)

    # Stable ordering to keep runs deterministic
    return sorted(stset)


def viterbi_fill_amap(ref: str,
                      base_amap: Dict[int, str],
                      covered: Set[int],
                      nodep: Optional[NodePriors],
                      edgep: Optional[EdgePriors]) -> Dict[int, str]:
    """
    Conditional MAP completion of missing positions using the same Markov prior as allele_log_prior_markov():
      - Covered positions are hard constrained to the base_amap/ref state.
      - Uncovered positions are chosen to maximize prior probability (node+edge priors).
    Returns a new amap (pos->state string). Likelihood is still computed on covered ranges only.
    If nodep is None (NO-PRIOR mode), returns a shallow copy of base_amap (no behavior change).
    """
    if nodep is None:
        return dict(base_amap)

    L = len(ref)
    if L == 0:
        return dict(base_amap)

    allowed: List[List[str]] = [
        _allowed_states_at_pos(ref, base_amap, covered, pos, nodep) for pos in range(L)
    ]

    # DP: maps state at pos -> best log-prob up to pos
    dp: Dict[str, float] = {}
    for s in allowed[0]:
        dp[s] = math.log(max(nodep.p(0, s), 1e-300))

    back: List[Dict[str, str]] = []  # back[pos][v] = best_u at pos

    for pos in range(0, L - 1):
        dp2: Dict[str, float] = {}
        bp: Dict[str, str] = {}

        for v in allowed[pos + 1]:
            best = -math.inf
            best_u = None

            # transition from any u at pos
            for u, lp_u in dp.items():
                if edgep is not None and edgep.has_from(pos, u):
                    pp = edgep.p(pos, u, v)
                else:
                    pp = nodep.p(pos + 1, v)
                cand = lp_u + math.log(max(pp, 1e-300))
                if cand > best:
                    best = cand
                    best_u = u

            dp2[v] = best
            # best_u should exist; still guard
            bp[v] = best_u if best_u is not None else allowed[pos][0]

        back.append(bp)
        dp = dp2

    # Termination
    end_state = max(dp.items(), key=lambda kv: kv[1])[0]

    # Backtrack
    states = [None] * L
    states[L - 1] = end_state
    for pos in range(L - 2, -1, -1):
        states[pos] = back[pos][states[pos + 1]]

    # Convert to amap: only store positions where state differs from ref
    filled: Dict[int, str] = {}
    for pos, s in enumerate(states):
        if s is None:
            continue
        if s != ref[pos].upper():
            filled[pos] = s

    return filled


def amap_to_sequence(ref: str, amap: Dict[int, str]) -> str:
    """Materialize a full-length sequence from amap/ref, then remove '-' deletions."""
    out = []
    for pos in range(len(ref)):
        out.append(allele_state_at(ref, amap, pos))
    return "".join(out).replace("-", "")



# ----------------------------
# Evidence
# ----------------------------

@dataclass
class EvidenceItem:
    weight: float           # used for scoring (effective count)
    ranges: List[Tuple[int, int]]
    state_map: Dict[int, str]
    mut_info: str


def build_evidence_items(sample_df: pd.DataFrame, ref: str, min_naive: int = 2) -> Tuple[List[EvidenceItem], Set[str], Set[int]]:
    """
    Build evidence items similar to v2 but use effective counts (DiversityIndexSplit + 2*NaiveDiversityIndexSplit)
    rather than normalized weights, to make pi estimation more stable.
    """
    df = sample_df.copy()
    # same filter as v2 (NaiveDiversityIndex >= min_naive)
    if "NaiveDiversityIndex" in df.columns:
        df = df[df["NaiveDiversityIndex"].fillna(0).astype(float) >= float(min_naive)].copy()
    if df.empty:
        return [], set(), set()

    # cluster_score_raw as in v2
    df["cluster_score_raw"] = df["DiversityIndex"].fillna(0).astype(float) + 2.0 * df["NaiveDiversityIndex"].fillna(0).astype(float)

    evidences: List[EvidenceItem] = []
    mut_infos: Set[str] = set()
    covered: Set[int] = set()
    L = len(ref)

    for r in df.itertuples(index=False):
        mi = str(getattr(r, "MutationInfo", "")).strip()
        if mi == "" or mi.lower() == "nan":
            mi = "Germline"
        mut_infos.add(mi)

        state_map = allele_map_from_tokens(ref, parse_mutation_info(mi))

        # split parts
        obs_parts = str(getattr(r, "ObservedRanges", "")).split("|") if hasattr(r, "ObservedRanges") else []
        div_parts = split_pipe_list(getattr(r, "DiversityIndexSplit", "")) if hasattr(r, "DiversityIndexSplit") else []
        naive_parts = split_pipe_list(getattr(r, "NaiveDiversityIndexSplit", "")) if hasattr(r, "NaiveDiversityIndexSplit") else []

        mlen = max(len(obs_parts), len(div_parts), len(naive_parts))
        obs_parts += [""] * (mlen - len(obs_parts))
        div_parts += [None] * (mlen - len(div_parts))
        naive_parts += [None] * (mlen - len(naive_parts))

        for i in range(mlen):
            ranges = parse_ranges_part(obs_parts[i])
            if not ranges:
                continue
            d = 0.0 if div_parts[i] is None else float(div_parts[i])
            n = 0.0 if naive_parts[i] is None else float(naive_parts[i])
            # Effective evidence strength for this split
            w = max(0.0, d + 2.0 * n)
            if w <= 0:
                continue

            # mark covered
            for (a, b) in ranges:
                a2 = max(0, int(a)); b2 = min(L, int(b))
                for pos0 in range(a2, b2):
                    covered.add(pos0)

            evidences.append(EvidenceItem(weight=w, ranges=ranges, state_map=state_map, mut_info=mi))

    return evidences, mut_infos, covered


def log_p_e_given_allele(e: EvidenceItem, ref: str, amap: Dict[int, str], eps: float) -> float:
    """Log-likelihood for evidence e given allele. Compare only within evidence ranges."""
    L = len(ref)
    lp = 0.0
    log_match = math.log(max(1.0 - eps, 1e-300))
    # mismatch assumes 4-way
    log_mis = math.log(max(eps / 4.0, 1e-300))
    for (a, b) in e.ranges:
        a2 = max(0, int(a)); b2 = min(L, int(b))
        for pos in range(a2, b2):
            obs = e.state_map.get(pos, ref[pos].upper())
            tru = amap.get(pos, ref[pos].upper())
            lp += (log_match if obs == tru else log_mis)
    return lp


# ----------------------------
# Candidate generation (merge compatible mutation token sets)
# ----------------------------

def amap_signature(amap: Dict[int, str]) -> Tuple[Tuple[int, str], ...]:
    """Stable signature of an allele for deduplication."""
    return tuple(sorted(amap.items(), key=lambda kv: kv[0]))

def check_compatibility_and_merge(ref: str, toks_a: List[MutationToken], toks_b: List[MutationToken]) -> Optional[List[MutationToken]]:
    """
    1. Check if compatible (using allele map state).
    2. If compatible, merge: A + (B-A).
    3. Sort merged list by Position (Stable sort preserves relative order).
    """
    amap_a = allele_map_from_tokens(ref, toks_a)
    amap_b = allele_map_from_tokens(ref, toks_b)
    
    # Conflict check for overlap mutation pos
    pos_intersect = set(amap_a.keys()).intersection(amap_b.keys())
    for p in pos_intersect:
        if amap_a[p] != amap_b[p]:
            return None # Conflict
            
    # Merge Logic: Start with A, append unique items from B
    merged = list(toks_a)
    for t in toks_b:
        if t not in merged:
            merged.append(t)
            
    # --- CORE CHANGE: STABLE SORT BY POSITION ---
    # Maintains relative order for same-pos tokens
    merged.sort(key=lambda x: x[0])
    
    return merged


def build_candidate_token_lists(ref: str, mut_infos: Set[str], max_candidates: int = 40) -> List[List[MutationToken]]:
    """
    Build candidates using List Merge Strategy.
    Replaces the old build_candidate_token_sets.
    """
    # 1. Base candidates (parsed as lists)
    base_lists: List[List[MutationToken]] = []
    for mi in mut_infos:
        base_lists.append(parse_mutation_info(mi))

    # 2. Dedup base
    uniq_cands: List[List[MutationToken]] = []
    seen_sig = set()
    
    for toks in base_lists:
        amap = allele_map_from_tokens(ref, toks)
        sig = amap_signature(amap)
        if sig in seen_sig:
            continue
        seen_sig.add(sig)
        uniq_cands.append(toks)

    # 3. BFS Closure (Merge compatible lists)
    all_cands = list(uniq_cands)
    
    qi = 0
    while qi < len(all_cands) and len(all_cands) < max_candidates:
        A_toks = all_cands[qi]
        qi += 1
        
        for B_toks in uniq_cands:
            if len(all_cands) >= max_candidates:
                break
            
            # Try to merge A and B
            merged_toks = check_compatibility_and_merge(ref, A_toks, B_toks)
            
            if merged_toks is not None:
                # Check if this new allele is strictly new
                amap = allele_map_from_tokens(ref, merged_toks)
                sig = amap_signature(amap)
                if sig not in seen_sig:
                    seen_sig.add(sig)
                    all_cands.append(merged_toks)

    return all_cands


# ----------------------------
# Genotype scoring (K=1 or K=2)
# ----------------------------

def estimate_pi_grid(logp1: np.ndarray, logp2: np.ndarray, weights: np.ndarray, pi_min: float) -> Tuple[float, float]:
    """
    Estimate pi by simple grid search maximizing:
      sum_i w_i * log(pi*exp(logp1_i) + (1-pi)*exp(logp2_i))
    Returns (best_pi, best_ll)
    """
    # Use stable computation in log domain
    pis = np.linspace(0.0, 1.0, 201)
    best_pi = 0.5
    best_ll = -np.inf
    for pi in pis:
        if pi_min > 0 and (pi < pi_min or pi > 1.0 - pi_min):
            # still evaluate for likelihood
            pass
        # log(pi*p1 + (1-pi)*p2) = logsumexp([log(pi)+logp1, log(1-pi)+logp2])
        if pi <= 0.0:
            lmix = logp2
        elif pi >= 1.0:
            lmix = logp1
        else:
            a = math.log(pi) + logp1
            b = math.log(1.0 - pi) + logp2
            m = np.maximum(a, b)
            lmix = m + np.log(np.exp(a - m) + np.exp(b - m))
        ll = float(np.sum(weights * lmix))
        if ll > best_ll:
            best_ll = ll
            best_pi = float(pi)
    return best_pi, best_ll


def consensus_fill_sequence(ref: str, amap: Dict[int, str], covered: Set[int], nodep: Optional[NodePriors]) -> str:
    """
    For covered positions: use allele state.
    For uncovered: use most probable node state (consensus) if available, else ref.
    """
    L = len(ref)
    out = []
    for pos in range(L):
        if pos in covered:
            out.append(amap.get(pos, ref[pos].upper()))
        else:
            if nodep is not None:
                out.append(nodep.argmax_state(pos, ref[pos].upper()))
            else:
                out.append(ref[pos].upper())
    return "".join(out).replace("-", "")


def infer_family(fam: str, sdf: pd.DataFrame, node_df: pd.DataFrame, edge_df: Optional[pd.DataFrame],
                 pop_id: Optional[str], ref: str, args) -> Dict[str, object]:
    evidences, mut_infos, covered = build_evidence_items(sdf, ref, min_naive=args.min_naive)
    if not evidences:
        return {"Family": fam, "K": 0, "alleles": "", "score": -np.inf, "completed_sequences": "", "note": "no_evidence"}
    # priors
    nodep = None
    edgep = None
    if node_df is not None:
        nodep = NodePriors(fam, node_df, pop_id, alpha=1e-6)
        if edge_df is not None:
            edgep = EdgePriors(fam, edge_df, pop_id, alpha=1e-6)

    # candidate alleles (token sets)
    cand_token_lists = build_candidate_token_lists(ref, mut_infos, max_candidates=args.max_candidates)
    # build allele maps and priors
    cand = []
    for toks in cand_token_lists:
        amap = allele_map_from_tokens(ref, toks)
        mi = tokens_to_mutinfo(toks)
        amap_filled = viterbi_fill_amap(ref, amap, covered, nodep, edgep)
        lp = allele_log_prior_markov(ref, amap_filled, nodep, edgep)
        cand.append((mi, toks, amap, amap_filled, lp))

    # Precompute log P(e|allele) for each candidate
    W = np.array([e.weight for e in evidences], dtype=float)
    # guard
    W[W < 0] = 0.0
    N_eff = float(np.sum(W))
    if N_eff <= 0:
        return {"Family": fam, "K": 0, "alleles": "", "score": -np.inf, "completed_sequences": "", "note": "zero_weight"}

    log_e = {}  # mi -> np.array per evidence
    for mi, toks, amap, amap_filled, lp in cand:
        arr = np.array([log_p_e_given_allele(e, ref, amap, eps=args.eps) for e in evidences], dtype=float)
        log_e[mi] = arr

    # K=1 best
    best1 = None
    best1_score = -np.inf
    for mi, toks, amap, amap_filled, lp in cand:
        ll = float(np.sum(W * log_e[mi]))
        score = lp + ll
        if score > best1_score:
            best1_score = score
            best1 = (mi, amap, amap_filled, lp, ll)

    # If Kmax == 1, return pure
    if int(args.Kmax) <= 1:
        mi, amap, amap_filled, lp, ll = best1
        return {
            "Family": fam,
            "K": 1,
            "alleles": mi,
            "score": best1_score,
            "pi": "",
            "completed_sequences": (consensus_fill_sequence(ref, amap, covered, nodep) if nodep is None else amap_to_sequence(ref, amap_filled)),
            "note": "Kmax1"
        }

    # K=2 search over pairs (small candidates)
    best2_score = -np.inf
    best2 = None
    # light Occam penalty to avoid always choosing K=2 when evidence is weak
    occam = 0.5 * math.log(max(N_eff, 1.0))
    for i in range(len(cand)):
        mi1, t1, a1, a1_filled, lp1 = cand[i]
        for j in range(i + 1, len(cand)):
            mi2, t2, a2, a2_filled, lp2 = cand[j]
            # Estimate pi and hetero ll
            pi_hat, ll2 = estimate_pi_grid(log_e[mi1], log_e[mi2], W, pi_min=args.pi_min)
            if min(pi_hat, 1.0 - pi_hat) < args.pi_min:
                continue
            score2 = (lp1 + lp2) + ll2 - occam - float(args.penalty_K)
            if score2 > best2_score:
                best2_score = score2
                best2 = (mi1, mi2, a1, a2, a1_filled, a2_filled, pi_hat, lp1, lp2, ll2)

    # choose
    if best2 is not None and best2_score > best1_score:
        mi1, mi2, a1, a2, a1_filled, a2_filled, pi_hat, lp1, lp2, ll2 = best2
        if nodep is None:
            seqs = consensus_fill_sequence(ref, a1, covered, nodep) + "|" + consensus_fill_sequence(ref, a2, covered, nodep)
        else:
            seqs = amap_to_sequence(ref, a1_filled) + "|" + amap_to_sequence(ref, a2_filled)
        return {
            "Family": fam,
            "K": 2,
            "alleles": f"{mi1}|{mi2}",
            "score": best2_score,
            "pi": f"{pi_hat:.3f}|{1.0-pi_hat:.3f}",
            "completed_sequences": seqs,
            "note": "hetero"
        }
    else:
        mi, amap, amap_filled, lp, ll = best1
        return {
            "Family": fam,
            "K": 1,
            "alleles": mi,
            "score": best1_score,
            "pi": "",
            "completed_sequences": (consensus_fill_sequence(ref, amap, covered, nodep) if nodep is None else amap_to_sequence(ref, amap_filled)),
            "note": "homo"
        }


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pangenome_dir", default="")
    ap.add_argument("--sample_csv", required=True)
    ap.add_argument("--population_id", default="")
    ap.add_argument("--Kmax", type=int, default=2)
    ap.add_argument("--min_naive", type=int, default=2)
    ap.add_argument("--eps", type=float, default=0.01)
    ap.add_argument("--penalty_K", type=float, default=0.0,
                    help="Extra penalty for K=2 vs K=1 (kept for CLI compatibility; default 0.0).")
    ap.add_argument("--out", required=True)

    # optional but safe
    ap.add_argument("--pi_min", type=float, default=0.1)
    ap.add_argument("--max_candidates", type=int, default=40)

    args = ap.parse_args()

    if args.pangenome_dir:
        pan = Path(args.pangenome_dir)
        node_path = pan / "node_support_long.csv"
        edge_path = pan / "edge_support_long.csv"
        gfa_dir = pan / "graphs"

        if not node_path.exists():
            raise FileNotFoundError(f"Missing node weights: {node_path}")

        node_df = pd.read_csv(node_path)
        edge_df = pd.read_csv(edge_path) if edge_path.exists() else None
    else:
        print("Wait! Running in NO-PRIOR mode (Likelihood only).")
        node_df = None
        edge_df = None
        gfa_dir = None

    sample = pd.read_csv(args.sample_csv)

    pop_id = args.population_id.strip()
    pop_id = pop_id if pop_id != "" else None

    families = sorted(sample["Family"].dropna().astype(str).unique().tolist())
    out_rows = []

    # cap Kmax to 2 for this simple solver
    if args.Kmax > 2:
        args.Kmax = 2

    for fam in families:
        sdf = sample[sample["Family"].astype(str) == str(fam)].copy()
        # reference
        if gfa_dir is None:
            ref = fallback_ref_from_sample(sdf)
        else:
            gfa = gfa_dir / f"{fam}.gfa"
            ref = parse_gfa_ref(gfa, fam) if gfa.exists() else fallback_ref_from_sample(sdf)
        if ref == "":
            # Error row
            out_rows.append({
                "gene": fam, "score": -np.inf, 
                "allele_A": "Error: No Ref", 
                "v3_note": "no_ref"
            })
            continue

        res = infer_family(fam, sdf, node_df, edge_df, pop_id, ref, args)
        
        # --- Transform V3 output logic to V2 CSV format ---
        
        final_row = {
            "gene": res["Family"],
            "score": res["score"]
        }

        # Raw output from infer_family is pipe-delimited
        al_list = str(res["alleles"]).split("|")
        sq_list = str(res["completed_sequences"]).split("|")
        
        is_error = (res["score"] == -np.inf)

        # Handle Diploid Display:
        # If result is valid but K=1, duplicate into A and B (like A=B)
        if not is_error and len(al_list) == 1:
            al_list = al_list * 2
            sq_list = sq_list * 2

        labels = ["A", "B", "C", "D"] # Support polyploid slots if needed
        
        for i, label in enumerate(labels):
            if i >= len(al_list):
                break
            
            mut_str = al_list[i].strip()
            seq_str = sq_list[i] if i < len(sq_list) else ""
            
            # V2 Naming Convention
            if is_error:
                name = "N/A"
            elif not mut_str or mut_str.lower() == "germline":
                name = f"{fam}*01"
            else:
                name = f"{fam}*01_mut"
            
            final_row[f"allele_{label}"] = name
            final_row[f"seq_{label}"] = seq_str

        # Add V3 specific extras at the end
        final_row["v3_pi"] = res.get("pi", "")
        final_row["v3_note"] = res.get("note", "")
        final_row["v3_K"] = res.get("K", 0)
        final_row["v3_alleles_exact"] = res.get("alleles", "")
        
        out_rows.append(final_row)

    out_df = pd.DataFrame(out_rows)
    
    # Enforce Column Order: gene, allele_A, seq_A, allele_B, seq_B, score, [extras]
    desired_cols = ["gene", "allele_A", "seq_A", "allele_B", "seq_B", "score"]
    # Check if we have C/D
    for extra in ["C", "D"]:
        if f"allele_{extra}" in out_df.columns:
            desired_cols.extend([f"allele_{extra}", f"seq_{extra}"])
    
    # Append V3 specific columns
    v3_cols = [c for c in out_df.columns if c.startswith("v3_")]
    desired_cols.extend(v3_cols)
    
    # Filter to only existing columns
    final_cols = [c for c in desired_cols if c in out_df.columns]
    
    out_df = out_df[final_cols]
    out_df.to_csv(args.out, index=False)
    print(f"Wrote inferred genotypes (V3 logic, V2 format) to: {args.out}")


if __name__ == "__main__":
    main()