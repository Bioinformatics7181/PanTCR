#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TCR Allele & Variant Analysis Pipeline (Expert Version V2.3)
------------------------------------------------------------
Updates:
1. Parameterized Quality Threshold (--minQual).
2. Log Output: Includes FULL sequence used for reference selection.
3. Strict Coordinate Parsing & V/J Split Alignment Logic.
4. Comprehensive Auditing.
"""

# Compared to V5, it handles split multi-mutation info and accommodates observed sites.

import pandas as pd
import re
import argparse
import sys
import os
from collections import Counter, defaultdict

# ==========================================
# 1. Configuration & Utils
# ==========================================

VALID_AA = set('ACDEFGHIKLMNPQRSTVWY')
VALID_NT = set('ATCG')

v_index = {
    "TRBV1": 69,
    "TRBV10-1": 77,
    "TRBV10-2": 77,
    "TRBV10-3": 77,
    "TRBV11-1": 77,
    "TRBV11-2": 77,
    "TRBV11-3": 77,
    "TRBV12-1": 74,
    "TRBV12-2": 74,
    "TRBV12-3": 77,
    "TRBV12-4": 77,
    "TRBV12-5": 77,
    "TRBV13": 107,
    "TRBV14": 77,
    "TRBV15": 77,
    "TRBV16": 77,
    "TRBV17": 77,
    "TRBV18": 77,
    "TRBV19": 77,
    "TRBV2": 77,
    "TRBV20-1": 62,
    "TRBV21-1": 78,
    "TRBV22-1": 77,
    "TRBV23-1": 77,
    "TRBV24-1": 77,
    "TRBV25-1": 77,
    "TRBV26": 77,
    "TRBV27": 77,
    "TRBV28": 77,
    "TRBV29-1": 65,
    "TRBV3-1": 77,
    "TRBV3-2": 77,
    "TRBV30": 71,
    "TRBV4-1": 77,
    "TRBV4-2": 77,
    "TRBV4-3": 77,
    "TRBV5-1": 77,
    "TRBV5-2": 77,
    "TRBV5-3": 77,
    "TRBV5-4": 77,
    "TRBV5-5": 77,
    "TRBV5-6": 77,
    "TRBV5-7": 77,
    "TRBV5-8": 77,
    "TRBV6-1": 77,
    "TRBV6-2": 77,
    "TRBV6-3": 77,
    "TRBV6-4": 77,
    "TRBV6-5": 77,
    "TRBV6-6": 77,
    "TRBV6-7": 77,
    "TRBV6-8": 77,
    "TRBV6-9": 77,
    "TRBV7-1": 77,
    "TRBV7-2": 77,
    "TRBV7-3": 77,
    "TRBV7-4": 77,
    "TRBV7-5": 77,
    "TRBV7-6": 77,
    "TRBV7-7": 77,
    "TRBV7-8": 77,
    "TRBV7-9": 77,
    "TRBV8": 74,
    "TRBV9": 77
}

def str2bool(v):
    """
    Used for parsing boolean values in command line arguments.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def natural_sort_key(s):
    """Sorts gene names naturally (TRBV2 before TRBV10)."""
    if pd.isna(s): return (0, "")
    return tuple([int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))])

def load_reference(filepath):
    """
    Loads reference sequences from FASTA or TSV.
    Returns: {GeneFamily: {AlleleName: Sequence_String}}
    """
    seqs = {}
    try:
        if filepath.endswith(('.fasta', '.fa')):
            curr_h, curr_s = None, []
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if curr_h:
                            fam = curr_h.split('*')[0]
                            if fam not in seqs: seqs[fam] = {}
                            seqs[fam][curr_h] = "".join(curr_s).replace('.', '').upper()
                        parts = line[1:].split('|')
                        curr_h = parts[1] if len(parts) > 1 else line[1:]
                        curr_s = []
                    else:
                        curr_s.append(line)
                if curr_h:
                    fam = curr_h.split('*')[0]
                    if fam not in seqs: seqs[fam] = {}
                    seqs[fam][curr_h] = "".join(curr_s).replace('.', '').upper()
        else: # TSV
            df = pd.read_csv(filepath, sep='\t', low_memory=False)
            df.columns = [c.strip().lower() for c in df.columns]
            for _, row in df.iterrows():
                allele = str(row['allele']).strip()
                seq = str(row['sequence']).upper().replace('.', '')
                if not allele or not seq: continue
                fam = allele.split('*')[0]
                if fam not in seqs: seqs[fam] = {}
                seqs[fam][allele] = seq
    except Exception as e:
        print(f"[Error] Failed to load reference {filepath}: {e}")
        sys.exit(1)
    return seqs

def check_protein_validity(aa_seq):
    if pd.isna(aa_seq) or len(str(aa_seq)) == 0: return False
    return set(str(aa_seq).upper()).issubset(VALID_AA)

def check_dna_validity(dna_seq):
    """Checks if sequence contains only A, T, C, G."""
    if pd.isna(dna_seq) or len(str(dna_seq)) == 0: return False
    return set(str(dna_seq).upper()).issubset(VALID_NT)

def parse_alignment_string(align_str):
    """
    Parses MiXCR alignment string.
    Format example: 62|352|375|0|290|ST64A|...
    Indices: 
      0: refStart (on Ref)
      1: refEnd (on Ref)
      2: refLen (Total Ref Len context)
      3: targetStart (on Read)
      4: targetEnd (on Read)
      5: MutationString
    """
    if pd.isna(align_str) or str(align_str).strip() == "":
        return None
    
    align_str = str(align_str).strip()

    blocks_raw = [b.strip() for b in align_str.split(",") if b.strip()]    

    min_refStart = None
    max_refEnd = None
    ref_len = None
    t_start = None
    t_end = None
    muts_str = ""
    observed_ranges = []

    for b in blocks_raw:
        parts = b.split('|')
        if len(parts) < 5:
            return None
        
        try:
            refStart =  int(parts[0])
            refEnd =  int(parts[1])
            refLen =  int(parts[2])
            tStart = int(parts[3])
            tEnd = int(parts[4])
            mutStr = parts[5] if len(parts) > 5 else ""
        except ValueError:
            return None

        if ref_len is None:
            ref_len = refLen
        elif ref_len != refLen:
            return None

        t_start = tStart if t_start is None else t_start
        t_end = tEnd if t_end is None else t_end
        min_refStart = refStart if min_refStart is None else min(min_refStart, refStart)
        max_refEnd   = refEnd   if max_refEnd   is None else max(max_refEnd, refEnd)
        muts_str += mutStr    
        observed_ranges.append((refStart, refEnd))

    if not observed_ranges:
        return None

    return {
        'refStart': min_refStart, 
        'refEnd': max_refEnd, 
        'refLen': ref_len, 
        'tStart': t_start, 
        'tEnd': t_end, 
        'mutStr': muts_str,
        'ranges': observed_ranges
    }


def parse_mutations(mut_string, ref_start, ref_end, offset):
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
    
    # Traverse all matches in the string to guarantee strict input order
    for m in combined_pattern.finditer(mut_string):
        
        # Determine the mutation type based on which named group matched
        if m.group('sub_pos'):
            # --- Substitution Match ---
            raw_pos = int(m.group('sub_pos'))
            if raw_pos < ref_start or raw_pos >= ref_end:
                continue
            rel_pos = raw_pos - ref_start - offset
            
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
            if raw_pos < ref_start or raw_pos >= ref_end:
                continue
            rel_pos = raw_pos - ref_start - offset
            
            muts.append({
                'type': 'D', 
                'pos': rel_pos, 
                'seq': deleted_seq
            })
            
        elif m.group('ins_pos'):
            # --- Insertion Match ---
            raw_pos = int(m.group('ins_pos'))
            if raw_pos < ref_start or raw_pos >= ref_end:
                continue
            rel_pos = raw_pos - ref_start - offset
            
            muts.append({
                'type': 'I', 
                'pos': rel_pos, 
                'seq': m.group('ins_seq')
            })
            
    # Sort by position descending for safe application
    # Sorting disabled to avoid tampering with default order of S/I at same position
    # Just output in the default order
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

def ref2target(ref, mutations):
    # ref + mutations -> target
    ref = list(ref)
    # Use reverse order for transformation
    mutations_desc = mutations[::-1]
    for m in mutations_desc:
        p = m['pos']
        if m['type'] == 'S':
            ref[p] = m['target']
        elif m['type'] == 'D':
            del ref[p : p + len(m['seq'])]
        elif m['type'] == 'I':
            ref[p: p] = list(m['seq'])

    return "".join(ref)

# ==========================================
# 2. Core Analyzer Class
# ==========================================

class VariantAnalyzer:
    def __init__(self, input_file, gene_type, ref_file, minClns, minQual, min_count, filter_productive):
        self.gene_type = gene_type.upper()
        self.minQual = int(minQual) # Store quality threshold parameter
        self.min_count = int(min_count)
        self.filter_productive = filter_productive # Boolean switch

        if self.gene_type not in ['V', 'J']:
            raise ValueError("Gene type must be 'V' or 'J'")
            
        print(f"Loading Data: {input_file}...")
        self.df = pd.read_csv(input_file, sep='\t', low_memory=False)
        clone_count = len(self.df)

        if clone_count < minClns:
            print("\n" + "="*60)
            print(f"[CRITICAL WARNING] Low Clone Count Detected: {clone_count}, fewer than {minClns} clones.")
            print(f"Statistical analysis of somatic hypermutations requires sufficient depth.")
            print("="*60 + "\n")
            print("[ABORT] Exiting program to prevent unreliable results.")
            sys.exit(0) # Exit cleanly
        
        print(f"Loading Reference: {ref_file}...")
        self.ref_dict = load_reference(ref_file)
        
        if self.gene_type == 'V':
            self.cols = {'hit': 'allVHitsWithScore', 'align': 'allVAlignments', 'pair': 'allJHitsWithScore', 'pairAlign': 'allJAlignments'}
        else:
            self.cols = {'hit': 'allJHitsWithScore', 'align': 'allJAlignments', 'pair': 'allVHitsWithScore', 'pairAlign': 'allVAlignments'}
            
        self.stats = {
            'total_input': len(self.df),
            'filters': defaultdict(int),
            'final_count': 0,
            'variant_types': defaultdict(int),
            'mismatch_counts': 0
        }
        
        self.ref_selection_log = [] 
        self.long_read_logs = []    
        self.bounds_logs = []
        self.error_detail_logs = []

    def preprocess_and_filter(self):
        """Split Multi-locus and Apply QC Filters."""
        print(f"Preprocessing and filtering (MinQual={self.minQual})...")
        df = self.df.copy()
        init_cnt = len(df)
        
        # 1. ReadCount
        if 'readCount' in df.columns:
            df = df[df['readCount'] >= self.min_count]
        self.stats['filters'][f'ReadCount < {self.min_count}'] = init_cnt - len(df)
        curr = len(df)

        # 2. Inference Existence
        df = df[df[self.cols['hit']].notna() & df[self.cols['pair']].notna()]
        self.stats['filters']['Missing V/J Inference'] += (curr - len(df))
        curr = len(df)

        # 3. MinQual (Parameterized)
        regions = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
        qual_cols = [f'minQual{r}' for r in regions]
        mask_keep = pd.Series(True, index=df.index)
        for col in qual_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                # Use self.minQual variable
                bad = (df[col].notna()) & (df[col] < self.minQual)
                mask_keep = mask_keep & (~bad)
        df = df[mask_keep]
        # Update statistics key to reflect actual threshold used
        self.stats['filters'][f'Region Quality < {self.minQual}'] = curr - len(df)
        curr = len(df)

        # 4. Valid Target Sequence
        # if 'targetSequences' in df.columns:
        #      # Check if targetSequences only contains ATCG
        #      valid_dna = df['targetSequences'].apply(check_dna_validity)
        #      df = df[valid_dna]
        # self.stats['filters']['Ambiguous DNA (Contains N/Other)'] = curr - len(df)
        # curr = len(df)

        # 5. Productive CDR3
        if 'aaSeqCDR3' in df.columns and self.filter_productive:
            valid_cdr3 = df['aaSeqCDR3'].apply(check_protein_validity)
            df = df[valid_cdr3]
        self.stats['filters']['Non-Productive (Invalid AA)'] = curr - len(df)
        curr = len(df)

        # 6. Split Multi-Locus
        # MiXCR logic: Hits are comma-separated, Alignments are semicolon-separated.
        # We must pair them strictly (Hit_1 <-> Align_1, Hit_2 <-> Align_2).
        
        hit_col = self.cols['hit']
        align_col = self.cols['align']
        
        # A. Split strings into lists
        # Hits: "TRBV6-2(score),TRBV6-3(score)" -> ["TRBV6-2(score)", "TRBV6-3(score)"]
        hits_split = df[hit_col].astype(str).str.split(',')
        
        # Alignments: "Align1;Align2" -> ["Align1", "Align2"]
        # Note: Using regex=False for speed, assumes standard MiXCR format
        aligns_split = df[align_col].astype(str).str.split(';')
        
        # B. Zip them together to preserve the index relationship
        # Creates a temporary column of tuples: [(Hit1, Align1), (Hit2, Align2)]
        # using list comprehension is faster than apply for this operation
        df['paired_data'] = [
            list(zip(h, a)) if len(h) == len(a) else list(zip(h, h)) # Fallback safe-guard (rare)
            for h, a in zip(hits_split, aligns_split)
        ]
        
        # Filter out cases where length mismatch might have caused data loss (optional but safe)
        # In standard MiXCR output, len(hits) == len(aligns) should always hold.
        # If mismatch, the 'zip' above truncates to the shortest list automatically.
        
        # C. Explode the paired data
        df_exploded = df.explode('paired_data')
        
        # Remove empty entries if any
        df_exploded = df_exploded.dropna(subset=['paired_data'])
        
        # D. Unpack the tuples back into specific columns
        # Index 0 is Hit, Index 1 is Alignment
        df_exploded['current_hit_str'] = df_exploded['paired_data'].apply(lambda x: x[0])
        df_exploded[align_col] = df_exploded['paired_data'].apply(lambda x: x[1])  # CRITICAL: Overwrite alignment column
        
        # E. Clean Gene Names
        # Remove score: "TRBV6-2*00(3591)" -> "TRBV6-2*00"
        df_exploded['AnalyzedGene'] = df_exploded['current_hit_str'].apply(lambda x: x.split('(')[0])
        
        # Extract Family: "TRBV6-2*00" -> "TRBV6-2" (or "TRBV6" depending on strictness, usually split by *)
        df_exploded['AnalyzedFamily'] = df_exploded['AnalyzedGene'].apply(lambda x: x.split('*')[0])
        
        self.stats['filters']['Multi-Locus Expansion (Added)'] = len(df_exploded) - len(df)
        
        # Update self.df to the exploded version
        self.df = df_exploded
        self.stats['final_count'] = len(self.df)

    def resolve_reference(self, group_df, candidates, family_name):
        """
        Selects reference allele using biologically anchored alignment (V=Left, J=Right).
        """
        best_row = None
        max_len = -1
        anchor_info = None 
        
        # A. Find Anchor (Based on target sequence length ONLY if clean alignment)
        for _, row in group_df.iterrows():
            info = parse_alignment_string(str(row[self.cols['align']]))
            # Improvement: Strict screening, must have alignment info to serve as reference Anchor
            if not info: continue 

            if not check_dna_validity(row['targetSequences']): continue
            
            curr_len = info['refEnd'] - info['refStart']

            if self.gene_type == 'V':
                if info['refStart'] < v_index[family_name]:
                    curr_len = info['refEnd'] - v_index[family_name]
            
            if curr_len > max_len:
                max_len = info['refEnd']
                best_row = row
                anchor_info = {
                    'seq': row['targetSequences'][info['tStart']:info['tEnd']], 
                    'id': row['cloneId'],
                    'mut': info['mutStr'],
                    'rstart': info['refStart'],
                    'rend': info['refEnd'],
                    'rlen': info['refLen']
                }

        if not anchor_info:
            return candidates[0] 
            
        # B. Match Candidates (V=Left-Align, J=Right-Align)
        selected_allele = None
        perfect_matches = []

        # Process mutation sequence
        # Obtain reference sequence based on target sequence and mutation info

        raw_mutations = parse_mutations(anchor_info['mut'], anchor_info['rstart'], anchor_info['rend'], 0)

        anchor_seq = target2ref(anchor_info['seq'], raw_mutations)

        if self.gene_type == 'V':
            if anchor_info['rstart'] < v_index[family_name]:
                # If target contains Leader sequence, excision
                anchor_seq = anchor_seq[v_index[family_name] - anchor_info['rstart']:]
         
        for cand in candidates:
            if cand.split('*')[0] not in self.ref_dict: continue
            ref_full = self.ref_dict[cand.split('*')[0]][cand]

            if self.gene_type == 'V':
                if anchor_info['rstart'] > v_index[family_name]:
                    ref_full = ref_full[anchor_info['rstart'] - v_index[family_name]: ]

            else:
                if anchor_info['rend'] - anchor_info['rlen'] < 0:
                    ref_full = ref_full[: anchor_info['rend'] - anchor_info['rlen']]
            
            compare_len = min(len(anchor_seq), len(ref_full))
            if compare_len == 0: continue

            # --- V GENE: LEFT ALIGN ---
            if self.gene_type == 'V':
                if anchor_seq[:compare_len] == ref_full[:compare_len]:
                    perfect_matches.append(cand)

            # --- J GENE: RIGHT ALIGN ---
            else:
                if anchor_seq[-compare_len:] == ref_full[-compare_len:]:
                    perfect_matches.append(cand)
        
        # Decision Logic & Error Reporting
        match_status = "Perfect"
        note = "OK"

        if perfect_matches:
            perfect_matches.sort(key=natural_sort_key)
            selected_allele = perfect_matches[0]
            # perfect_matches list preserved all candidates, used during log output
        else:
            # Fallback Logic
            raw_hit = str(best_row['AnalyzedGene'])
            
            # Scenario 1: Fallback to MiXCR call
            if raw_hit in candidates:
                selected_allele = raw_hit
                match_status = "Fallback_MiXCR"
                note = f"ALIGNMENT FAILED. Used MiXCR raw hit {raw_hit}."
            
            # Scenario 2: Force First Candidate (Severe Error)
            else:
                selected_allele = candidates[0]
                match_status = "Fallback_First"
                note = f"ALIGNMENT FAILED & MiXCR hit invalid. FORCED selection of first candidate {selected_allele}."

        # Log Selection with detailed status AND candidate list
        self.ref_selection_log.append({
            'Family': family_name,
            'Selected': selected_allele,
            'AllCandidates': perfect_matches, # New: store all matching items
            'ByCloneID': anchor_info['id'],
            'AnchorSeq': anchor_seq,
            'Status': match_status,
            'Note': note
        })
        
        return selected_allele

    def run_analysis(self):
        print(f"Running analysis...")
        results = []
        groups = self.df.groupby('AnalyzedFamily')
        
        for family, group in groups:
            if family == 'nan' or family not in self.ref_dict: continue
            
            # 1. Determine Reference
            candidates = sorted(list(self.ref_dict[family].keys()), key=natural_sort_key)
            ref_allele = self.resolve_reference(group, candidates, family)
            ref_seq_full = self.ref_dict[family][ref_allele] # V: includes Leader/J: full J region
            
            # Determine mature gene absolute start coordinate (0-based)
            mature_start_abs = 0
            if self.gene_type == 'V':
                if family in v_index:
                    mature_start_abs = v_index[family]
                else:
                    # If V gene not in v_index, skip or give warning (MiXCR results should include)
                    self.error_detail_logs.append(f"Family {family} missing from v_index. Skipping group.")
                    continue
            
            # Mature gene absolute end coordinate (based on ref_seq_full length)
            mature_end_abs = mature_start_abs + len(ref_seq_full)
            
            # 2. Process Clones
            for (align_str, full_target_seq), subgroup in group.groupby([self.cols['align'], 'targetSequences']):
                align_info = parse_alignment_string(align_str)
                if not align_info: continue
                
                if self.gene_type == 'J':
                    mature_end_abs = align_info['refLen']
                    mature_start_abs = mature_end_abs - len(ref_seq_full)

                alignment_status = "Valid"
                has_indel = False
                has_snp = False
                
                # --- A. Mutation Coordinate Conversion & Boundary Filtering ---
                # MiXCR mutStr reports coordinate P_mixcr based on 1-based, relative to full reference absolute.
                # parse_mutations internal already converted to 0-based, but still absolute.
                # MiXCR mutStr reports coordinate P_mixcr based on 1-based, relative to full reference absolute.
                # parse_mutations internal converts it to (P_mixcr - 1) and offset relative to align_info['refStart']
                # To simplify, we use target2ref back-deduction, temporarily not using offset (set to 0),
                # so that pos in raw_mutations is the offset relative to align_info['refStart'] (0-based).
                
                raw_mutations = parse_mutations(align_info['mutStr'], align_info['refStart'], align_info['refEnd'], 0) 
                valid_mutations = []

                gap_index = []
                # Convert coordinates, filter boundaries
                for m in raw_mutations:
                    # Calculate mutation absolute coordinate on full reference sequence (0-based)
                    p_abs = m['pos'] + align_info['refStart']

                    # Check if mutation is within mature gene range [mature_start_abs, mature_end_abs)
                    if p_abs < mature_start_abs or p_abs >= mature_end_abs:
                        self.bounds_logs.append(f"ID:{subgroup['cloneId'].iloc[0]} | IgnoredMut:{m['type']}@{p_abs} | MatureRefRange:[{mature_start_abs}:{mature_end_abs}]")
                        continue
                    
                    # Convert coordinate: relative coordinate from mature gene start (0-based)
                    m['pos'] = p_abs - mature_start_abs

                    if m['type'] == 'S' and m['target'] == 'N':
                        # Handle gap situation
                        gap_index.append(m['pos'])
                    else:
                        valid_mutations.append(m)

                gap_index.sort()

                # Record observed range
                observed_ranges = ""
                # Ensure observation is within range
                for r in align_info['ranges']:
                    ob_start = r[0]
                    ob_end = r[1]
                    if ob_start < mature_start_abs:
                        ob_start = mature_start_abs
                    if ob_end > mature_end_abs:
                        ob_end = mature_end_abs
                    ob_start -= mature_start_abs  # Calculate offset based on used reference frame
                    ob_end -= mature_start_abs
                    for g in gap_index:
                        if ob_start >= ob_end:
                            break
                        if g == ob_start:
                            ob_start = g + 1
                        if ob_start < g < ob_end:
                            observed_ranges += f"{ob_start},{g};"
                            ob_start = g + 1
                    if ob_start < ob_end:
                        observed_ranges += f"{ob_start},{ob_end};"
                observed_ranges = observed_ranges[: -1]

                # Check if there is mutation type
                has_indel = any(m['type'] in ['D', 'I'] for m in valid_mutations)
                has_snp = any(m['type'] == 'S' for m in valid_mutations)
                
                # --- B. Sequence Generation ---
                
                # Generate final sequence: use mature gene reference + filtered mutations
                final_output_seq = ref2target(ref_seq_full, valid_mutations)

                # --- C. Result Classification & Formatting ---
                
                category = "Germline"
                if alignment_status != "Valid":
                    category = "Error"
                elif has_indel:
                    category = "Indel_Variant"
                elif has_snp:
                    category = "PointMutation"
                
                self.stats['variant_types'][category] += len(subgroup)
                
                mut_out_list = []
                for m in valid_mutations:
                    p = m['pos'] # Already coordinate relative to mature gene
                    if m['type'] == 'S':
                        mut_out_list.append(f"S{m['ref']}{p}{m['target']}")
                    elif m['type'] == 'D':
                        mut_out_list.append(f"D{m['seq']}{p}")
                    elif m['type'] == 'I':
                        mut_out_list.append(f"I{p}{m['seq']}")
                
                mut_string_final = ";".join(mut_out_list) if mut_out_list else "Germline"
                
                pair_align_col_name = self.cols.get('pairAlign')

                div_feature_set = set()
                naive_feature_set = set()

                if 'aaSeqCDR3' in subgroup.columns and pair_align_col_name in subgroup.columns:
                    # 1. Preprocess column data (Vectorized operations for speed)
                    # Extract paired gene name
                    pair_genes_series = subgroup[self.cols['pair']].apply(lambda x: str(x).split('(')[0].split(',')[0])
                    # Extract CDR3 length
                    cdr3_lens_series = subgroup['aaSeqCDR3'].str.len().fillna(0).astype(int)
                    # Extract paired alignment string
                    pair_align_series = subgroup[pair_align_col_name]

                    # 2. Calculate standard diversity (Direct packing)
                    div_feature_set = set(zip(cdr3_lens_series, pair_genes_series))

                    # 3. Calculate Naive Diversity - Must check row by row
                    # Use zip to traverse three columns, avoiding overhead of apply
                    for c_len, p_gene, p_align in zip(cdr3_lens_series, pair_genes_series, pair_align_series):
                        
                        # Basic check
                        if pd.isna(p_align) or p_align == "":
                            continue
                        
                        # Parse alignment string
                        p_info = parse_alignment_string(p_align)

                        # Core decision logic:
                        # 1. p_info successfully parsed (Not None)
                        # 2. mutStr is empty (Indicates no mutation, i.e., Naive)
                        if p_info and (not p_info['mutStr']):
                            naive_feature_set.add((c_len, p_gene))

                results.append({
                    'GeneType': self.gene_type,
                    'Family': family,
                    'ReferenceAllele': ref_allele,
                    'Sequence': final_output_seq,
                    'MutationInfo': mut_string_final,
                    'CloneCount': len(subgroup),
                    'CloneCountSplit': len(subgroup),
                    'DiversityFeatures': div_feature_set,
                    'DiversityIndexSplit': len(div_feature_set),
                    'NaiveDiversityFeatures': naive_feature_set,
                    'NaiveDiversityIndexSplit': len(naive_feature_set),
                    'Status': category,
                    'Validation': alignment_status,
                    'ObservedRanges': observed_ranges
                })
                
        # --- AGGREGATION STEP --- (Unchanged)
        # ... (Aggregation logic)
        raw_df = pd.DataFrame(results)
        if raw_df.empty: return raw_df
        
        group_keys = ['GeneType', 'Family', 'ReferenceAllele', 'Sequence', 'MutationInfo', 'Status', 'Validation']
        agg_rules = {
            'CloneCount': 'sum',
            'DiversityFeatures': lambda x: set().union(*x),
            'NaiveDiversityFeatures': lambda x: set().union(*x),
            'ObservedRanges': lambda x: '|'.join(x),
            'CloneCountSplit': lambda x: '|'.join(map(str, x)),
            'DiversityIndexSplit': lambda x: '|'.join(map(str, x)),
            'NaiveDiversityIndexSplit': lambda x: '|'.join(map(str, x)),
        }
        aggregated_df = raw_df.groupby(group_keys, as_index=False).agg(agg_rules)
        aggregated_df['DiversityIndex'] = aggregated_df['DiversityFeatures'].apply(len)
        aggregated_df['NaiveDiversityIndex'] = aggregated_df['NaiveDiversityFeatures'].apply(len)
        aggregated_df = aggregated_df.drop(columns=['DiversityFeatures', 'NaiveDiversityFeatures'])
        
        return aggregated_df

    def write_outputs(self, result_df, prefix):
        # 1. Sequence CSV
        csv_file = f"{prefix}_sequences.csv"
        if not result_df.empty:
            out_df = result_df[result_df['Status'] != 'Error'].copy()
            out_df['Family_Sort'] = out_df['Family'].apply(natural_sort_key)
            out_df = out_df.sort_values(by=['Family_Sort', 'DiversityIndex'], ascending=[True, False])
            
            cols = ['GeneType', 'Family', 'ReferenceAllele', 'Sequence', 'MutationInfo', 'CloneCount', 'DiversityIndex', 'NaiveDiversityIndex', 'Status', 'ObservedRanges', 'CloneCountSplit', 'DiversityIndexSplit', 'NaiveDiversityIndexSplit']
            out_df[cols].to_csv(csv_file, index=False)
            print(f"Sequence CSV saved: {csv_file}")

        # 2. Comprehensive Report
        report_file = f"{prefix}_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== TCR Variant Analysis Expert Report ===\n")
            
            # A. Filtering Stats
            f.write("\n1. Filtering Statistics\n")
            init = self.stats['total_input']
            f.write(f"   Total Input Clones: {init}\n")
            for k, v in self.stats['filters'].items():
                pct = (v / init * 100) if init > 0 else 0
                f.write(f"   - Filtered by {k}: {v} ({pct:.2f}%)\n")
            
            final = self.stats['final_count']
            final_pct = (final / init * 100) if init > 0 else 0
            f.write(f"   => Final Valid Entries (Post-Split): {final} ({final_pct:.2f}% of input)\n")

            # B. Ref Selection (Updated Formatting)
            f.write("\n2. Reference Allele Selection Log:\n")
            f.write("   (Logic: Longest clean alignment coverage serves as anchor. V=Left, J=Right)\n")
            ref_log_df = pd.DataFrame(self.ref_selection_log)
            
            if not ref_log_df.empty:
                for fam, sub in ref_log_df.groupby('Family'):
                    # Extract info
                    sel = sub['Selected'].iloc[0]
                    evidence = sub['ByCloneID'].iloc[0]
                    seq = sub['AnchorSeq'].iloc[0]
                    status = sub['Status'].iloc[0]
                    note = sub['Note'].iloc[0]
                    candidates = sub['AllCandidates'].iloc[0]
                    
                    # Improved format 1: handle multiple candidate display
                    candidates_str = ""
                    if isinstance(candidates, list) and len(candidates) > 1:
                        # Exclude the selected one (candidates[0]), display others
                        others = ", ".join(candidates[1:])
                        candidates_str = f", Others including: {others}"

                    # Highlight Fallbacks
                    if status == "Perfect":
                        # Output: - TRBV5-6: Selected TRBV5-6*01, Others including: TRBV5-6*02
                        f.write(f"   - {fam}: Selected {sel}{candidates_str}\n")
                    else:
                        f.write(f"   [CRITICAL WARNING] {fam}: Selected {sel} ({status})\n")
                        f.write(f"      -> Reason: {note}\n")
                    
                    # Improved format 2: ID moved after Sequence
                    f.write(f"      Used Sequence: {seq} (Clone ID: {evidence})\n")
            
            # C. Variant Stats
            f.write("\n3. Variant Classification Statistics:\n")
            total_processed = sum(self.stats['variant_types'].values())
            f.write(f"   Total Clones Analyzed: {total_processed}\n")
            for k, v in self.stats['variant_types'].items():
                pct = (v / total_processed * 100) if total_processed > 0 else 0
                f.write(f"   - {k}: {v} ({pct:.2f}%)\n")
            
            f.write(f"   - Mismatches/Errors (excluded from output): {self.stats['mismatch_counts']}\n")
            
            # D. Warning Logs
            f.write("\n4. Long Read Truncation Logs (Clone > Ref):\n")
            if self.long_read_logs:
                for l in self.long_read_logs: # No limit
                    f.write(f"   {l}\n")
            else:
                f.write("   None.\n")
            
            f.write("\n5. Out-of-Bounds Mutation Logs (Ignored):\n")
            if self.bounds_logs:
                for l in self.bounds_logs: # No limit
                    f.write(f"   {l}\n")
            else:
                f.write("   None.\n")

            f.write("\n6. Error/Mismatch Detail Logs (Excluded from Output):\n")
            if self.error_detail_logs:
                for l in self.error_detail_logs:
                    f.write(f"   {l}\n")
            else:
                f.write("   None.\n")

        print(f"Report saved: {report_file}")

# ==========================================
# 3. Entry Point
# ==========================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TCR Variant Analysis Expert V2.3")
    parser.add_argument("--input", required=True)
    parser.add_argument("--gene", required=True, choices=['V', 'J', 'v', 'j'])
    parser.add_argument("--ref", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--minClns", default=10)
    parser.add_argument("--minQual", default=0, type=int, help="Minimum quality score threshold (default: 0)")
    parser.add_argument('--min_count', type=int, default=1, help='Minimum count threshold (Default: 1)')
    parser.add_argument('--filter_productive', type=str2bool, default=False, help='Filter non-productive sequences (Default: True)')
    
    args = parser.parse_args()
    
    analyzer = VariantAnalyzer(args.input, args.gene, args.ref, int(args.minClns), args.minQual, args.min_count, args.filter_productive)
    analyzer.preprocess_and_filter()
    res_df = analyzer.run_analysis()
    analyzer.write_outputs(res_df, args.prefix)