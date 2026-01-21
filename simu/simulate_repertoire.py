#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import pickle
import argparse
import random
import re
import sys

# ================= 1. Sequence Processing Helper Functions =================

def get_v_sequence_before_cdr3(full_seq):
    """
    Truncates V gene sequence, keeping the part before the last Cys.
    Strategy: Find the last Cys (TGT/TGC), keep all sequences before its start.
    """
    matches = list(re.finditer(r'(TGT|TGC)', full_seq))
    if not matches:
        return full_seq # Defensive programming
    return full_seq[:matches[-1].start()]

def smart_stitch_j(cdr3_seq, j_germline_seq):
    """
    Calculates the J gene fragment needed to be stitched after CDR3.
    
    Logic:
    1. Find the start of the FGXG Motif in J Germline.
    2. Use the fragment starting from the Motif as the 'candidate fragment'.
    3. Check if the end of CDR3 already contains part of the 'candidate' (Overlap).
    4. Return: Candidate fragment - Overlapping part.
    """
    if not cdr3_seq or not j_germline_seq:
        return ""

    # 1. Define Motif Regex (covers F/W - G - X - G)
    # (TTC|TTT|TGG) : Phe/Trp
    # GG.           : Gly
    # [ATGC]{3}     : Any amino acid
    # GG.           : Gly
    motif_pattern = r'(TTC|TTT|TGG)GG.[ATGC]{3}GG.'
    
    # 2. Search for Motif in J Germline
    j_match = re.search(motif_pattern, j_germline_seq)
    
    if not j_match:
        # If Motif is not found (extremely rare), conservative strategy:
        # Assume CDR3 hasn't covered FR4 yet, return whole sequence or half
        # To prevent repetition, return empty or last 30bp
        # Depends on simulation aggressiveness. Here return last 30bp as fallback.
        return j_germline_seq[-30:] if len(j_germline_seq) > 30 else j_germline_seq

    # Get the start index of the Motif in the J gene
    j_motif_start_index = j_match.start()
    
    # 'Candidate fragment': J gene sequence starting from the Motif (FR4 + 3'UTR)
    j_candidate = j_germline_seq[j_motif_start_index:]
    
    # 3. Calculate overlap between CDR3 tail and J candidate head
    overlap_len = 0
    
    # Only check overlap within Motif length range (e.g., max 15bp)
    # This step prevents accidental matching deep inside CDR3
    max_check = min(len(cdr3_seq), 15) 
    
    for i in range(max_check, 0, -1):
        # Slice the last i bases of CDR3
        cdr3_tail = cdr3_seq[-i:]
        # Slice the first i bases of the J candidate
        j_head = j_candidate[:i]
        
        if cdr3_tail == j_head:
            overlap_len = i
            break # Found maximum match length, stop immediately
            
    # 4. Return the part to be appended
    # I.e.: J candidate fragment - the first i bases that overlapped
    j_part_to_append = j_candidate[overlap_len:]
    
    return j_part_to_append

# ================= 2. Gene Family Parsing =================

def get_gene_family(gene_name):
    """ TRBV5-1 -> TRBV5; TRBJ2-7 -> TRBJ2 """
    if pd.isna(gene_name): return ""
    base = gene_name.split('*')[0]
    if '-' in base:
        return base.split('-')[0]
    return base

# ================= 3. 5-Tier Fallback Query Logic =================

def get_cdr3_tiered_v_family(target_v, target_j, db):
    """
    Returns: (cdr3_sequence, method_name)
    """
    lookup = db['lookup']
    
    # --- Level A: Target J exists ---
    if target_j in lookup:
        j_group = lookup[target_j]
        
        # Tier 1: Exact Match
        if target_v in j_group and j_group[target_v]:
            return random.choice(j_group[target_v]), "Tier1_Exact"
        
        # Tier 2: V-Family Fallback
        target_v_fam = get_gene_family(target_v)
        same_fam_vs = [v for v in j_group.keys() if get_gene_family(v) == target_v_fam]
        
        if same_fam_vs:
            fallback_v = random.choice(same_fam_vs)
            return random.choice(j_group[fallback_v]), "Tier2_V_Family"
            
        # Tier 3: J-Fallback
        available_vs = list(j_group.keys())
        if available_vs:
            fallback_v = random.choice(available_vs)
            return random.choice(j_group[fallback_v]), "Tier3_J_Fallback"

    # --- Level B: Target J missing ---
    
    # Tier 4: J-Family Fallback
    target_j_fam = get_gene_family(target_j)
    same_fam_js = [j for j in lookup.keys() if get_gene_family(j) == target_j_fam]
    
    if same_fam_js:
        fallback_j = random.choice(same_fam_js)
        fallback_v = random.choice(list(lookup[fallback_j].keys()))
        return random.choice(lookup[fallback_j][fallback_v]), "Tier4_J_Family"
        
    # Tier 5: Global Fallback
    if not lookup: return None, "Fail"
    random_j = random.choice(list(lookup.keys()))
    random_v = random.choice(list(lookup[random_j].keys()))
    return random.choice(lookup[random_j][random_v]), "Tier5_Global"

# ================= 4. Main Simulation Loop =================

def simulate(genotype_file, dict_file, n_clones, n_reads, alpha, output_prefix, seed):
    # --- Set Random Seed (Core Modification) ---
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
        print(f"Random Seed set to: {seed}")
    else:
        print("No Random Seed set. Results will be stochastic.")

    # 1. Load Resources
    print(f"Loading Genotype: {genotype_file}")
    geno_df = pd.read_csv(genotype_file)
    
    print(f"Loading Dictionary: {dict_file}")
    with open(dict_file, 'rb') as f:
        cdr3_db = pickle.load(f)

    # Determine available gene pool
    my_v_genes = geno_df[geno_df['gene'].str.startswith('TRBV')]['gene'].unique()
    my_j_genes = geno_df[geno_df['gene'].str.startswith('TRBJ')]['gene'].unique()
    
    if len(my_v_genes) == 0:
        sys.exit("Error: No TRBV genes found in genotype file.")

    repertoire_data = []
    unique_check = set()
    
    print(f"Generating {n_clones} unique clones...")
    
    count = 0
    attempts = 0
    max_attempts = n_clones * 200
    
    while count < n_clones:
        attempts += 1
        if attempts > max_attempts:
            print("Warning: Max attempts reached.")
            break
            
        # A. Randomly pair V and J
        v_gene = np.random.choice(my_v_genes)
        j_gene = np.random.choice(my_j_genes)
        
        # B. Get CDR3
        selected_cdr3, method = get_cdr3_tiered_v_family(v_gene, j_gene, cdr3_db)
        
        if not selected_cdr3: continue
            
        # C. Online De-duplication
        if (v_gene, j_gene, selected_cdr3) in unique_check: continue
            
        # D. Instantiation (Allele Selection)
        v_rows = geno_df[geno_df['gene'] == v_gene]
        j_rows = geno_df[geno_df['gene'] == j_gene]
        
        chrom_v = random.choice(['A', 'B'])
        v_allele = v_rows.iloc[0][f'allele_{chrom_v}']
        v_seq = v_rows.iloc[0][f'seq_{chrom_v}']
        
        chrom_j = random.choice(['A', 'B'])
        j_allele = j_rows.iloc[0][f'allele_{chrom_j}']
        j_seq = j_rows.iloc[0][f'seq_{chrom_j}']
        
        # E. Smart Stitching
        v_part = get_v_sequence_before_cdr3(v_seq)
        j_part = smart_stitch_j(selected_cdr3, j_seq)
        
        full_seq = v_part + selected_cdr3 + j_part
        
        # F. Save
        unique_check.add((v_gene, j_gene, selected_cdr3))
        
        repertoire_data.append({
            'clone_id': f"clone_{count+1:06d}",
            'v_gene': v_gene,
            'j_gene': j_gene,
            'v_allele': v_allele,
            'j_allele': j_allele,
            'cdr3_nt': selected_cdr3,
            'full_sequence': full_seq,
            'method': method
        })
        count += 1
        
    print(f"Successfully generated {count} clones.")
    
    # 2. Statistics
    rep_df = pd.DataFrame(repertoire_data)
    print("\n[Method Statistics]")
    print(rep_df['method'].value_counts())
    
    # 3. Abundance Simulation
    # Random shuffle (Controlled by np.random.seed)
    rep_df = rep_df.sample(frac=1, random_state=seed).reset_index(drop=True)
    
    ranks = np.arange(1, len(rep_df) + 1)
    weights = 1 / np.power(ranks, alpha)
    probs = weights / np.sum(weights)
    
    # Simulate reads (Controlled by np.random.seed)
    rep_df['frequency'] = probs
    rep_df['read_count'] = np.random.multinomial(n_reads, probs)
    
    # 4. Output
    csv_out = f"{output_prefix}_repertoire.csv"
    rep_df.to_csv(csv_out, index=False)
    
    fasta_out = f"{output_prefix}_transcripts.fasta"
    with open(fasta_out, 'w') as f:
        for _, row in rep_df.iterrows():
            if row['read_count'] > 0:
                header = (f">{row['clone_id']}|{row['v_allele']}|{row['j_allele']}|"
                          f"reads={row['read_count']}|method={row['method']}")
                f.write(f"{header}\n{row['full_sequence']}\n")
                
    print(f"\nDone. Output files:\n  CSV: {csv_out}\n  FASTA: {fasta_out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate TCR repertoire with Hybrid Strategy")
    
    parser.add_argument("--genotype", required=True, help="Simulated Genotype CSV")
    parser.add_argument("--dict", required=True, help="CDR3 Dictionary PKL")
    parser.add_argument("-nc", "--num_clones", type=int, default=10000, help="Number of clones to generate")
    parser.add_argument("-nr", "--num_reads", type=int, default=100000, help="Number of reads to generate")
    parser.add_argument("--alpha", type=float, default=1.0, help="Zipfian distribution alpha")
    parser.add_argument("-o", "--output", required=True, help="Output file prefix")
    
    # Added Seed Parameter
    parser.add_argument("-s", "--seed", type=int, default=None, help="Random seed (int) for reproducibility")
    
    args = parser.parse_args()
    
    simulate(args.genotype, args.dict, args.num_clones, args.num_reads, args.alpha, args.output, args.seed)