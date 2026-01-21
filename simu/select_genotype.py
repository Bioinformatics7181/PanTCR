#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import os
import sys

def parse_gene_name(allele_name):
    """
    Extract gene name from allele name.
    e.g.: 'TRBV5-1*01' -> 'TRBV5-1'
    """
    if pd.isna(allele_name): return "Unknown"
    if '*' in allele_name:
        return allele_name.split('*')[0]
    return allele_name

def str_to_bool(val):
    """
    Transfer the str value (TRUE/FALSE) to Python bool
    """
    if isinstance(val, bool):
        return val
    s = str(val).strip().upper()
    return s in ("T", "TRUE")

def select_alleles(input_file, population, seed, output_dir, prefix, functionality, exdefault_val):
    # 1. Smart Load (CSV vs TSV)
    print(f"Loading database from {input_file}...")
    
    try:
        if input_file.lower().endswith('.tsv') or input_file.lower().endswith('.txt'):
            df = pd.read_csv(input_file, sep='\t')
        else:
            # Default to CSV
            df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # 2. Check basic columns
    required_cols = ['allele', 'sequence']
    for col in required_cols:
        if col not in df.columns:
            # Handle case sensitivity (e.g., Sequence vs sequence)
            found = False
            for existing_col in df.columns:
                if existing_col.lower() == col.lower():
                    df.rename(columns={existing_col: col}, inplace=True)
                    found = True
                    break
            if not found:
                print(f"Error: Required column '{col}' missing from input file.")
                sys.exit(1)

    # --- Feature: Filter by Functionality ---
    if str_to_bool(functionality):
        # Find column (case insensitive)
        func_col = None
        for col in df.columns:
            if col.lower() == 'functionality':
                func_col = col
                break
        
        if func_col:
            print(f"Filtering alleles with {func_col} = '{functionality}'...")
            initial_count = len(df)
            # Force string comparison to avoid type mismatch
            df = df[df[func_col].astype(str) == str(functionality)]
            print(f"  Kept {len(df)} out of {initial_count} alleles.")
            
            if len(df) == 0:
                print("Warning: No alleles left after functionality filtering!")
        else:
            print(f"Warning: --functionality specified as '{functionality}', but column 'functionality' not found. Filtering skipped.")

    # --- Feature: Filter by is_default ---
    if str_to_bool(exdefault_val):
        # Find column (case insensitive, usually 'is_default')
        default_col = None
        for col in df.columns:
            if col.lower() == 'is_default':
                default_col = col
                break
        
        if default_col:
            print(f"Applying priority filtering based on {default_col} corresponding to '{exdefault_val}'...")

            if 'gene' not in df.columns:
                df['gene'] = df['allele'].apply(parse_gene_name)
                
            temp_bool_col = '__is_default_normalized__'
            df[temp_bool_col] = df[default_col].apply(str_to_bool)
            
            filtered_dfs = []
            
            for gene, group in df.groupby('gene'):
                has_special_allele = (~group[temp_bool_col]).any()
                
                if has_special_allele:
                    kept_rows = group[group[temp_bool_col] == False]
                else:
                    kept_rows = group
                    
                filtered_dfs.append(kept_rows)
                
            if filtered_dfs:
                df = pd.concat(filtered_dfs, ignore_index=True)
                df.drop(columns=[temp_bool_col], inplace=True)
                print(f"After is_default filtering: {len(df)} alleles remaining.")
            else:
                df = pd.DataFrame(columns=df.columns)
        else:
            print(f"Warning: --exdefault specified as '{exdefault_val}', but column 'is_default' not found. Filtering skipped.")

    # 3. Determine Sampling Strategy (Frequency vs Uniform)
    use_uniform = False
    
    if not population:
        print("No population specified (-p). Using Uniform Sampling (equal probability).")
        use_uniform = True
        pop_label = "UNI"
    elif population not in df.columns:
        print(f"Population column '{population}' not found in file. Fallback to Uniform Sampling.")
        use_uniform = True
        pop_label = "UNI"
    else:
        print(f"Using frequency data from population: {population}")
        use_uniform = False
        pop_label = population

    # 4. Preprocessing: Extract Gene Name
    if 'gene' not in df.columns:
        df['gene'] = df['allele'].apply(parse_gene_name)
    
    # 5. Set Random Seed
    np.random.seed(seed)
    print(f"Simulating genotype... Seed: {seed}")
    
    genotype_rows = []
    
    # 6. Group by Gene and Sample Diploid
    grouped = df.groupby('gene')
    
    total_genes = 0
    
    for gene, group in grouped:
        alleles = group['allele'].values
        sequences = group['sequence'].values
        
        # Calculate probabilities
        if use_uniform:
            # Uniform sampling with replacement
            probs = np.ones(len(alleles)) / len(alleles)
        else:
            # Frequency-based sampling
            counts = group[population].values
            # Handle non-numeric or NaN
            counts = pd.to_numeric(counts, errors='coerce')
            counts = np.nan_to_num(counts, nan=0.0)
            
            total_count = np.sum(counts)
            
            if total_count > 0:
                probs = counts / total_count
            else:
                # If all counts are 0, fallback to uniform for this gene
                probs = np.ones(len(alleles)) / len(alleles)
            
        # Sample (Size=2, Replace=True for diploid)
        indices = np.arange(len(alleles))
        
        # Handle case where filtering left no alleles for a specific gene (unlikely but possible)
        if len(indices) == 0:
            continue

        chosen_indices = np.random.choice(indices, size=2, replace=True, p=probs)
        
        idx_A = chosen_indices[0]
        idx_B = chosen_indices[1]
        
        # Construct output row
        row = {
            'gene': gene,
            'allele_A': alleles[idx_A],
            'seq_A': sequences[idx_A],
            'allele_B': alleles[idx_B],
            'seq_B': sequences[idx_B]
        }
        genotype_rows.append(row)
        total_genes += 1
        
    # 7. Generate Result DataFrame
    result_df = pd.DataFrame(genotype_rows)
    
    if result_df.empty:
        print("Error: No genotype generated. Please check your filters.")
        return

    # Construct output filename
    # Incorporate exdefault info into filename if used, to avoid overwriting
    exdefault_suffix = ""
    if str_to_bool(exdefault_val):
        exdefault_suffix = f"_exdefaultT"
        
    output_filename = f"{prefix}_{pop_label}{exdefault_suffix}_seed{seed}.csv"
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, output_filename)
    
    # Save
    result_df.to_csv(output_path, index=False)
    print(f"Successfully generated genotype file: {output_path}")
    print(f"Total genes processed: {total_genes}")

def main():
    parser = argparse.ArgumentParser(description="Simulate a diploid genotype (Supports CSV/TSV, Frequency/Uniform).")
    
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the allele file (CSV or TSV). Must contain 'allele' and 'sequence' columns."
    )
    parser.add_argument(
        "-p", "--pop", 
        default=None, 
        help="Target population column (e.g., EAS). If omitted or not found, uses uniform sampling."
    )
    parser.add_argument(
        "-s", "--seed", 
        default=42, 
        type=int, 
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "-o", "--output_dir", 
        default=".", 
        help="Directory to save the output CSV"
    )
    parser.add_argument(
        "--prefix", 
        default="genotype", 
        help="Prefix for the output filename (default: 'genotype')."
    )

    parser.add_argument(
        "-f", "--functionality", 
        default="", 
        help="Filter sequences by functionality column (e.g., 'F', 'ORF'). Default is '' (no filter)."
    )

    # New Argument
    parser.add_argument(
        "-d", "--exdefault", 
        default="", 
        help="Filter by 'is_default' column. Accepts 'T', 'True' or 'F', 'False'. Default is '' (no filter)."
    )
    
    args = parser.parse_args()
    
    select_alleles(args.input, args.pop, args.seed, args.output_dir, args.prefix, args.functionality, args.exdefault)

if __name__ == "__main__":
    main()
