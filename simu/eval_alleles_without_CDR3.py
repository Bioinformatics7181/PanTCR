import pandas as pd
import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Evaluate inferred alleles against genotype ground truth.")
    parser.add_argument('--gt', required=True, help='Genotype CSV (Ground Truth)')
    parser.add_argument('--infer', required=True, help='Inference CSV (Predictions)')
    parser.add_argument('--pmtr', required=True, help='pmTR split CSV (for trimmed sequences)')
    parser.add_argument('--index', required=True, help='TRB Index CSV (for V/J coordinates)')
    parser.add_argument('--gene_type', required=True, choices=['V', 'J'], help='Gene type to evaluate (V or J)')
    parser.add_argument('--intersect', action='store_true', help='If set, only evaluate genes present in BOTH ground truth and predictions.')
    parser.add_argument('--out_prefix', default='report', help='Output prefix for reports')
    return parser.parse_args()

def sequences_match(x_seq, y_seq, gene_type):
    """
    Compare two sequences according to gene type.
    - V: Left-aligned (compare prefix)
    - J: Right-aligned (compare suffix)
    Returns True if they match (exactly for the overlapping region).
    """
    x_seq = x_seq.strip().upper()
    y_seq = y_seq.strip().upper()
    compare_len = min(len(x_seq), len(y_seq))
    
    if compare_len == 0:
        return False
    
    if gene_type == 'V':
        return x_seq[:compare_len] == y_seq[:compare_len]
    elif gene_type == 'J':
        return x_seq[-compare_len:] == y_seq[-compare_len:]
    else:
        return False

def main():
    # 1. Parse arguments
    try:
        args = parse_args()
    except:
        print("Warning: Running without command line arguments (Interactive Mode)")
        class Args:
            gt = 'genotype_AFR_seed0.csv'
            infer = 'infer_AFR_seed0.csv'
            pmtr = 'pmTR_TRB_V_J_split_2_label.csv'
            index = 'TRB_index.csv'
            gene_type = 'J' 
            intersect = False
            out_prefix = 'report_J_test'
        args = Args()

    print(f"Loading data for gene type: {args.gene_type}...")
    if args.intersect:
        print("Mode: Intersection Only (evaluating common genes only)")
    else:
        print("Mode: Standard (evaluating all genes in Ground Truth)")

    # 2. Load Data
    genotype_df = pd.read_csv(args.gt)
    infer_df = pd.read_csv(args.infer)
    pmTR_df = pd.read_csv(args.pmtr)
    trb_index_df = pd.read_csv(args.index)

    # 3. Filter Data (Based on gene_type)
    filter_str = f"TRB{args.gene_type}"
    
    if 'gene' not in genotype_df.columns or 'gene' not in infer_df.columns:
        print("Error: 'gene' column missing in inputs.")
        sys.exit(1)

    # Initial filtering for V or J rows
    genotype_filtered = genotype_df[genotype_df['gene'].astype(str).str.contains(filter_str)].copy()
    infer_filtered = infer_df[infer_df['gene'].astype(str).str.contains(filter_str)].copy()

    # 4. (Optional) Intersection Mode: Filter for common genes
    if args.intersect:
        gt_genes = set(genotype_filtered['gene'])
        infer_genes = set(infer_filtered['gene'])
        common_genes = gt_genes.intersection(infer_genes)
        
        print(f"  - Genes in GT: {len(gt_genes)}")
        print(f"  - Genes in Infer: {len(infer_genes)}")
        print(f"  - Common Genes: {len(common_genes)}")
        
        if len(common_genes) == 0:
            print("Warning: No common genes found between Ground Truth and Prediction. Stopping.")
            return

        genotype_filtered = genotype_filtered[genotype_filtered['gene'].isin(common_genes)]
        infer_filtered = infer_filtered[infer_filtered['gene'].isin(common_genes)]

    print(f"  - Ground Truth rows (final): {len(genotype_filtered)}")
    print(f"  - Inference rows (final): {len(infer_filtered)}")

    if len(genotype_filtered) == 0:
        print(f"No ground truth data remains after filtering. Exiting.")
        return

    # ==============================================================================
    # 5. Build Set X (Ground Truth)
    # ==============================================================================
    print("Building Set X (Ground Truth)...")
    # Assumption: pmTR_df['sequence_trimmed'] is already correctly trimmed for the specific gene type
    # (e.g., CDR3 removed for V, or only FR4 kept for J)
    allele_to_seq = dict(zip(pmTR_df['allele'], pmTR_df['sequence_trimmed']))

    X_dict = {}  # sequence -> set of gene types (alleles)

    for _, row in genotype_filtered.iterrows():
        for col in ['allele_A', 'allele_B']:
            if pd.isna(row[col]): continue
            
            a_name = row[col]
            if a_name in allele_to_seq:
                seq = str(allele_to_seq[a_name]).strip().upper()
                if seq:
                    if seq not in X_dict:
                        X_dict[seq] = set()
                    X_dict[seq].add(a_name)

    X_seqs = list(X_dict.keys())
    print(f"  - Total unique sequences in X: {len(X_seqs)}")

    # ==============================================================================
    # 6. Build Set Y (Predictions)
    # ==============================================================================
    print("Building Set Y (Predictions)...")
    
    # Convert columns to numeric, handling errors
    cols_to_numeric = ['FR1Begin', 'CDR3Begin', 'FR4Begin']
    for col in cols_to_numeric:
        if col in trb_index_df.columns:
            trb_index_df[col] = pd.to_numeric(trb_index_df[col], errors='coerce')
    
    index_map = {}
    for _, row in trb_index_df.iterrows():
        gene = row['baseGene']
        
        # Store coordinates based on Gene Type
        if args.gene_type == 'V':
            if pd.notnull(row['FR1Begin']) and pd.notnull(row['CDR3Begin']):
                # V gene: Keep region between FR1Begin and CDR3Begin
                index_map[gene] = {
                    'start': int(row['FR1Begin']),
                    'end': int(row['CDR3Begin'])
                }
        elif args.gene_type == 'J':
            # J gene: Keep region from FR4Begin onwards
            start_val = 0 # Default to 0 if missing
            
            if 'FR4Begin' in row and pd.notnull(row['FR4Begin']):
                start_val = int(row['FR4Begin'])
            
            # NOTE: We do not check CDR3End anymore. 
            # If FR4Begin is missing, start_val remains 0.
            
            index_map[gene] = {
                'start': start_val,
                'end': None # J takes everything until the end
            }

    Y_dict = {}  # sequence -> set of gene types (genes)

    skipped_genes = 0
    for _, row in infer_filtered.iterrows():
        gene = row['gene']
        
        if gene in index_map:
            coords = index_map[gene]
            
            for col in ['seq_A', 'seq_B']:
                if pd.isna(row[col]): continue
                
                full_seq = str(row[col])
                trimmed_seq = ""

                # === Core Trimming Logic ===
                if args.gene_type == 'V':
                    # V Logic: Slice [0 : CDR3Begin-FR1Begin]
                    # (Assuming full_seq starts relative to Index 0)
                    cut_point = coords['end'] - coords['start']
                    if len(full_seq) >= cut_point:
                        trimmed_seq = full_seq[:cut_point].strip().upper()
                
                elif args.gene_type == 'J':
                    # J Logic: Slice [FR4Begin : ]
                    # Discard the CDR3 part at the beginning, keep FR4
                    start_point = coords['start']
                    if len(full_seq) > start_point:
                        trimmed_seq = full_seq[start_point:].strip().upper()

                if trimmed_seq:
                    if trimmed_seq not in Y_dict:
                        Y_dict[trimmed_seq] = set()
                    Y_dict[trimmed_seq].add(gene)
        else:
            skipped_genes += 1

    if skipped_genes > 0:
        print(f"  - Warning: {skipped_genes} rows skipped due to missing gene index info.")

    Y_seqs = list(Y_dict.keys())
    print(f"  - Total unique sequences in Y: {len(Y_seqs)}")

    # ==============================================================================
    # 7. Validation
    # ==============================================================================
    print("Validating...")
    TP_list = []
    FP_list = []
    FN_list = []
    matched_X_indices = set()

    for y_seq in Y_seqs:
        is_match = False
        y_genes = ",".join(sorted(list(Y_dict[y_seq])))
        
        for i, x_seq in enumerate(X_seqs):
            # Comparison Logic:
            # For V: X and Y ending at CDR3Begin without CDR3.
            # For J: X and Y starting at FR4Begin without CDR3.
            
            if sequences_match(x_seq, y_seq, args.gene_type):
                is_match = True
                matched_X_indices.add(i)
        
        if is_match:
            TP_list.append({'Sequence': y_seq, 'GeneType': y_genes, 'Status': 'TP'})
        else:
            FP_list.append({'Sequence': y_seq, 'GeneType': y_genes, 'Status': 'FP'})

    # Identify False Negatives
    for i, x_seq in enumerate(X_seqs):
        if i not in matched_X_indices:
            x_genes = ",".join(sorted(list(X_dict[x_seq])))
            FN_list.append({'Sequence': x_seq, 'GeneType': x_genes, 'Status': 'FN'})

    tp_count = len(TP_list)
    recall_count = len(matched_X_indices)
    
    precision = tp_count / len(Y_seqs) if len(Y_seqs) > 0 else 0
    recall = recall_count / len(X_seqs) if len(X_seqs) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"\nResults ({args.gene_type}, Intersect={args.intersect}):")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1-score:  {f1:.4f}")

    # ==============================================================================
    # 8. Generate Reports
    # ==============================================================================
    summary_df = pd.DataFrame([
        {'Metric': 'GeneType', 'Value': args.gene_type},
        {'Metric': 'IntersectOnly', 'Value': args.intersect},
        {'Metric': 'Precision', 'Value': precision},
        {'Metric': 'Recall', 'Value': recall},
        {'Metric': 'F1-score', 'Value': f1}
    ])
    summary_file = f"{args.out_prefix}_summary.csv"
    summary_df.to_csv(summary_file, index=False)

    detailed_df = pd.DataFrame(TP_list + FP_list + FN_list)
    if not detailed_df.empty:
        detailed_df = detailed_df[['Sequence', 'GeneType', 'Status']]
    else:
        detailed_df = pd.DataFrame(columns=['Sequence', 'GeneType', 'Status'])
        
    detailed_file = f"{args.out_prefix}_detailed.csv"
    detailed_df.to_csv(detailed_file, index=False)
    
    print(f"Reports saved to {summary_file} and {detailed_file}")

if __name__ == "__main__":
    main()