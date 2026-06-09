import argparse
import pandas as pd
import re
import sys

def parse_args():
    """Configure and parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Extract gene names from TSV, find sequences in IMGT ref, and output to specific CSV."
    )
    parser.add_argument(
        '--clones', 
        required=True, 
        help="Input file path: MiXCR generated clones file (.tsv)"
    )
    parser.add_argument(
        '--ref', 
        required=True, 
        help="Reference database file path: containing allele sequences (.csv)"
    )
    parser.add_argument(
        '--gene',
        choices=['V', 'J', 'all'],
        default='all',
        help="Target gene segment type (V, J or all) to read from corresponding allXHitsWithScore columns."
    )
    parser.add_argument(
        '--output', 
        required=True, 
        help="Output file path: CSV file to save results"
    )
    return parser.parse_args()

def extract_genes_from_tsv(df, target_segment='all'):
    """Extract unique gene names from TSV allVHitsWithScore columns (based on *00 format)"""
    unique_genes = set()

    col_map = {
        'V': 'allVHitsWithScore',
        'J': 'allJHitsWithScore'
    }
    
    # Check if columns exist
    columns_to_check = []
    if target_segment == 'all':
        columns_to_check = list(col_map.values())
    elif target_segment in col_map:
        columns_to_check = [col_map[target_segment]]

    found_cols = []
    
    for col in columns_to_check:
        # Check if column exists
        if col not in df.columns:
            continue
        found_cols.append(col)

        for hit_str in df[col].dropna():
            # Example hit_str: "TRBV1*00(1234),TRBV2*00(500)"
            hits = hit_str.split(',')
            for h in hits:
                # Remove score in brackets to get allele name
                allele_name = h.split('(')[0]
                # Only process entries ending with *00
                if allele_name.endswith('*00'):
                    # Extract gene name (remove *00)
                    gene = allele_name.split('*')[0]
                    unique_genes.add(gene)
    
    if not found_cols:
        print(f"Warning: No columns related to '{target_segment}' found in TSV (Expected: {columns_to_check}).")
        # Decide whether to sys.exit(1) if columns are mandatory
        # Here we only warn and allow processing other existing columns if 'all'
    return unique_genes

def build_reference_lookup(df_ref):
    """Build reference lookup dictionary: Gene -> List of {allele, seq}"""
    lookup = {}
    
    required_cols = ['Gene', 'Selected Allele', 'Sequence']
    for col in required_cols:
        if col not in df_ref.columns:
            print(f"Error: Reference CSV missing required column '{col}'.")
            sys.exit(1)

    for _, row in df_ref.iterrows():
        g = row['Gene']
        a = row['Selected Allele']
        s = row['Sequence']
        
        if g not in lookup:
            lookup[g] = []
        lookup[g].append({'allele': a, 'seq': s})
    
    return lookup

def main():
    args = parse_args()
    
    print(f"Reading clones file: {args.clones} ...")
    try:
        df_clones = pd.read_csv(args.clones, sep='\t')
    except Exception as e:
        print(f"Failed to read TSV: {e}")
        sys.exit(1)

    print(f"Reading reference library: {args.ref} ...")
    try:
        df_ref = pd.read_csv(args.ref)
    except Exception as e:
        print(f"Failed to read reference CSV: {e}")
        sys.exit(1)

    print(f"Extracting gene list (Mode: {args.gene})...")
    unique_genes = extract_genes_from_tsv(df_clones, target_segment=args.gene)
    print(f"Found {len(unique_genes)} unique genes in clones file (based on *00 matching).")

    ref_lookup = build_reference_lookup(df_ref)

    results = []
    sorted_genes = sorted(list(unique_genes), key=lambda x: (int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 0, x))

    for gene in sorted_genes:
        best_match = None
        if gene in ref_lookup:
            candidates = ref_lookup[gene]
            match = next((c for c in candidates if c['allele'].endswith('*01')), None)
            if not match and candidates:
                match = candidates[0]
            best_match = match
        
        if best_match:
            results.append({
                'gene': gene,
                'allele_A': best_match['allele'],
                'seq_A': best_match['seq'],
                'allele_B': best_match['allele'], # Homozygous format: B = A
                'seq_B': best_match['seq']
            })
        else:
            print(f"Warning: No sequence found for gene {gene} in reference library.")
            results.append({
                'gene': gene,
                'allele_A': None,
                'seq_A': None,
                'allele_B': None,
                'seq_B': None
            })

    df_out = pd.DataFrame(results, columns=['gene', 'allele_A', 'seq_A', 'allele_B', 'seq_B'])
    df_out.to_csv(args.output, index=False)
    print(f"Processing complete! Results saved to: {args.output}")

if __name__ == "__main__":
    main()