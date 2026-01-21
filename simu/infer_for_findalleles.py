import argparse
import csv
import json
import re
import sys

# ---------------- Core Functions ----------------

def parse_mutations(mut_string: str):
    """
        Parses mutations (Compact Format) and adjusts coordinates.
        Formats handled:
          - Substitution: ST6A (Ref T at 6 -> A)
          - Deletion:     DA28 (Ref A at 28 deleted)
          - Insertion:    I28C (Insert C at 28)
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

    return muts

def ref2target(ref, mutations):
    """
    Target sequence reconstruction via reference and mutation list.
    """
    ref = list(ref)
    
    mutations_desc = mutations[::-1]
    
    for m in mutations_desc:
        p = m['pos']

        if p >= len(ref) or p < 0:
            continue

        if m['type'] == 'S':
            ref[p] = m['target']
        elif m['type'] == 'D':
            del ref[p : p + len(m['seq'])]
        elif m['type'] == 'I':
            ref[p: p] = list(m['seq'])

    return "".join(ref)


# ---------------- Main----------------

def main():
    parser = argparse.ArgumentParser(description="Infer Allele Sequences from MiXCR Results")
    parser.add_argument('--tsv', required=True, help='Path to the .alleles.tsv file')
    parser.add_argument('--json', required=True, help='Path to the .customAlleles.json file')
    parser.add_argument('--ref', required=True, help='Path to the IMGT reference .csv file')
    parser.add_argument('--gene', choices=['V', 'J', 'all'], default='all', help='Specify which segment type to process: V, J, D or all (default)')
    parser.add_argument('--output', required=True, help='Path to the output .csv file')
    
    args = parser.parse_args()

    # 1. Load Reference (IMGT CSV)
    print(f"Loading reference sequences from {args.ref}...")
    ref_seqs = {}
    with open(args.ref, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('Gene')
            seq = row.get('Sequence')
            if name and seq:
                ref_seqs[name] = seq.upper()

    # 2. Load json
    print(f"Loading custom alleles from {args.json}...")
    custom_mutations_map = {}  # alleleName -> (parent, mutation_string)
    with open(args.json, 'r') as f:
        json_data = json.load(f)
        genes_data = []
        if isinstance(json_data, list):
            for entry in json_data:
                if 'genes' in entry and 'speciesNames' in entry:
                    if 'hsa' in entry['speciesNames']:
                        genes_data.extend(entry['genes'])
        else:
            if 'genes' in json_data and 'speciesNames' in json_data:
                if 'hsa' in json_data['speciesNames']:
                    genes_data = json_data['genes']

        for g in genes_data:
            name = g.get('name')
            mut = ""
            parent = ""
            
            if 'meta' in g and 'alleleMutations' in g['meta'] and 'alleleVariantOf' in g['meta']:
                mut = g['meta']['alleleMutations']
                parent = g['meta'].get('alleleVariantOf')
            elif 'alleleInfo' in g and 'mutations' in g['alleleInfo'] and 'parent' in g['alleleInfo']:
                mut = g['alleleInfo'].get('mutations', "")
                parent = g['alleleInfo'].get('parent')
            
            if name:
                custom_mutations_map[name] = (parent, mut)

    # 3. Read tsv files
    print(f"Processing TSV file {args.tsv} with gene {args.gene}...")
    gene_alleles = {}
    
    valid_statuses = {'DE_NOVO', 'FOUND_KNOWN_VARIANT', 'NOT_CHANGED_AFTER_SEARCH'}

    with open(args.tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            row_type = row.get('type', '')
            
            if args.gene == 'V' and row_type != 'Variable':
                continue
            if args.gene == 'J' and row_type != 'Joining':
                continue
            
            status = row['status']
            if status not in valid_statuses:
                continue

            gene_name = row['geneName']
            allele_name = row['alleleName']
            clones_count = float(row.get('clonesCount', 0))
            
            final_seq = ""

            if status == 'NOT_CHANGED_AFTER_SEARCH':
                if gene_name in ref_seqs:
                    final_seq = ref_seqs[gene_name]
                else:
                    print(f"[Warn] Ref sequence not found for {allele_name}")
                    continue

            else:
                parent = row.get('varianceOf')
                mut_str = row.get('mutations')
                
                if allele_name in custom_mutations_map:
                    json_parent, json_mut = custom_mutations_map[allele_name]
                    if json_parent: parent = json_parent
                    if json_mut: mut_str = json_mut

                if not parent:
                    parent = gene_name 

                ref_seq = None
                if gene_name in ref_seqs:
                    ref_seq = ref_seqs[gene_name]
                
                if not ref_seq:
                    print(f"[Error] Cannot find parent sequence for {allele_name} (parent: {parent})")
                    continue

                try:
                    muts = parse_mutations(mut_str)
                    final_seq = ref2target(ref_seq, muts)
                except Exception as e:
                    print(f"[Error] Failed to reconstruct {allele_name}: {e}")
                    continue

            if gene_name not in gene_alleles:
                gene_alleles[gene_name] = []
            
            gene_alleles[gene_name].append({
                'name': allele_name,
                'seq': final_seq,
                'count': clones_count
            })

    # 4. Output
    print(f"Writing output to {args.output}...")
    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['gene', 'allele_A', 'seq_A', 'allele_B', 'seq_B'])

        sorted_genes = sorted(gene_alleles.keys())
        
        for gene in sorted_genes:
            alleles = gene_alleles[gene]
            alleles.sort(key=lambda x: x['count'], reverse=True)
            
            top_alleles = alleles[:2]
            
            if not top_alleles:
                continue

            row_data = [gene]
            
            # Allele A
            row_data.append(top_alleles[0]['name'])
            row_data.append(top_alleles[0]['seq'])
            
            # Allele B
            if len(top_alleles) >= 2:
                row_data.append(top_alleles[1]['name'])
                row_data.append(top_alleles[1]['seq'])
            else:
                row_data.append(top_alleles[0]['name'])
                row_data.append(top_alleles[0]['seq'])
            
            writer.writerow(row_data)

    print("Done.")

if __name__ == '__main__':
    main()