#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize public bulk PanTCR inference calls into candidate allele tables.

It trims inferred V/J sequences to the mature evaluation interval with
``TRB_index.csv``, matches trimmed sequences against the clean pmTR reference,
and reports matched/unmatched candidate rows consumed by experiment 05.
"""

import argparse
import glob
import os
import sys

import pandas as pd

def load_trb_index(filepath):
    """
    加载 TRB_index.csv
    逻辑更新：
    V 基因：计算长度差值 cut_pos = CDR3Begin - FR1Begin (即 0:end-start)
    J 基因：保留 FR4Begin 用于后缀截断
    """
    print(f"[INFO] Loading Index file: {filepath}")
    if not os.path.exists(filepath):
        sys.exit(f"[ERROR] Index file not found: {filepath}")

    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read index file: {e}")

    index_map = {}
    for _, row in df.iterrows():
        fam = str(row['baseGene']).strip()
        gene_type = str(row['geneType']).strip().upper()
        
        cut_pos = None
        if gene_type == 'V':
            # 逻辑：因为序列从 FR1 开始，所以截取范围是 [0, CDR3Begin - FR1Begin)
            start = row.get('FR1Begin', 0)
            end = row.get('CDR3Begin')
            if pd.notna(end):
                cut_pos = int(float(end)) - int(float(start))
        elif gene_type == 'J':
            # 逻辑：保留 [FR4Begin:]
            val = row.get('FR4Begin')
            if pd.notna(val):
                cut_pos = int(float(val))
        
        if cut_pos is not None:
            index_map[fam] = {'type': gene_type, 'cut_pos': cut_pos}
    return index_map

def trim_sequence(seq, gene_type, cut_pos):
    """
    根据基因类型和切割位点截断序列。
    V: seq[:cut_pos] -> [0 : end-start]
    J: seq[cut_pos:] -> [FR4_begin : ]
    """
    if not seq or pd.isna(seq): return None
    seq = str(seq).strip().upper()
    
    if gene_type == 'V':
        # 保留从 0 到差值长度的部分
        return seq[:cut_pos] if cut_pos <= len(seq) else seq
    elif gene_type == 'J':
        # 保留从 FR4Begin 开始到结尾的部分
        return seq[cut_pos:] if cut_pos < len(seq) else ""
    return None

def build_processed_reference(ref_path, index_map):
    """
    预处理参考库：生成 {trimmed_seq: [allele_info]}
    """
    print(f"[INFO] Pre-processing Reference alleles from: {ref_path}")
    if not os.path.exists(ref_path):
        # 如果没有参考文件，我们可以返回空字典，意味着所有序列都会变成 Unmatched
        print(f"[WARN] Reference file not found: {ref_path}. All sequences will be treated as Unmatched.")
        return {'V': {}, 'J': {}}
    
    try:
        ref_df = pd.read_csv(ref_path)
    except Exception as e:
        sys.exit(f"[ERROR] Failed to read reference file: {e}")

    processed_ref = {'V': {}, 'J': {}}
    
    count_loaded = 0
    for _, row in ref_df.iterrows():
        allele = str(row['allele']).strip()
        full_seq = str(row['sequence']).strip().upper()
        family = allele.split('*')[0] if '*' in allele else allele
        
        if family in index_map:
            info = index_map[family]
            trimmed = trim_sequence(full_seq, info['type'], info['cut_pos'])
            
            if trimmed:
                target_dict = processed_ref[info['type']]
                if trimmed not in target_dict:
                    target_dict[trimmed] = []
                
                target_dict[trimmed].append({
                    'allele': allele,
                    'family': family,
                    'is_new_noCDR3': row.get('is_new_noCDR3', ''),
                    'is_default_noCDR3': row.get('is_default_noCDR3', '')
                })
                count_loaded += 1
    
    print(f"[INFO] Loaded and trimmed {count_loaded} reference alleles.")
    return processed_ref

def process_samples_and_aggregate(file_list, gene_type, index_map, ref_lookup):
    """
    处理样本文件。
    逻辑更新：
    1. 截断序列。
    2. 尝试匹配参考库。
    3. 如果匹配不到，则记录预测的家族 (pred_family) 并标记为 Unmatched。
    4. 无论是否匹配，都按 'trimmed_sequence' 进行聚合。
    """
    # 结果字典结构: 
    # trimmed_sequence -> {
    #    'sample_count': int, 
    #    'files': set(), 
    #    'matches': list (参考库信息),
    #    'pred_families': set (样本中预测的家族，用于 Unmatched 情况)
    # }
    results_map = {}
    
    for fpath in file_list:
        fname = os.path.basename(fpath)
        try:
            df = pd.read_csv(fpath)
            seen_in_this_file = set()

            if 'gene' not in df.columns:
                continue

            for _, row in df.iterrows():
                pred_family = str(row['gene']).strip()
                if pred_family not in index_map:
                    continue
                
                info = index_map[pred_family]
                if info['type'] != gene_type:
                    continue

                cut_pos = info['cut_pos']
                
                # 获取序列
                raw_seqs = []
                if 'seq_A' in row and pd.notna(row['seq_A']): raw_seqs.append(row['seq_A'])
                if 'seq_B' in row and pd.notna(row['seq_B']): raw_seqs.append(row['seq_B'])
                
                for rs in raw_seqs:
                    trimmed = trim_sequence(rs, gene_type, cut_pos)
                    
                    if not trimmed:
                        continue
                    
                    # 初始化聚合记录
                    if trimmed not in results_map:
                        results_map[trimmed] = {
                            'sample_count': 0, 
                            'files': set(), 
                            'matches': [],
                            'pred_families': set()
                        }
                    
                    # 尝试匹配参考库
                    ref_matches = ref_lookup[gene_type].get(trimmed)
                    
                    if ref_matches:
                        # 如果匹配，存入参考库信息（可能会覆盖之前的空列表，或者追加）
                        # 注意：对于同一个 trimmed 序列，参考库的匹配结果是固定的。
                        # 我们只需要赋值一次即可，但为了简单，这里直接赋值
                        results_map[trimmed]['matches'] = ref_matches
                    else:
                        # 如果没匹配，记录样本中预测的家族
                        results_map[trimmed]['pred_families'].add(pred_family)
                        
                    # 计数逻辑
                    if trimmed not in seen_in_this_file:
                        results_map[trimmed]['sample_count'] += 1
                        results_map[trimmed]['files'].add(fname)
                        seen_in_this_file.add(trimmed)
                            
        except Exception as e:
            print(f"[WARN] Error processing file {fname}: {e}")
            
    return results_map

def main():
    parser = argparse.ArgumentParser(description="TCR Allele Analysis (Matched & Unmatched)")
    parser.add_argument('--input_dir', required=True, help="Folder with V/ and J/ subfolders")
    parser.add_argument('--index', default='TRB_index.csv', help="Path to TRB_index.csv")
    parser.add_argument('--ref', default='pmTR_TRB_V_J_cleaned.csv', help="Path to pmTR reference CSV")
    parser.add_argument('--output', default='candidate_allele_summary.csv', help="Output CSV filename")
    
    args = parser.parse_args()
    
    index_map = load_trb_index(args.index)
    ref_lookup = build_processed_reference(args.ref, index_map)
    
    all_rows = []
    
    for gtype in ['V', 'J']:
        subdir = os.path.join(args.input_dir, gtype)
        if not os.path.exists(subdir):
            continue
            
        csv_files = glob.glob(os.path.join(subdir, "*.csv"))
        print(f"[INFO] Processing {len(csv_files)} {gtype} files...")
        
        type_results = process_samples_and_aggregate(csv_files, gtype, index_map, ref_lookup)
        
        for trimmed_seq, data in type_results.items():
            matches = data['matches']
            pred_families = data['pred_families']
            
            row = {
                'sequence': trimmed_seq,
                'sample_counts': data['sample_count'],
                'filenames': "|".join(sorted(list(data['files'])))
            }
            
            if matches:
                # === 情况 A: 匹配到了参考库 ===
                # 使用参考库的信息
                unique_families = sorted(list(set(m['family'] for m in matches)))
                unique_alleles = sorted(list(set(m['allele'] for m in matches)))
                unique_is_new = sorted(list(set(str(m['is_new_noCDR3']) for m in matches)))
                unique_is_default = sorted(list(set(str(m['is_default_noCDR3']) for m in matches)))
                
                row['family'] = "|".join(unique_families)
                row['allele'] = "|".join(unique_alleles)
                row['is_new_noCDR3'] = "|".join(unique_is_new)
                row['is_default_noCDR3'] = "|".join(unique_is_default)
                row['match_status'] = 'Matched'
                
            else:
                # === 情况 B: 没有匹配到 ===
                # 使用样本预测的家族
                row['family'] = "|".join(sorted(list(pred_families)))
                row['allele'] = 'Unmatched' 
                row['is_new_noCDR3'] = ''
                row['is_default_noCDR3'] = ''
                row['match_status'] = 'Unmatched'
            
            all_rows.append(row)

    if all_rows:
        res_df = pd.DataFrame(all_rows)
        # 调整列顺序，增加 match_status 方便筛选
        cols = ['family', 'allele', 'sequence', 'sample_counts', 'match_status', 'is_new_noCDR3', 'is_default_noCDR3', 'filenames']
        # 保证列存在
        for c in cols:
            if c not in res_df.columns: res_df[c] = ''
            
        res_df = res_df[cols]
        res_df.to_csv(args.output, index=False)
        print(f"[SUCCESS] Summary saved to: {args.output}")
        print(f"[INFO] Total unique sequences: {len(res_df)}")
        print(f"[INFO] Matched: {len(res_df[res_df['match_status']=='Matched'])}")
        print(f"[INFO] Unmatched: {len(res_df[res_df['match_status']=='Unmatched'])}")
    else:
        print("[WARN] No data found.")

if __name__ == "__main__":
    main()
