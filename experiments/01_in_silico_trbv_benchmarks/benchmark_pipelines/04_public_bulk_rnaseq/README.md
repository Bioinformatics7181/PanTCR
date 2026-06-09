# Public Bulk RNA-seq Application Pipeline

Purpose: generate public bulk RNA-seq application outputs consumed by `experiments/05_public_bulk_candidate_nonimgt`.

Canonical output paths:

- `experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/04_public_bulk_rnaseq/generated/BULK-1`
- `experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/04_public_bulk_rnaseq/generated/BULK-2`

Package-local drivers:

- `run_sra_to_clones.sh`: driver from public SRA accessions to MiXCR
  `*.clones_TRB.tsv` files.
- `run_from_clones.sh`: driver from existing MiXCR `*.clones_TRB.tsv` files to
  the mutation and PanTCR inference tree expected by experiment 05.

Reusable implementation code is under `code/pantcr/` and
`code/experiments/05_public_bulk_candidate_nonimgt/`.

Workflow:

1. Read SRA accessions from `SRR_Acc_List.txt`.
2. Download SRA reads with `prefetch` and convert to FASTQ with `fasterq-dump`.
3. Run `mixcr analyze rna-seq` and collect `*.clones_TRB.tsv` into a sample directory.
4. Build V/J mutation-observation CSVs with `code/pantcr/collect_mutation.py`.
5. Build or reuse the configured graph roots and place the selected graph under `generated/{BULK}/pang/semi/{V,J}` for downstream evidence auditing.
6. Run PanTCR Bayesian inference for all samples.
7. Summarize candidate non-IMGT alleles with
   `code/experiments/05_public_bulk_candidate_nonimgt/summarize_public_bulk_candidates.py`
   and write `candidate_allele_summary.csv`.

Required tools:

- SRA Toolkit (`prefetch`, `fasterq-dump`)
- `mixcr`
- Python dependencies from `requirements.txt`

Dataset source:

The public bulk cohorts are listed in the manuscript Data Availability
statement and Supplementary Table S1: BULK-1 is GEO `GSE112087`, and BULK-2 is
GEO `GSE107011`. Prepare `input/BULK-1/SRR_Acc_List.txt` and
`input/BULK-2/SRR_Acc_List.txt` from the GEO/SRA run tables for those studies,
then use the commands below to download reads and rebuild the PanTCR result
trees.

Downstream use:

`experiments/05_public_bulk_candidate_nonimgt/run_rebuild.sh` reads `generated/BULK-1` and `generated/BULK-2` from this pipeline. Public-bulk raw reads and large intermediate files follow the input and generated result layout shown above.

Example rerun from public SRA accessions through PanTCR inference:

```bash
bash run_sra_to_clones.sh \
  --bulk-label BULK-1 \
  --acc-list input/BULK-1/SRR_Acc_List.txt

bash run_from_clones.sh \
  --bulk-label BULK-1 \
  --clones-dir generated/BULK-1/clones
```

`run_from_clones.sh` writes `candidate_allele_summary.csv` under each bulk
output root by default. Use `--summary-output` to override the summary filename.
