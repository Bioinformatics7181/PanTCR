# 08 Pseudo-Bulk RNA-seq Benchmark

Purpose: support Figure 4E and Supplementary Tables S17-S19 using paired 10x Gene Expression and targeted V(D)J datasets.

Implementation scripts:

- `code/experiments/08_pseudo_bulk_rnaseq_benchmark/scbulk_fine_evidence_analysis.py`: precursor truth-status and retained-evidence diagnostic analysis for pseudo-bulk PanTCR results.
- `code/experiments/08_pseudo_bulk_rnaseq_benchmark/audit_fig4e_nonreference_variants.py`: rebuilds the final S19 retained-evidence versus graph-imputed recovered-allele examples.
- `code/experiments/08_pseudo_bulk_rnaseq_benchmark/summarize_figure4_observed_splits.py`: rebuilds the final S18/S19 source CSVs. S18 uses unfiltered MiXCR-derived observed-region coverage (`min_naive=None`) for method comparison; S19 uses retained PanTCR evidence from `audit_fig4e_nonreference_variants.py`.
- `code/experiments/08_pseudo_bulk_rnaseq_benchmark/summarize_figure4_evidence_s18_s19.mjs`: writes compact CSV/Markdown S18/S19 result-summary artifacts from `generated/final_two_table_audit/` after the upstream evidence analyses have been rerun.

Evaluation setup: targeted V(D)J calls are used as matched proxy truth for
evaluation. The Gene Expression-derived TCR observations serve as fragmented
input for PanTCR inference, while the graph prior comes from the AIRR-1-derived
training resource.

Key parameters: the Figure 4E/S19 retained-evidence audit uses the pseudo-bulk PanTCR V-gene results generated with the package default pseudo-bulk settings and retained PanTCR observations filtered with `NaiveDiversityIndex >= 2`. Supplementary Table S18 uses unfiltered MiXCR-derived observed-region coverage for the method-comparison observed-site stratum, so singleton/low-diversity observations can define physical callable coverage there without being counted as retained PanTCR evidence.

Outputs: rerun products are written under `generated/evidence_analysis`, `generated/fig4e_nonreference_variant_audit`, `generated/final_two_table_audit`, and `generated/report`. The files `generated/final_two_table_audit/method_comparison_by_observed_defining_site.csv` and `generated/final_two_table_audit/pantcr_recovered_scbulk_examples.csv` are the regenerated source tables for Supplementary Tables S18 and S19.

Inputs: pseudo-bulk per-dataset source materials default to `experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/03_pseudo_bulk_rnaseq/generated/per_dataset_results`, with references from `data/ref`. The semi-synthetic truth-status input is expected at `experiments/07_semi_synthetic_airr_benchmark/generated/evidence_analysis/per_truth_call_status.csv`; if it is missing, `run_rebuild.sh` generates it from `experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr/generated/results`.

## Run

First generate the upstream semi-synthetic and pseudo-bulk result roots:

```bash
cd experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/02_semi_synthetic_airr
bash run_all_airr_benchmarks.sh

cd ../03_pseudo_bulk_rnaseq
bash run_all_datasets.sh
cd ../../../..
```

Then rebuild the Figure 4E and S18/S19 analysis sources:

```bash
bash experiments/08_pseudo_bulk_rnaseq_benchmark/run_rebuild.sh
```

The pseudo-bulk/scBulk upstream pipeline is documented in
`experiments/01_in_silico_trbv_benchmarks/benchmark_pipelines/03_pseudo_bulk_rnaseq/README.md`.
