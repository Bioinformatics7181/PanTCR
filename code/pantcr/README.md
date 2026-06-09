# PanTCR Application Components

This directory contains the user-facing PanTCR application layer. It is separated from the paper-specific `experiments/` directory so that users can run allele inference without mixing the core workflow with benchmark evaluation code.

## Main Components

- `collect_mutation.py`: converts a MiXCR-processed rearranged TCR clone table into gene-assigned, coordinate-projected V/J observation tables.
- `build_pangenome_graph.py`: builds population-aware local pan-allele graphs from multiple observation tables and sample metadata.
- `infer_genotype_bayes.py`: performs graph-guided Bayesian V/J allele inference for one sample.
- `audit_inference_evidence.py`: annotates inferred allele sequences after inference, separating retained-evidence-covered regions from graph-imputed regions and classifying default/reference-relative changes.
- `run_sample.py`: convenience wrapper for one-sample inference using a pre-trained graph.

Allele-recovery scoring, result builders, and benchmark audits live under `experiments/`.

## Default User Workflow

1. Run MiXCR with an RNA-seq or suitable repertoire workflow to obtain a `clones_TRB.tsv` file.
2. Run PanTCR on the MiXCR clone table:

```bash
python code/pantcr/run_sample.py \
  --clones data/examples/sc1/SC1.clones_TRB.tsv \
  --sample-id SC1 \
  --out-dir data/examples/sc1/output \
  --graph-root data/pretrained_graphs/trb_airr1 \
  --reference-dir data/ref
```

The command writes V and J observation tables, inference CSVs, and evidence/imputation audit CSVs.

## Building a Custom Graph

For a cohort-specific graph, first run `collect_mutation.py` for each graph-construction sample, keeping V and J outputs in separate directories. Then run:

```bash
python code/pantcr/build_pangenome_graph.py \
  --in_dir population_observations/V \
  --out_dir custom_graph/V \
  --metadata metadata.csv \
  --write_paths

python code/pantcr/build_pangenome_graph.py \
  --in_dir population_observations/J \
  --out_dir custom_graph/J \
  --metadata metadata.csv \
  --write_paths
```

The resulting `custom_graph/V` and `custom_graph/J` directories can be passed to `run_sample.py` through `--graph-root custom_graph` if both subdirectories are present.

The metadata CSV contains at least `sample_id` and `population_id`. Add a `filename` column when the observation-table filename differs from `sample_id`:

```csv
sample_id,population_id,filename
sample_001,AFR,sample_001.V.csv
sample_002,EUR,sample_002.V.csv
```
