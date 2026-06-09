# Main Benchmark Pipelines

This directory is the upstream benchmark/application workflow hub for results
reused by later audit folders. It centralizes the benchmark-generating steps
used by the downstream experiment drivers.

Canonical generated outputs:

- `02_semi_synthetic_airr/generated/results`: semi-synthetic AIRR benchmark result root used by experiment 07.
- `03_pseudo_bulk_rnaseq/generated/per_dataset_results`: pseudo-bulk/scBulk benchmark result root used by experiment 08.
- `04_public_bulk_rnaseq/generated/BULK-1` and `generated/BULK-2`: public bulk application result roots used by experiment 05.

Workflow coverage:

- In silico TRBV simulation is implemented by the parent `01_in_silico_trbv_benchmarks` scripts.
- Semi-synthetic AIRR processing is implemented under `02_semi_synthetic_airr/`. It derives proxy genotype labels from the full-length MiXCR FindAlleles outputs (`*.alleles.tsv` and `*.customAlleles.json`) supplied in that pipeline's input tree, then runs the degraded-read benchmark.
- Pseudo-bulk/scBulk processing is implemented under
  `03_pseudo_bulk_rnaseq/`; place full raw 10x inputs under the subfolder
  `input/` paths described there.
- Public bulk SRA/MiXCR/mutation/graph/inference processing is implemented
  under `04_public_bulk_rnaseq/`; reruns require SRA Toolkit, MiXCR, and the
  public accession list.
