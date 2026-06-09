# Shared Benchmark Utilities

This directory contains reusable utilities shared by multiple manuscript
benchmarks.

## Files

- `evaluate_allele_calls.py`: evaluates inferred allele calls against genotype
  labels. Default `sequence` mode is used for PanTCR, MiXCR-default, and
  FindAlleles-style outputs; `mixcr-all` mode is reserved for the
  expanded-reference MiXCR-all baseline.
- `clones2alleles.py`: converts MiXCR clone-table-derived material into a
  baseline-compatible allele-call format.
- `infer_for_findalleles.py`: helper for formatting FindAlleles baseline
  inference outputs for evaluation.
- `summarize_eval_results.py`: summarizes per-sample evaluation files into
  method-level tables used by downstream experiment builders.
- `common_pantcr_io.py`: shared parsing and path helpers used by the
  experiment-specific analysis scripts.

These utilities are called by the numbered experiment drivers under
`../../../experiments/`.
