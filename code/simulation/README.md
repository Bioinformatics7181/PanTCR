# Simulation Code

This directory contains the reusable simulation components used by the
truth-known in silico benchmarks and leave-target workflows.

## Files

- `select_genotype.py`: samples diploid V/J genotypes from the cleaned pmTR TRB
  allele catalog. It also supports the Scenario C introduced-variant setting
  through mutation parameters passed by the experiment drivers.
- `simulate_repertoire.py`: generates rearranged transcript-like V-CDR3-J
  molecules from simulated genotypes and an empirical CDR3 dictionary.
- `simulate_bias_advanced.py`: applies molecule-level 5' truncation and
  degradation settings before read simulation.

These scripts are called by the bash drivers in
`../../experiments/01_in_silico_trbv_benchmarks/`.
The top-level experiment README documents the manuscript-facing scenario
parameters and output locations.
