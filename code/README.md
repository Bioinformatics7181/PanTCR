# Code

This directory contains the reusable Python implementations used by PanTCR and
by the manuscript reproduction workflows. Numbered experiment folders under
`../experiments/` call these files through their bash drivers.

## Directory Guide

- `pantcr/`: user-facing PanTCR components for MiXCR clone-table processing,
  population-aware graph construction, Bayesian allele inference, and
  post-inference evidence auditing.
- `simulation/`: genotype sampling and controlled synthetic repertoire
  simulation used by the in silico and leave-target workflows.
- `code/experiments/00_benchmark_utils/`: shared benchmark utilities, including
  MiXCR/FindAlleles converters, label-based allele-call evaluation, and common
  summary helpers.
- numbered subdirectories under `code/experiments/`: experiment-specific
  builders for tables, evidence audits, and source summaries. These scripts are
  tied to manuscript reproduction and are not part of the core PanTCR
  application API.
