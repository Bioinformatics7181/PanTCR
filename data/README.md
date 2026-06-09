# Data

This directory contains lightweight data used by the PanTCR example and by
package-local reproduction workflows.

## Directory Guide

- `ref/`: shared TRB references, mature-region coordinate files, metadata, the
  cleaned pmTR allele catalog, the MiXCR-all custom library, and the compressed
  CDR3 dictionary used by simulation workflows.
- `examples/sc1/`: a small MiXCR-processed TRB clone table for demonstrating
  one-sample PanTCR inference.
- `pretrained_graphs/trb_airr1/`: bundled pretrained TRB graph material used by
  the example command in the root README.

For paper reproduction, place large raw reads and full MiXCR intermediate trees
under the `input/` directory documented by the relevant experiment folder.
Reruns write outputs under that experiment's `generated/` or `results/`
directory.

Obtain large public inputs from the sources reported in the manuscript Data
Availability statement and Supplementary Table S1. The paper lists AIRR-1 as ArrayExpress `E-MTAB-13593`, AIRR-2 as ENA
`PRJEB28370`, SC-1 through SC-8 as official 10x Genomics datasets, BULK-1 as
GEO `GSE112087`, and BULK-2 as GEO `GSE107011`.
