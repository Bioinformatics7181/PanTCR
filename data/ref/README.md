# Shared Reference Files

This directory contains the reference files shared by the PanTCR example,
simulation workflows and baseline converters.

## Files

- `pmTR_TRB_V_J_cleaned.csv`: canonical cleaned pmTR TRB V/J allele catalog for
  reproduction workflows.
- `IMGT_TRB_default.csv`: MiXCR-default allele reference mapping used by
  coordinate and default/non-default reporting steps.
- `IMGT_TRBV_pro.tsv` and `IMGT_TRBJ_pro.tsv`: TRBV/TRBJ reference sequences
  used during observation extraction.
- `TRB_index.csv`: mature-region coordinate intervals used for V/J evaluation
  trimming.
- `metadata.csv`: sample metadata used by bundled graph construction examples.
- `hsa_custom_library.json`: custom MiXCR library used for the MiXCR-all
  expanded-reference baseline.
- `cdr3_dict.zip`: compressed empirical CDR3 dictionary for in silico
  repertoire simulation.

Before running simulation workflows, unpack the CDR3 dictionary from the
package root:

```bash
unzip data/ref/cdr3_dict.zip -d data/ref/
```

This creates:

```text
data/ref/cdr3_dict.pkl
```

Use `pmTR_TRB_V_J_cleaned.csv` as the package allele catalog for reproduction
workflows.
