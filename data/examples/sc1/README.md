# SC1 Example

This example contains a processed MiXCR TRB clone table and a matched genotype label file used for demonstration and optional evaluation.

Input files:

- `SC1.clones_TRB.tsv`: MiXCR-processed rearranged TRB clone table.
- `SC1.genotype.csv`: matched genotype/proxy-truth label used only for evaluation examples.

Run PanTCR with the bundled AIRR-1-derived pre-trained graph:

```bash
python code/pantcr/run_sample.py \
  --clones data/examples/sc1/SC1.clones_TRB.tsv \
  --sample-id SC1 \
  --out-dir data/examples/sc1/output \
  --graph-root data/pretrained_graphs/trb_airr1 \
  --reference-dir data/ref
```

Main outputs:

- `data/examples/sc1/output/SC1.V_sequences.csv` and `data/examples/sc1/output/SC1.J_sequences.csv`: coordinate-projected observations retained before inference filtering.
- `data/examples/sc1/output/infer_SC1.V.csv` and `data/examples/sc1/output/infer_SC1.J.csv`: PanTCR inferred V/J allele sequences.
- `data/examples/sc1/output/audit_SC1.V.summary.csv` and `data/examples/sc1/output/audit_SC1.J.summary.csv`: retained-evidence-covered and graph-imputed intervals for each inferred allele sequence.
- `data/examples/sc1/output/audit_SC1.V.site_detail.csv` and `data/examples/sc1/output/audit_SC1.J.site_detail.csv`: default/reference-relative change-level evidence audit.

The label file is optional for PanTCR inference and shows the benchmark-style evaluation input used in the manuscript.
