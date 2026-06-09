# 03 Leave Target Allele Out

Purpose: support Supplementary Tables S12 and S13 by testing recovery of selected real TRBV target alleles when the exact target allele path is absent from the graph-prior construction resource.

Driver and implementation:

- `run_leaveout_allele_specific.sh`: driver for allele-specific leave-target simulations. It calls the central implementations under `code/`.
- `code/experiments/03_leave_target_allele_out/prepare_leaveout_experiment.py`: prepares target-specific leave-target-allele-out inputs.
- `code/experiments/03_leave_target_allele_out/summarize_leaveout_target_recovery.py`: per-target recovery summarizer using the final observed-region compatible partial-recovery definition.
- `code/experiments/03_leave_target_allele_out/summarize_leaveout_allele_recovery.py`: builds the final S12/S13 source tables from the generated leaveout summary directory.
- `supplementary_table_s12_overall.csv`, `supplementary_table_s13_by_target.csv`: final lightweight source tables matching Supplementary Tables S12 and S13.
- `pmtr_default_rule_audit.csv`: audit of default/reference rule handling for the selected targets.

Scope: related local states or same-gene paths may remain available; the stress test removes the exact target allele path, not all related sequence context.

Inputs: full reruns require the simulation reference files under `data/ref/` plus an uncompressed `cdr3_dict.pkl`. From the package root, unpack the bundled archive with `unzip data/ref/cdr3_dict.zip -d data/ref/`, or set `PANTCR_SIM_REF_ROOT` to a full-provenance directory containing `cdr3_dict.pkl`.

## Run

Run the leave-target-allele-out workflow from the package root:

```bash
bash experiments/03_leave_target_allele_out/run_leaveout_allele_specific.sh
```

The workflow prepares target-specific simulation inputs, removes the exact
target allele path from the graph-prior construction resource, reruns inference,
and summarizes exact and observed-region compatible recovery.

Final Supplementary Tables S12/S13 use the observed-region compatible partial definition. Under this definition PanTCR has 75 exact recoveries and one observed-region compatible partial recovery, and PanTCR-NP has 74 exact recoveries and one observed-region compatible partial recovery. Stricter diagnostics that require at least two covered defining positions are reported as QC diagnostics.

Outputs: full reruns write S12/S13 source CSVs under `runs/results/leaveout/<expr>/summary/`.
