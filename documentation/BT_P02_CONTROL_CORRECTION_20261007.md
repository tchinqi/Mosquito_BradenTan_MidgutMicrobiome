# BT_P02 negative-control correction — 2026-10-07

## Problem
Script 49 designated `BT_CLX_120` as the BT_P02 negative control. The
demultiplexing sheet (`01_demux_metadata/demux_metadata_audit.csv`) shows that
BT specimen numbers follow well order and that BT_CLX_120 is in well **H3**.
The final well (**H12**) of BT_P02 is `BT_ALB_192`, which was the blank.

## Fix
`53_swap_BT_P02_negative_control.py`
- backup: `metadata/sample_metadata_before_BTP02_swap.tsv`
- BT_CLX_120 restored from `sample_metadata_before_6control_fix.tsv` (Lab-reared, H3)
- BT_ALB_192 set to Control (H12)
- validation: 380 Wild / 190 Lab-reared / 6 Control — PASS

Re-run afterwards with conda Python 3.11.16
(`/hpctmp/e0958567/conda_envs/ngspeciesid-0.3.1/bin/python`; the system
`python3` resolves to STAR-CCM+ Python 3.6 and fails on UTF-8 input):
45 → 50 → 51 → 52 (`--run $RUN`).

## Script 52 changes
- `EXP_CONTROLS`: BT_CLX_120 → BT_ALB_192
- `EXP_GENUS_MEANS`: Elizabethkingia (25.62, 15.43, 45.99) → (25.59, 15.43, 45.91);
  Wolbachia (20.49, 25.96, 9.55) → (20.63, 25.96, 9.97)
- pre-swap copy: `development_and_superseded/final_reanalysis_nonpipeline/52_verify_final_R_inputs_pre_BTP02_swap.py`

## Result
PASS 44 / FAIL 0 / WARN 1.
- Biological classified reads: 1,422,797; controls: 10,199; total 1,432,996 (unchanged)
- Lab species: *Ae. albopictus* 165, *Culex* 25
- WARN: OTU_441 (*Staphylococcus epidermidis*) has zero reads in all 570
  biological samples — present only in negative controls (expected handling
  contaminant; previously hidden because BT_ALB_192 was treated as a sample).

## Downstream
R scripts 00, 01, 02, 06–10 updated (expected totals and group sizes) and must be re-run.
