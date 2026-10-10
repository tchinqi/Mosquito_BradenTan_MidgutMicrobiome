# 09 - Metadata, negative controls and R inputs (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `45_join_MC_GG97_table_to_metadata.py` | Joins the OTU table to sample metadata (origin, plate, well, collection site). |
| `48_audit_BT_candidate_negative_controls.py` | Audits read support of candidate BT negative-control wells before correction. |
| `49_correct_BT_negative_controls.py` | Labels the BT negative controls in the metadata (BT_P02 assignment later corrected by script 53); writes a backup first. |
| `50_prepare_MC_GG97_plotting_tables_6controls.py` | Writes the R input tables: 570 biological / 380 wild / 190 lab / 345 location / 6 control matrices x 328 OTUs, and the OTU-to-genus map. |
| `51_summarise_MC_GG97_wild_vs_lab_genera_final.py` | Mean genus relative abundance, wild vs lab (cross-check used by R script 01). |
| `52_verify_final_R_inputs.py` | Read-only check of every R input file against the agreed totals (PASS/FAIL/WARN). |
| `53_swap_BT_P02_negative_control.py` | Corrects script 49: BT_ALB_192 (well H12) is the BT_P02 control; BT_CLX_120 restored as a lab sample. Then re-run 45, 50, 51, 52. |
