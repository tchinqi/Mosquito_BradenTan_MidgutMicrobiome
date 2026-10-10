# 04 - Abundance-support filtering: LC and MC datasets (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `16_filter_final_racon_postracon_thresholds.py` | Applies per-sample support thresholds (0.1%, 0.2%, 0.5%, 5, 10, 20 reads) to the 22,528 pooled Racon consensuses. The 0.1% set is the final LC dataset. |
| `17_build_postracon_threshold_per_sample_table.py` | Per-sample table of consensuses and supporting reads retained under each threshold. |
| `18_build_LC_max0p5pct_20reads.py` | Builds the final MC dataset: consensuses supported by >= max(0.5% of the sample's reads, 20 reads). |

**Naming note.** Historical file/folder names swap LC and MC: '0p1pct' / 'MC_0p1pct' = final LC (0.1%); 'LC_max0p5pct_20reads' = final MC (max(0.5%, 20 reads)).
