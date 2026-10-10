# 08 - Final OTU x sample table (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `44_make_MC_GG97_classified_OTU_sample_table.py` | Builds the supporting-read table of the 328 MC OTUs classified by Greengenes2 at >= 97% across the 576 wells (reads per consensus summed per OTU). |

**Naming note.** Historical file/folder names swap LC and MC: '0p1pct' / 'MC_0p1pct' = final LC (0.1%); 'LC_max0p5pct_20reads' = final MC (max(0.5%, 20 reads)).
