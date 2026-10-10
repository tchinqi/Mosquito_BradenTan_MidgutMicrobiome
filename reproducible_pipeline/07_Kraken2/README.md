# 07 - Kraken2 cross-check of unresolved OTUs (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `41_extract_final_GG80_rescue_FASTAs.py` | Extracts OTUs unassigned or assigned only to d__Bacteria by Greengenes2 at 80% (LC 4,884; MC 188). |
| `42_kraken2_FINAL_GG80_rescue_default_vs_conf097.pbs` | Classifies the extracted OTUs with Kraken2 (Standard-8 database) at confidence 0 (default) and 0.97. |
| `43_summarise_FINAL_GG80_Kraken_rescue.py` | Summarises Kraken2 results by deepest rank and taxon; Fig. 1B, Tables S2C-D. |

**Naming note.** Historical file/folder names swap LC and MC: '0p1pct' / 'MC_0p1pct' = final LC (0.1%); 'LC_max0p5pct_20reads' = final MC (max(0.5%, 20 reads)).
