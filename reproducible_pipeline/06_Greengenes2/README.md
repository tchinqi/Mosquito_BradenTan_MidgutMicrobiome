# 06 - Greengenes2 reference preparation and classification (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `12_gg2_extract_v3v4.pbs` | Extracts the 338F-806R region (350-550 bp, both orientations) from the Greengenes2 2022.10 full-length backbone with QIIME 2. |
| `20_greengenes2_MC_LC_final_uppercase.pbs` | Classifies LC and MC OTU centroids with classify-consensus-vsearch at >= 97% identity, query coverage 0.80, both strands. |
| `24_greengenes2_MC_LC_default80_finalnames.pbs` | Same classification at >= 80% identity (query coverage 0.80), with the corrected LC/MC folder names. |
| `36_compare_GG2_80_97_by_OTU.py` | Pairs each OTU's 80% and 97% assignments (unchanged / reassigned / assigned->unassigned / unassigned->assigned); Table S2B. |
| `37_make_corrected_GG2_Priya_tables.py` | Summary tables of the 80% vs 97% pairing for LC and MC. |
| `40_audit_GG2_80_taxonomic_depth.py` | Counts the deepest rank (species/genus/family/above/unassigned) reached at 80% identity; Fig. 1A / Table S2A. |
| `RECONSTRUCTED_uppercase_fasta.py` | RECONSTRUCTED: converts OTU centroid FASTA to upper case before classification; validated byte-for-byte. |

**Naming note.** Historical file/folder names swap LC and MC: '0p1pct' / 'MC_0p1pct' = final LC (0.1%); 'LC_max0p5pct_20reads' = final MC (max(0.5%, 20 reads)).
