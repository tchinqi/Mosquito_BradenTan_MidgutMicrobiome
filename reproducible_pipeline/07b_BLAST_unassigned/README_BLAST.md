# 07b – BLAST of OTUs left unresolved by Greengenes2 (10 Oct 2026)

Supplementary check (Methods 2.9; Table S2E). Not part of the abundance table:
the final dataset still uses only the 328 MC OTUs classified by Greengenes2 at ≥97%.

## Input
`MC_GG80_Unassigned_plus_domain_only.fasta` – the 188 MC OTUs that Greengenes2 at
≥80% identity left unassigned or assigned only to d__Bacteria (same set as the
Kraken2 cross-check, script 41; Atlas9
`06_taxonomy/kraken2_final/FINAL_GG80_unassigned_domain_rescue/`).

## Searches (NCBI web BLAST, 10 Oct 2026)
1. Nucleotide BLAST, database **core_nt**, program **megablast**, expect threshold
   1e-10, maximum 10 target sequences; all 188 OTUs. Downloaded as
   *Download All → Single-file JSON* (`web_blast.json`).
2. The 83 OTUs with no megablast hit (`nohit_83.fasta`) were re-searched with
   **blastn** (somewhat similar sequences), same database and settings
   (`web_blast_nohit.json`).

The JSON downloads are not stored here (large; reproducible from the inputs above).

## Summary
`summarise_web_blast.py` (Python 3, run on Atlas9 in `08_unassigned_blast/`):
- best hit per OTU = highest bit score; % identity of its best HSP;
- query coverage = union of all HSPs to that subject / OTU length;
- category from the full NCBI lineage of the hit's taxid (BLAST+
  `taxonomy4blast.sqlite3`): Human, Culicidae, other arthropod, other vertebrate,
  Fungi, Viridiplantae, other eukaryote, Archaea, Virus, Bacteria (split at 97%
  identity), no hit;
- `possible_chimera` = best hit covers < 70% of the OTU;
- reads and samples from
  `07_analysis/MC_GG2_97_classified_OTU_table/MC_867_OTUs_x_576_samples_supporting_reads.tsv`
  (all 576 wells).

```
python3 summarise_web_blast.py web_blast.json        # -> unassigned_188_top_hits_megablast.tsv
python3 summarise_web_blast.py web_blast_nohit.json  # -> nohit_83_top_hits_blastn.tsv
```
(The script always writes `unassigned_188_top_hits_web.tsv`; the outputs were
renamed as above. In the blastn file, OTUs not in that search appear as "No hit".)

## Result (88,915 supporting reads in the 188 OTUs; 5.3% of all MC reads)
| Category | OTUs | % reads |
|---|---|---|
| Mosquito (Culicidae), megablast | 90 | 21.5 |
| Mosquito (Culicidae), blastn | 15 | 22.5 |
| *Elizabethkingia* / *Wolbachia* / *Sphingomonas* (98.1–99.8% identity) | 10 | 27.9 |
| Other bacteria | 3 | 0.6 |
| Fish (weak, 68–91% identity) | 3 | 1.0 |
| Virus (partial, Mampikony meta-like virus) | 4 | 1.1 |
| No significant hit | 63 | 25.4 |

96% of mosquito-assigned reads had *Aedes albopictus* as the best hit.
