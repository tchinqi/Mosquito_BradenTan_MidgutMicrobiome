# LSM4288 FYP — Master Code Recovery Ledger

Date: 2026-10-06 (updated 2026-10-07: BT_P02 negative-control correction; updated 2026-10-08: R analysis stage added)

## Classification key

- ORIGINAL — historical production source code preserved.
- RECOVERED DEPENDENCY — exact historical external/modified dependency recovered.
- RECONSTRUCTED — original interactive source unavailable; workflow reconstructed from preserved outputs/logs.
- RECONSTRUCTED — BYTE-IDENTICAL — reconstruction reproduces historical outputs byte-for-byte.
- SUPERSEDED — historical code/output preserved but not part of the final analysis.
- CURRENT WORK — analysis/code is still being developed and is not historical recovery.

---

## Pipeline provenance

| Stage | Status | Evidence / notes |
|---|---|---|
| Demultiplexing sheet construction | ORIGINAL | `scripts/make_demux_braden.py`; identical to dated 2026-08-20 FINAL snapshot |
| Minibar demultiplexing | ORIGINAL | `scripts/02_run_minibar.pbs`; identical to dated 2026-08-20 FINAL snapshot |
| Intermediate `03_no_gaps_fastq` preparation | NOT SEPARATELY VERIFIED | Input to production length-filter script is preserved, but exact production command creating this intermediate directory has not been independently audited |
| 350–550 bp read filtering | ORIGINAL | `scripts/07_length_filter.pbs`; SeqKit `seq -g -m 350 -M 550`; identical to dated 2026-08-20 FINAL snapshot |
| NGSpeciesID within-sample consensus workflow | ORIGINAL | Final min2 baseline/rescue scripts preserved in final reanalysis scripts |
| Modified NGSpeciesID Racon-failure continuation logic | RECOVERED DEPENDENCY | Exact patched `ngspeciesid_min2_continue` source recovered and snapshotted |
| Final Racon consensus pooling + provenance table | RECONSTRUCTED — BYTE-IDENTICAL | Reconstructed script reproduces historical pooled FASTA and provenance TSV exactly |
| Post-Racon LC threshold: 0.1% per sample | ORIGINAL | Historical threshold workflow preserved; 18,751 retained consensuses |
| Post-Racon MC threshold: max(0.5%, 20 reads) | ORIGINAL | Historical threshold workflow preserved; 7,375 retained consensuses |
| Final LC VSEARCH97 clustering | RECONSTRUCTED — BYTE-IDENTICAL | 7,642 OTUs; reconstructed centroid FASTA and UC file exactly match historical outputs |
| Final MC VSEARCH97 clustering | RECONSTRUCTED — BYTE-IDENTICAL | 867 OTUs; reconstructed centroid FASTA and UC file exactly match historical outputs |
| Uppercase OTU centroid conversion | RECONSTRUCTED — BYTE-IDENTICAL | LC and MC uppercase FASTAs exactly match historical outputs |
| Greengenes2 2022.10 V3–V4 reference extraction | ORIGINAL | `scripts/12_gg2_extract_v3v4.pbs`; exact match to dated 2026-08-22 GG2_COMPLETE snapshot |
| Greengenes2 taxonomy classification | ORIGINAL | Historical QIIME2 classification scripts preserved in final code snapshot |
| Kraken2 rescue / cross-check classification | ORIGINAL | Historical final Kraken2 scripts preserved in final code snapshot |
| MC + GG2 ≥97% abundance/taxonomy table construction | ORIGINAL | Final Python scripts preserved in final code snapshot |
| Historical Wild/Lab origin rule | ORIGINAL LOGIC RECOVERED | Historical bash history: BT → Lab; AG/ZX → Wild |
| Historical 576-sample metadata state | PRESERVED / TRACEABLE | `sample_metadata_location_576.tsv` is byte-identical to `sample_metadata_before_6control_fix.tsv` |
| Six-negative-control correction | ORIGINAL | `scripts/49_correct_BT_negative_controls.py`; exactly two BT H12 sample rows changed. **Partly superseded by script 53:** script 49 wrongly took BT_CLX_120 as the BT_P02 control |
| BT_P02 negative-control correction (2026-10-07) | ORIGINAL | `scripts/53_swap_BT_P02_negative_control.py`; per `01_demux_metadata/demux_metadata_audit.csv`, BT_P02 H12 = BT_ALB_192 (now control); BT_CLX_120 (well H3) restored to Lab-reared sample; backup `metadata/sample_metadata_before_BTP02_swap.tsv`; validation PASS |
| Final sample/plotting input preparation | ORIGINAL | Scripts 45, 50, 51 re-run after script 53 (code unchanged); script 52 expectations updated (controls; genus means) — pre-swap copy in `development_and_superseded/`; 52 result PASS 44 / FAIL 0 / WARN 1 (OTU_441, *Staphylococcus epidermidis*, present only in controls); final biological set = 380 Wild + 190 Lab-reared |
| Earlier four-control plotting workflow | SUPERSEDED | Scripts 46–47 and associated outputs retained only for provenance |
| Final collection-location matching | FINAL | Locations as in `sample_metadata_current.tsv`; 345 Wild samples assigned (Mandai Zoo 194, Kranji Marshes 88, Bird Paradise 39, River Safari 24); 35 without a confirmed location excluded from site analyses (`R_analysis/02_locations.R`) |
| Species metadata integration | FINAL (done in R, not on Atlas) | `R_analysis/11a_build_host_species.R` matches wild specimen IDs to COI identifications (Agabus P1: BLAST vs GBIF/BOLD; Zi Xin: tree placement) and uses `01_demux_metadata/demux_metadata_audit.csv` for plate/well; output `wild_host_species_COI.tsv` (380 rows: 219 *Ae. albopictus*, 36 *Mansonia indiana*, 23 *Culex* spp., 2 *Anopheles*, 1 *Aedeomyia*, 4 unidentified, 95 pending AG_P02). COI source spreadsheets not redistributed. Lab species from sample ID (BT_ALB / BT_CLX) |
| Final R statistics and plotting scripts | FINAL | Snapshot 2026-10-08 in `R_analysis/` (00–11, 11a); inputs = script 50/51 outputs + `sample_metadata_current.tsv` (SHA-256 834d8162…); see `R_analysis/README_R.md`. Draft all-in-one script `LSM4288_2.11_community_analyses.R` SUPERSEDED and not included |

---

## Important naming note

Historical directory/script labels reversed the later LC/MC terminology.

Final terminology:

- LC = 0.1% per-sample support threshold; less conservative; 7,642 OTUs.
- MC = max(0.5% of sample reads, 20 reads); more conservative; 867 OTUs.

Historical names such as `MC_0p1pct` correspond to current LC.
Historical names such as `LC_max0p5pct_20reads` correspond to current MC.

---

## Final sample scope

- 576 total midgut samples across six 96-well plates
- 380 Wild
- 190 Lab-reared
- 6 negative controls
- 570 biological samples for Wild/Lab analysis
- AG_Repro samples excluded from the final midgut dataset

---

## Software environment snapshot

File:
`ATLAS_SOFTWARE_VERSIONS_20261006.txt`

SHA-256:
`08dbfe1066bc17c2dde29397010ebe0794c3049b92ac8bfdfe441d2b4db675ec`

Recorded software:
- QIIME2 / q2cli 2026.7.0
- Python 3.12.13 in QIIME2 environment
- NGSpeciesID 0.3.1
- Python 3.11.16 in NGSpeciesID environment
- Racon 1.5.0
- SPOA 4.1.5
- VSEARCH 2.31.0
- Kraken2 2.17.1
- SeqKit 2.13.0

---

## Repository policy

Historical originals must remain unmodified.

Reconstructed code must be stored separately and explicitly labelled as reconstructed.

Superseded scripts should remain available for provenance but must not be presented as final analysis code.

Location, species and R analysis code were finalised on 2026-10-08 and are included in `R_analysis/`.
