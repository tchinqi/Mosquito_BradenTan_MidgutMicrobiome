# LSM4288 FYP – R analysis scripts (final, 10 Oct 2026)

Run from the R project root, in this order. All scripts `source("00_setup.R")`; seed 4288.

| Script | Report item |
|---|---|
| 00_setup.R | paths, input checks (570 bio samples, 1,422,797 reads), helpers |
| 05_classification.R | Figure 1 – Greengenes2 vs Kraken2 |
| 01_wild_vs_lab.R | Figure 2 – wild vs lab |
| 07_core_microbiome.R | Figure 3 – core microbiome |
| 02_locations.R | Figure 4 – wild collection sites (file: Fig5_locations) |
| 06_lab_species.R | Lab *Ae. albopictus* vs lab *Culex* (file: Fig_lab_species) |
| 11a_build_host_species.R | builds `data/raw/species/wild_host_species_COI.tsv` (wild host species from COI) |
| 11_host_species.R | Figure 6 – origin × species; Figure 7 – wild host species |
| 03, 04, 08, 09, 10 | supplementary figures/tables (driver genera, *Wolbachia*, controls, depth vs richness, heatmap, contamination sensitivity/source) |
| 12_discussion_checks.R | Discussion checks: *Wolbachia* excluding *Elizabethkingia* (Fig. S5), depth-adjusted and rarefied richness (Table S4B), source of negative-control reads by well and tag (Table S3) |

## Host species tagging (not done on ATLAS9)
ATLAS9 produced the OTU table and sample metadata only. Wild host species were assigned
afterwards in R (11a) by matching the specimen ID of each wild 16S sample to COI
identifications: Agabus P1 (GBIF/BOLD BLAST, 99.4–100% identity) and Zi Xin (tree
placement). Specimen IDs are normalised (uppercase, non-alphanumerics removed). The
ATLAS9 demultiplexing audit (`demux_metadata_audit.csv`) supplies plate/well. Result:
219 *Ae. albopictus*, 36 *Mansonia indiana*, 23 *Culex* spp., 2 *Anopheles*, 1 *Aedeomyia*,
4 unidentified, 95 pending (AG_P02). Lab species come from the sample ID (BT_ALB / BT_CLX).

Superseded: LSM4288_2.11_community_analyses.R (single-file draft; replaced by 01–11).

## Changes since the 08 Oct 2026 release
- 01: adds Figure S3 (mean top-8 bars + per-mosquito profiles; `FigS3_genus_composition`).
- 07: Fig. 3A labels show one decimal place where a non-core value would round to 90 (e.g. 89.7), so labels never contradict the ≥90% outline.
- 12: new script (see table above).
