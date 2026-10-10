# Historical final-reanalysis scripts not promoted to the final pipeline

These 32 scripts are preserved for provenance but are not part of the
curated final reproducible pipeline.

## Development / parameter exploration / testing

- 04_build_abundance_threshold_plan.sh
- 05_ngspeciesid_racon_threshold.pbs
- 05b_ngspeciesid_racon_abs5_midgut.pbs
- 05b_ngspeciesid_racon_abs5_midgut_4gb_backup.pbs
- 11_ngspeciesid_racon_noabund_midgut.pbs
- 13_test_rescue_one_sample.pbs
- 15_rescue_0p1_0p2_thresholds.pbs
- 28_kraken2_sanity_test.pbs
- 35_kraken2_full_LC_MC_conf80_conf97.pbs

## Intermediate / supporting analyses

- 06_build_common_midgut_manifest.sh
- 07_pool_and_cluster_common_midgut.sh
- 08_build_per_sample_threshold_table.sh
- 09_pool_and_cluster_abs5_common570.sh
- 10_build_final_priya_midgut_table.sh
- 21_make_full_GG2_taxon_tables.py
- 21_summarise_GG2_extract_unassigned_MC_LC.py
- 22_GG2_species_genus_summary_MC_LC.py
- 23_make_MC_LC_genus_species_comparison.py
- 25_GG2_family_97pct_corrected_MC_LC.py
- 26_GG2_family_genus_species_97pct_corrected.py
- 29_kraken2_final_MC_LC.pbs
- 31_make_Priya_GG2_hierarchical_tables.py
- 32_make_Priya_Kraken_hierarchical_tables.py
- 34_kraken2_GG2_default80_unassigned.pbs

## Superseded / replaced by later final workflow

- 19_greengenes2_MC_LC_final.pbs
  Replaced by the uppercase Greengenes2 workflow.

- 27_GG2_complete_80_vs_97_tables.py
  Replaced by later corrected GG2 comparison/table scripts.

- 30_make_final_GG2_Kraken_combined_tables.py
  Earlier GG2/Kraken aggregation workflow; replaced by later final rescue workflow.

- 33_extract_GG2_default80_unassigned_for_Kraken.py
  Earlier GG80 rescue extraction; final workflow also accounts for domain-only assignments.

- 38_kraken2_GG2_default80_unassigned_rescue.pbs
- 39_summarise_GG2_unassigned_Kraken_rescue.py
  Earlier Kraken rescue workflow; replaced by scripts 41–43.

- 46_prepare_MC_GG97_plotting_tables.py
- 47_summarise_MC_GG97_wild_vs_lab_genera.py
  Superseded four-negative-control plotting workflow.
  Final workflow uses six negative controls via scripts 49–51.

## Important note

Files in this directory are retained to preserve the development history.
Their presence does not mean they should be rerun as part of the final
analysis.

The authoritative upstream workflow is the curated `reproducible_pipeline/`.
