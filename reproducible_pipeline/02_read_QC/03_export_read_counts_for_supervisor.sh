#!/usr/bin/env bash
# =============================================================================
# Script 03 - 03_export_read_counts_for_supervisor.sh
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Read QC
# Purpose : Exports the per-well read-count table as CSV (used to set per-
#           sample abundance-ratio cut-offs).
# Input   : sample_read_counts_with_retention.tsv
# Output  : FINAL_sample_read_counts_for_abundance_ratio.csv
# Run on  : Atlas9 (bash)
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================
set -euo pipefail

PROJECT="/hpctmp/e0958567/Braden_16S_20260819"
RUN="$PROJECT/20_FINAL_reanalysis_20260903"

INPUT="$RUN/01_demux_length_qc/sample_read_counts_with_retention.tsv"
OUTPUT="$RUN/01_demux_length_qc/FINAL_sample_read_counts_for_abundance_ratio.csv"

echo "=========================================="
echo "STEP 03: EXPORT READ-COUNT TABLE"
echo "=========================================="
echo "Started: $(date -Is)"
echo

awk -F'\t' 'BEGIN {
    OFS=","
}
NR==1 {
    print "Sample ID",
          "Sequences demultiplexed",
          "Sequences passing 350-550 bp",
          "Retention (%)"
    next
}
{
    print $1,$2,$3,$4
}
' "$INPUT" > "$OUTPUT"

echo "Output:"
echo "$OUTPUT"

echo
echo -n "Samples exported: "
awk 'END {print NR-1}' "$OUTPUT"

echo
echo "Preview:"
head -6 "$OUTPUT"

echo
echo "Completed: $(date -Is)"
