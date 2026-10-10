#!/usr/bin/env bash
# =============================================================================
# Script 02 - 02_summarise_read_depths.sh
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Read QC
# Purpose : Adds per-well retention (% of reads kept by the length filter) and
#           writes a read-depth summary.
# Input   : sample_read_counts.tsv
# Output  : sample_read_counts_with_retention.tsv; read_depth_summary.tsv
# Run on  : Atlas9 (bash)
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================
set -euo pipefail

PROJECT="/hpctmp/e0958567/Braden_16S_20260819"
RUN="$PROJECT/20_FINAL_reanalysis_20260903"

INPUT="$RUN/01_demux_length_qc/sample_read_counts.tsv"
OUTDIR="$RUN/01_demux_length_qc"
SUMMARY="$OUTDIR/read_depth_summary.tsv"
PER_SAMPLE="$OUTDIR/sample_read_counts_with_retention.tsv"

echo "=========================================="
echo "STEP 02: SUMMARISE PER-SAMPLE READ DEPTH"
echo "=========================================="
echo "Started: $(date -Is)"
echo

awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {
    print $0, "retention_percent"
    next
}
{
    if ($2>0)
        retention=100*$3/$2
    else
        retention=0

    printf "%s\t%s\t%s\t%.2f\n", $1,$2,$3,retention
}
' "$INPUT" > "$PER_SAMPLE"

python3 - "$INPUT" "$SUMMARY" <<'PY'
import csv
import sys
import statistics

inp = sys.argv[1]
out = sys.argv[2]

demux = []
filtered = []

with open(inp) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        demux.append(int(row["demultiplexed_reads"]))
        filtered.append(int(row["reads_350_550"]))

def q(values, p):
    x = sorted(values)
    pos = (len(x)-1)*p
    lo = int(pos)
    hi = min(lo+1, len(x)-1)
    frac = pos-lo
    return x[lo]*(1-frac)+x[hi]*frac

with open(out, "w") as f:
    f.write("metric\tdemultiplexed_reads\treads_350_550\n")
    f.write(f"minimum\t{min(demux)}\t{min(filtered)}\n")
    f.write(f"Q1\t{q(demux,0.25):.1f}\t{q(filtered,0.25):.1f}\n")
    f.write(f"median\t{statistics.median(demux):.1f}\t{statistics.median(filtered):.1f}\n")
    f.write(f"mean\t{statistics.mean(demux):.1f}\t{statistics.mean(filtered):.1f}\n")
    f.write(f"Q3\t{q(demux,0.75):.1f}\t{q(filtered,0.75):.1f}\n")
    f.write(f"maximum\t{max(demux)}\t{max(filtered)}\n")

print("Samples:", len(demux))
print("Summary written to:", out)
PY

echo
echo "===== READ DEPTH SUMMARY ====="
column -t -s $'\t' "$SUMMARY"

echo
echo "Completed: $(date -Is)"
