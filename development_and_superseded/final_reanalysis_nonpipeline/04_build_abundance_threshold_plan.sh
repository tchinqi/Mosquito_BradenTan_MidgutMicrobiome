#!/usr/bin/env bash
set -euo pipefail

PROJECT="/hpctmp/e0958567/Braden_16S_20260819"
RUN="$PROJECT/20_FINAL_reanalysis_20260903"

INPUT="$RUN/01_demux_length_qc/sample_read_counts.tsv"
OUT="$RUN/02_abundance_ratio/threshold_plan.tsv"

echo "=========================================="
echo "STEP 04: ABUNDANCE THRESHOLD PLAN"
echo "=========================================="
echo "Started: $(date -Is)"
echo

awk -F'\t' 'BEGIN {
    OFS="\t"
    print "sample_id",
          "demultiplexed_reads",
          "reads_350_550",
          "cutoff_0p1pct",
          "cutoff_0p2pct",
          "cutoff_0p5pct",
          "cutoff_abs10",
          "cutoff_abs20"
}

NR>1 {
    n=$3

    c01=int(n*0.001)
    c02=int(n*0.002)
    c05=int(n*0.005)

    # Technical minimum: singleton clusters do not proceed to Racon.
    if (c01 < 2) c01=2
    if (c02 < 2) c02=2
    if (c05 < 2) c05=2

    print $1,$2,$3,c01,c02,c05,10,20
}
' "$INPUT" > "$OUT"

echo "Samples:"
awk 'END {print NR-1}' "$OUT"

echo
echo "===== EFFECTIVE CUTOFF RANGES ====="

awk -F'\t' '
NR==2 {
    for(i=4;i<=8;i++) {
        min[i]=$i
        max[i]=$i
    }
}
NR>1 {
    for(i=4;i<=8;i++) {
        if($i < min[i]) min[i]=$i
        if($i > max[i]) max[i]=$i
    }
}
END {
    print "0.1%:", min[4], "-", max[4], "reads"
    print "0.2%:", min[5], "-", max[5], "reads"
    print "0.5%:", min[6], "-", max[6], "reads"
    print "Absolute 10:", min[7], "-", max[7], "reads"
    print "Absolute 20:", min[8], "-", max[8], "reads"
}
' "$OUT"

echo
echo "Preview:"
column -t -s $'\t' "$OUT" | head -10

echo
echo "Completed: $(date -Is)"
