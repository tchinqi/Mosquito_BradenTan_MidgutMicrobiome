#!/usr/bin/env bash
set -euo pipefail

PROJECT="/hpctmp/e0958567/Braden_16S_20260819"
RUN="$PROJECT/20_FINAL_reanalysis_20260903"

MID="$RUN/manifests/all_576_midgut_samples.txt"
UPSTREAM="$RUN/01_demux_length_qc/sample_read_counts.tsv"

OUTDIR="$RUN/02_abundance_ratio/final_per_sample_comparison"
OUT="$OUTDIR/per_sample_threshold_counts_576.tsv"

mkdir -p "$OUTDIR"

SEQKIT=$(command -v seqkit)

if [[ -z "$SEQKIT" ]]; then
    echo "ERROR: seqkit not found in PATH" >&2
    exit 1
fi

echo "============================================================"
echo "STEP 08: PER-SAMPLE ABUNDANCE-THRESHOLD COUNTS"
echo "============================================================"
echo "Started: $(date -Is)"
echo "SeqKit: $($SEQKIT version 2>&1 | head -1)"
echo

printf \
"sample_id\tdemultiplexed_reads\treads_350_550\tracon_0p1pct\tracon_0p2pct\tracon_0p5pct\tracon_abs10\tracon_abs20\tstatus_0p1pct\tstatus_0p2pct\tstatus_0p5pct\tstatus_abs10\tstatus_abs20\n" \
> "$OUT"

count_racon_sequences () {

    local TEST="$1"
    local SAMPLE="$2"

    local DIR="$RUN/03_ngspeciesid_racon/$TEST/samples/$SAMPLE"
    local TMP

    TMP=$(mktemp)

    find "$DIR" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -path '*/racon_cl_id_*/racon_polished_it_2.fasta' \
        -size +0c \
        -exec cat {} \; > "$TMP"

    if [[ ! -s "$TMP" ]]; then
        echo 0
        rm -f "$TMP"
        return
    fi

    "$SEQKIT" stats -T "$TMP" 2>/dev/null \
        | awk -F'\t' 'NR==2 {print $4+0}'

    rm -f "$TMP"
}

while IFS= read -r SAMPLE
do
    [[ -n "$SAMPLE" ]] || continue

    DEMUX=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $2; exit}' "$UPSTREAM")

    LEN350=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $3; exit}' "$UPSTREAM")

    C01=$(count_racon_sequences "0p1pct" "$SAMPLE")
    C02=$(count_racon_sequences "0p2pct" "$SAMPLE")
    C05=$(count_racon_sequences "0p5pct" "$SAMPLE")
    C10=$(count_racon_sequences "abs10" "$SAMPLE")
    C20=$(count_racon_sequences "abs20" "$SAMPLE")

    S01=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $9; exit}' \
        "$RUN/03_ngspeciesid_racon/0p1pct/sample_summary.tsv")

    S02=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $9; exit}' \
        "$RUN/03_ngspeciesid_racon/0p2pct/sample_summary.tsv")

    S05=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $9; exit}' \
        "$RUN/03_ngspeciesid_racon/0p5pct/sample_summary.tsv")

    S10=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $9; exit}' \
        "$RUN/03_ngspeciesid_racon/abs10/sample_summary.tsv")

    S20=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $9; exit}' \
        "$RUN/03_ngspeciesid_racon/abs20/sample_summary.tsv")

    printf \
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" \
    "$DEMUX" \
    "$LEN350" \
    "$C01" \
    "$C02" \
    "$C05" \
    "$C10" \
    "$C20" \
    "$S01" \
    "$S02" \
    "$S05" \
    "$S10" \
    "$S20" \
    >> "$OUT"

done < "$MID"

echo
echo "===== QC ====="

echo -n "Midgut wells: "
awk 'END {print NR-1}' "$OUT"

echo
echo "===== COLUMN TOTALS ====="

awk -F'\t' '
NR>1 {
    demux += $2
    len   += $3
    p01   += $4
    p02   += $5
    p05   += $6
    a10   += $7
    a20   += $8
}
END {
    print "Demultiplexed reads:", demux
    print "350-550 reads:", len
    print "0.1% Racon-polished sequences:", p01
    print "0.2% Racon-polished sequences:", p02
    print "0.5% Racon-polished sequences:", p05
    print "10-read Racon-polished sequences:", a10
    print "20-read Racon-polished sequences:", a20
}
' "$OUT"

echo
echo "Completed: $(date -Is)"
