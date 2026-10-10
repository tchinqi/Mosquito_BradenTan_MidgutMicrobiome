#!/usr/bin/env bash
set -euo pipefail

RUN="/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903"

MID="$RUN/manifests/all_576_midgut_samples.txt"
UPSTREAM="$RUN/01_demux_length_qc/sample_read_counts.tsv"

OUTDIR="$RUN/02_abundance_ratio/final_per_sample_comparison"
DATA="$OUTDIR/PRIYA_midgut_576_per_sample_six_thresholds.tsv"
FINAL="$OUTDIR/PRIYA_midgut_576_FINAL_for_excel.tsv"

mkdir -p "$OUTDIR"

SEQKIT=$(command -v seqkit)

if [[ -z "$SEQKIT" ]]; then
    echo "ERROR: seqkit not found in PATH" >&2
    exit 1
fi

count_racon () {

    local TEST="$1"
    local SAMPLE="$2"
    local DIR="$RUN/03_ngspeciesid_racon/$TEST/samples/$SAMPLE"

    if [[ ! -d "$DIR" ]]; then
        echo 0
        return
    fi

    mapfile -t FILES < <(
        find "$DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -path '*/racon_cl_id_*/racon_polished_it_2.fasta' \
            -size +0c \
            | sort
    )

    if [[ "${#FILES[@]}" -eq 0 ]]; then
        echo 0
        return
    fi

    cat "${FILES[@]}" \
    | "$SEQKIT" stats -T - 2>/dev/null \
    | awk -F'\t' 'NR==2 {print $4+0}'
}

echo "============================================================"
echo "STEP 10: FINAL PRIYA 576-MIDGUT TABLE"
echo "============================================================"
echo "Started: $(date -Is)"
echo "SeqKit: $($SEQKIT version 2>&1 | head -1)"
echo

printf \
"Sample ID\tSequences demultiplexed\tSequences passing 350-550 bp\t0.1%% filter\t0.2%% filter\t0.5%% filter\t5 reads\t10 reads\t20 reads\n" \
> "$DATA"

N=0

while IFS= read -r SAMPLE
do
    [[ -n "$SAMPLE" ]] || continue

    N=$((N + 1))

    DEMUX=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $2; exit}' "$UPSTREAM")

    LEN=$(awk -F'\t' -v s="$SAMPLE" \
        'NR>1 && $1==s {print $3; exit}' "$UPSTREAM")

    P01=$(count_racon 0p1pct "$SAMPLE")
    P02=$(count_racon 0p2pct "$SAMPLE")
    P05=$(count_racon 0p5pct "$SAMPLE")
    A5=$(count_racon abs5 "$SAMPLE")
    A10=$(count_racon abs10 "$SAMPLE")
    A20=$(count_racon abs20 "$SAMPLE")

    printf \
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" "$DEMUX" "$LEN" \
    "$P01" "$P02" "$P05" "$A5" "$A10" "$A20" \
    >> "$DATA"

    if (( N % 50 == 0 )); then
        echo "Processed $N / 576 samples"
    fi

done < "$MID"

# Calculate totals from the actual 576 per-sample rows.
read -r DEMUX_T LEN_T P01_T P02_T P05_T A5_T A10_T A20_T < <(
    awk -F'\t' '
    NR>1 {
        b += $2
        c += $3
        d += $4
        e += $5
        f += $6
        g += $7
        h += $8
        i += $9
    }
    END {
        print b,c,d,e,f,g,h,i
    }
    ' "$DATA"
)

# Build table in the same general layout as Priya's sheet.
{
    printf \
    "OTU count (97%% clustering; common 570 wells)\t-\t-\t7533\t3894\t1434\t2781\t1431\t898\n"

    printf \
    "Total sequences after each step (all 576 wells)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$DEMUX_T" "$LEN_T" \
    "$P01_T" "$P02_T" "$P05_T" "$A5_T" "$A10_T" "$A20_T"

    cat "$DATA"

} > "$FINAL"

echo
echo "===== FINAL QC ====="

echo -n "Midgut sample rows: "
awk 'END {print NR-1}' "$DATA"

echo
echo "Totals:"
printf "Demultiplexed: %s\n" "$DEMUX_T"
printf "350-550 bp:    %s\n" "$LEN_T"
printf "0.1%%:          %s\n" "$P01_T"
printf "0.2%%:          %s\n" "$P02_T"
printf "0.5%%:          %s\n" "$P05_T"
printf "5 reads:       %s\n" "$A5_T"
printf "10 reads:      %s\n" "$A10_T"
printf "20 reads:      %s\n" "$A20_T"

echo
echo "Final Excel-ready table:"
echo "$FINAL"

echo
echo "Completed: $(date -Is)"
