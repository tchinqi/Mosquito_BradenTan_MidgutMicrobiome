#!/usr/bin/env bash
set -euo pipefail

PROJECT="/hpctmp/e0958567/Braden_16S_20260819"
RUN="$PROJECT/20_FINAL_reanalysis_20260903"

DEMUX="$RUN/inputs/02_demultiplexed"
FILTERED="$RUN/inputs/03b_length_filtered_fastq"
OUTDIR="$RUN/01_demux_length_qc"

OUT="$OUTDIR/sample_read_counts.tsv"

mkdir -p "$OUTDIR"

echo "=========================================="
echo "STEP 01: DEMULTIPLEXING / LENGTH QC"
echo "=========================================="
echo "Started: $(date -Is)"
echo
echo "Demultiplexed input:"
echo "$DEMUX"
echo
echo "350-550 bp filtered input:"
echo "$FILTERED"
echo
echo "Output:"
echo "$OUT"
echo

printf "sample_id\tdemultiplexed_reads\treads_350_550\n" > "$OUT"

for F in "$DEMUX"/*.fastq
do
    BASE=$(basename "$F")

    # Minibar non-sample outputs are not biological/sample wells.
    [[ "$BASE" == "sample_unk.fastq" ]] && continue
    [[ "$BASE" == "sample_Multiple_Matches.fastq" ]] && continue

    SAMPLE=${BASE%.fastq}
    FILTERED_FILE="$FILTERED/${SAMPLE}_filtered.fastq"

    DEMUX_N=$(awk 'END {print NR/4}' "$F")

    if [[ ! -f "$FILTERED_FILE" ]]; then
        echo "ERROR: filtered FASTQ not found for $SAMPLE" >&2
        exit 1
    fi

    FILTERED_N=$(awk 'END {print NR/4}' "$FILTERED_FILE")

    printf "%s\t%s\t%s\n" \
        "$SAMPLE" \
        "$DEMUX_N" \
        "$FILTERED_N" \
        >> "$OUT"
done

echo "=========================================="
echo "QC"
echo "=========================================="

N_SAMPLES=$(awk 'NR>1 {n++} END {print n+0}' "$OUT")
N_MISSING=$(awk -F'\t' 'NR>1 && $3=="NA" {n++} END {print n+0}' "$OUT")
N_INVALID=$(awk -F'\t' 'NR>1 && $3>$2 {n++} END {print n+0}' "$OUT")

echo "Samples represented: $N_SAMPLES"
echo "Missing filtered counts: $N_MISSING"
echo "Filtered count > demultiplexed count: $N_INVALID"

echo
echo "Total demultiplexed reads:"
awk -F'\t' 'NR>1 {s+=$2} END {print s}' "$OUT"

echo "Total reads passing 350-550 bp:"
awk -F'\t' 'NR>1 {s+=$3} END {print s}' "$OUT"

echo
echo "Completed: $(date -Is)"
