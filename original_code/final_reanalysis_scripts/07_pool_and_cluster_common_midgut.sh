#!/usr/bin/env bash
set -euo pipefail

PROJECT="/hpctmp/e0958567/Braden_16S_20260819"
RUN="$PROJECT/20_FINAL_reanalysis_20260903"

COMMON="$RUN/manifests/common_successful_midgut_samples.txt"
MASTER_SUMMARY="$RUN/04_global_vsearch97/common570_threshold_comparison.tsv"

VSEARCH="$HOME/bin/vsearch"

mkdir -p "$RUN/04_global_vsearch97"

printf \
"condition\tcommon_samples\tpooled_racon_consensuses\trepresented_supporting_reads\totu97_count\n" \
> "$MASTER_SUMMARY"

echo "============================================================"
echo "STEP 07: COMMON-570 RACON POOLING + VSEARCH 97%"
echo "============================================================"
echo "Started: $(date -Is)"
echo
echo "Common samples: $(wc -l < "$COMMON")"
echo "VSEARCH: $($VSEARCH --version 2>&1 | head -1)"
echo

for T in 0p1pct 0p2pct 0p5pct abs10 abs20
do
    echo
    echo "############################################################"
    echo "CONDITION: $T"
    echo "############################################################"

    RBASE="$RUN/03_ngspeciesid_racon/$T"
    SSUM="$RBASE/sample_summary.tsv"

    OUT="$RUN/04_global_vsearch97/$T/common570"
    mkdir -p "$OUT"

    POOLED="$OUT/racon_it2_common570.fasta"
    PROV="$OUT/consensus_provenance.tsv"
    PER_SAMPLE="$OUT/per_sample_racon_counts.tsv"

    CENTROIDS="$OUT/otu97_centroids.fasta"
    UC="$OUT/otu97_clusters.uc"

    : > "$POOLED"

    printf \
    "consensus_id\tsample_id\tcluster_id\tsupporting_reads\tsequence_length\tsource_file\n" \
    > "$PROV"

    printf \
    "sample_id\tfinal_racon_consensuses\trepresented_supporting_reads\n" \
    > "$PER_SAMPLE"

    GLOBAL_N=0
    GLOBAL_SUPPORT=0
    SAMPLE_N=0

    while IFS= read -r SAMPLE
    do
        [[ -n "$SAMPLE" ]] || continue

        SAMPLE_N=$((SAMPLE_N+1))

        STATUS=$(awk -F'\t' -v s="$SAMPLE" '
            NR>1 && $1==s {print $9; exit}
        ' "$SSUM")

        if [[ "$STATUS" != "SUCCESS" ]]; then
            echo "ERROR: $SAMPLE is not SUCCESS for $T" >&2
            exit 1
        fi

        DIR="$RBASE/samples/$SAMPLE"

        SAMPLE_CONS=0
        SAMPLE_SUPPORT=0

        for RDIR in "$DIR"/racon_cl_id_*
        do
            [[ -d "$RDIR" ]] || continue

            CLUSTER=$(basename "$RDIR")
            CLID=${CLUSTER#racon_cl_id_}

            FINAL="$RDIR/racon_polished_it_2.fasta"
            READS="$DIR/reads_to_consensus_${CLID}.fastq"

            # Carry forward only a non-empty final Racon iteration.
            if [[ ! -s "$FINAL" ]] || ! grep -q '^>' "$FINAL"; then
                echo "ERROR: invalid final Racon output:" >&2
                echo "$FINAL" >&2
                exit 1
            fi

            if [[ ! -s "$READS" ]]; then
                echo "ERROR: supporting-read FASTQ missing:" >&2
                echo "$READS" >&2
                exit 1
            fi

            SUPPORT=$(awk 'END {print NR/4}' "$READS")

            SEQ=$(awk '
                /^>/ {next}
                {
                    gsub(/[[:space:]]/, "")
                    seq=seq $0
                }
                END {print seq}
            ' "$FINAL")

            LEN=${#SEQ}

            if [[ "$LEN" -eq 0 ]]; then
                echo "ERROR: empty sequence in $FINAL" >&2
                exit 1
            fi

            GLOBAL_N=$((GLOBAL_N+1))
            SAMPLE_CONS=$((SAMPLE_CONS+1))

            GLOBAL_SUPPORT=$((GLOBAL_SUPPORT+SUPPORT))
            SAMPLE_SUPPORT=$((SAMPLE_SUPPORT+SUPPORT))

            CID=$(printf "%s_CONS_%06d" "$T" "$GLOBAL_N")

            printf ">%s\n%s\n" \
                "$CID" "$SEQ" >> "$POOLED"

            printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$CID" \
                "$SAMPLE" \
                "$CLUSTER" \
                "$SUPPORT" \
                "$LEN" \
                "$FINAL" \
                >> "$PROV"

        done

        printf "%s\t%s\t%s\n" \
            "$SAMPLE" \
            "$SAMPLE_CONS" \
            "$SAMPLE_SUPPORT" \
            >> "$PER_SAMPLE"

    done < "$COMMON"

    echo "Common wells processed: $SAMPLE_N"
    echo "Final Racon consensuses pooled: $GLOBAL_N"
    echo "Represented supporting reads: $GLOBAL_SUPPORT"

    FASTA_N=$(grep -c '^>' "$POOLED")

    if [[ "$FASTA_N" -ne "$GLOBAL_N" ]]; then
        echo "ERROR: pooled FASTA count mismatch." >&2
        exit 1
    fi

    echo
    echo "Running VSEARCH 97% clustering..."

    "$VSEARCH" \
        --cluster_fast "$POOLED" \
        --id 0.97 \
        --iddef 2 \
        --strand both \
        --centroids "$CENTROIDS" \
        --uc "$UC" \
        --relabel OTU_ \
        --threads 1

    OTUS=$(grep -c '^>' "$CENTROIDS")

    echo
    echo "97% OTUs: $OTUS"

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$T" \
        "$SAMPLE_N" \
        "$GLOBAL_N" \
        "$GLOBAL_SUPPORT" \
        "$OTUS" \
        >> "$MASTER_SUMMARY"

done

echo
echo "============================================================"
echo "FINAL COMMON-570 THRESHOLD COMPARISON"
echo "============================================================"

column -t -s $'\t' "$MASTER_SUMMARY"

echo
echo "Completed: $(date -Is)"
