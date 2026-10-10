#!/usr/bin/env bash
set -euo pipefail

RUN="/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903"

TEST="abs5"
COMMON="$RUN/manifests/common_successful_midgut_samples.txt"
BASE="$RUN/03_ngspeciesid_racon/$TEST/samples"

OUT="$RUN/04_global_vsearch97/$TEST/common570"

POOLED="$OUT/pooled_racon_polished.fasta"
PROV="$OUT/consensus_provenance.tsv"
PER_SAMPLE="$OUT/per_sample_racon_counts.tsv"
CENTROIDS="$OUT/otu97_centroids.fasta"
UC="$OUT/otu97_clusters.uc"

mkdir -p "$OUT"

: > "$POOLED"

printf \
"consensus_id\tsample_id\tcluster_id\tsupporting_reads\tsequence_length\tsource_file\n" \
> "$PROV"

printf \
"sample_id\tfinal_racon_consensuses\trepresented_supporting_reads\n" \
> "$PER_SAMPLE"

echo "============================================================"
echo "STEP 09: ABS5 — COMMON 570 — GLOBAL VSEARCH 97%"
echo "============================================================"
echo "Started: $(date -Is)"
echo "Common samples: $(wc -l < "$COMMON")"
echo

CONS_ID=0

while IFS= read -r SAMPLE
do
    [[ -n "$SAMPLE" ]] || continue

    SAMPLE_DIR="$BASE/$SAMPLE"

    SAMPLE_N=0
    SAMPLE_SUPPORT=0

    while IFS= read -r RACON
    do
        [[ -s "$RACON" ]] || continue

        CLUSTER_DIR=$(basename "$(dirname "$RACON")")
        CID="${CLUSTER_DIR#racon_cl_id_}"

        READS="$SAMPLE_DIR/reads_to_consensus_${CID}.fastq"

        if [[ -s "$READS" ]]; then
            SUPPORT=$(awk 'END {print NR/4}' "$READS")
        else
            SUPPORT=0
        fi

        SEQ=$(awk '
            /^>/ {next}
            {printf "%s",$0}
            END {printf "\n"}
        ' "$RACON")

        [[ -n "$SEQ" ]] || continue

        CONS_ID=$((CONS_ID + 1))
        SAMPLE_N=$((SAMPLE_N + 1))
        SAMPLE_SUPPORT=$((SAMPLE_SUPPORT + SUPPORT))

        ID=$(printf "%s_CONS_%06d" "$TEST" "$CONS_ID")
        LEN=${#SEQ}

        printf ">%s\n%s\n" "$ID" "$SEQ" >> "$POOLED"

        printf \
        "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$ID" "$SAMPLE" "$CID" "$SUPPORT" "$LEN" "$RACON" \
        >> "$PROV"

    done < <(
        find "$SAMPLE_DIR" \
            -mindepth 2 \
            -maxdepth 2 \
            -type f \
            -path '*/racon_cl_id_*/racon_polished_it_2.fasta' \
            | sort
    )

    printf "%s\t%s\t%s\n" \
        "$SAMPLE" "$SAMPLE_N" "$SAMPLE_SUPPORT" \
        >> "$PER_SAMPLE"

done < "$COMMON"

POOLED_N=$(grep -c '^>' "$POOLED")

echo "Pooled Racon-polished consensuses: $POOLED_N"

echo
echo "Running VSEARCH 97% clustering..."

"$HOME/bin/vsearch" \
    --cluster_fast "$POOLED" \
    --id 0.97 \
    --iddef 2 \
    --strand both \
    --centroids "$CENTROIDS" \
    --uc "$UC" \
    --relabel OTU_ \
    --threads 1

OTU_N=$(grep -c '^>' "$CENTROIDS")

echo
echo "============================================================"
echo "ABS5 COMMON-570 RESULT"
echo "============================================================"
echo "Common samples: 570"
echo "Pooled Racon consensuses: $POOLED_N"
echo "97% OTUs: $OTU_N"
echo "Completed: $(date -Is)"
