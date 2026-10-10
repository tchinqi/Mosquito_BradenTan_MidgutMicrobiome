#!/usr/bin/env bash
set -euo pipefail

PROJECT="/hpctmp/e0958567/Braden_16S_20260819"
RUN="$PROJECT/20_FINAL_reanalysis_20260903"

MID="$RUN/manifests/all_576_midgut_samples.txt"
OUT="$RUN/manifests/common_successful_midgut_samples.txt"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "============================================================"
echo "STEP 06: COMMON SUCCESSFUL MIDGUT SAMPLE SET"
echo "============================================================"
echo "Started: $(date -Is)"
echo

for T in 0p1pct 0p2pct 0p5pct abs10 abs20
do
    S="$RUN/03_ngspeciesid_racon/$T/sample_summary.tsv"

    awk -F'\t' '
        NR>1 && $9=="SUCCESS" {print $1}
    ' "$S" | sort > "$TMP/${T}_success.txt"

    echo "----- $T -----"
    echo -n "Successful midgut wells: "

    awk '
        NR==FNR {mid[$1]=1; next}
        ($1 in mid) {n++}
        END {print n+0}
    ' "$MID" "$TMP/${T}_success.txt"

    echo -n "Failed midgut wells: "

    awk '
        NR==FNR {success[$1]=1; next}
        !($1 in success) {print}
    ' "$TMP/${T}_success.txt" "$MID" \
    | tee "$TMP/${T}_failed.txt" \
    | wc -l

    if [[ -s "$TMP/${T}_failed.txt" ]]; then
        echo "Failed IDs:"
        cat "$TMP/${T}_failed.txt"
    fi

    echo
done

# A sample is retained only if SUCCESS under all five conditions.
awk '
    FILENAME ~ /0p1pct_success/ {a[$1]=1; next}
    FILENAME ~ /0p2pct_success/ {b[$1]=1; next}
    FILENAME ~ /0p5pct_success/ {c[$1]=1; next}
    FILENAME ~ /abs10_success/  {d[$1]=1; next}
    FILENAME ~ /abs20_success/  {e[$1]=1; next}
    END {
        for (x in a)
            if (x in b && x in c && x in d && x in e)
                print x
    }
' \
"$TMP/0p1pct_success.txt" \
"$TMP/0p2pct_success.txt" \
"$TMP/0p5pct_success.txt" \
"$TMP/abs10_success.txt" \
"$TMP/abs20_success.txt" \
| sort > "$TMP/common_all.txt"

# Restrict explicitly to the midgut manifest.
awk '
    NR==FNR {common[$1]=1; next}
    ($1 in common) {print}
' "$TMP/common_all.txt" "$MID" > "$OUT"

echo "============================================================"
echo -n "Original midgut wells: "
wc -l < "$MID"

echo -n "Common successful wells across all five filters: "
wc -l < "$OUT"

echo
echo "Completed: $(date -Is)"
