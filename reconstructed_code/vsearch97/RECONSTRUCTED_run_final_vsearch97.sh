#!/usr/bin/env bash
set -euo pipefail

# RECONSTRUCTED CODE
#
# The original final VSEARCH97 commands were executed interactively and were
# not preserved as standalone scripts.
#
# These commands were reconstructed from preserved inputs, VSEARCH logs,
# analysis documentation and historical outputs.
#
# Validation performed 2026-10-06:
# - LC centroid FASTA: byte-identical to historical output
# - LC UC file:        byte-identical to historical output
# - MC centroid FASTA: byte-identical to historical output
# - MC UC file:        byte-identical to historical output
#
# FINAL terminology:
# LC = 0.1% per-sample support threshold
# MC = max(0.5% of sample reads, 20 reads)
#
# Historical directory names sometimes used the opposite LC/MC labels.

RUN="${RUN:-/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903}"
VSEARCH="${VSEARCH:-$HOME/bin/vsearch}"
OUT="${1:-$PWD/reconstructed_vsearch97_outputs}"

LC_INPUT="$RUN/04_global_vsearch97/postracon_thresholds_576/0p1pct_postracon_filtered.fasta"
MC_INPUT="$RUN/04_global_vsearch97/postracon_thresholds_576/LC_max0p5pct_20reads/LC_max0p5pct_20reads_postracon_filtered.fasta"

mkdir -p \
    "$OUT/LC_0p1pct" \
    "$OUT/MC_max0p5pct_20reads"

echo "===== LC: 0.1% per-sample threshold ====="

"$VSEARCH" \
    --cluster_fast "$LC_INPUT" \
    --id 0.97 \
    --iddef 2 \
    --strand both \
    --centroids "$OUT/LC_0p1pct/otu97_centroids.fasta" \
    --uc "$OUT/LC_0p1pct/otu97_clusters.uc" \
    --relabel OTU_ \
    --threads 1

echo
echo "===== MC: max(0.5%, 20 reads) threshold ====="

"$VSEARCH" \
    --cluster_fast "$MC_INPUT" \
    --id 0.97 \
    --iddef 2 \
    --strand both \
    --centroids "$OUT/MC_max0p5pct_20reads/otu97_centroids.fasta" \
    --uc "$OUT/MC_max0p5pct_20reads/otu97_clusters.uc" \
    --relabel OTU_ \
    --threads 1

echo
echo "VSEARCH97 reconstruction complete."
