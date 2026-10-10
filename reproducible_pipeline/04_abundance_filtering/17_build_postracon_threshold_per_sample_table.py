#!/usr/bin/env python3
# =============================================================================
# Script 17 - 17_build_postracon_threshold_per_sample_table.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Abundance-support filtering
# Purpose : Per-sample table of consensuses and supporting reads retained
#           under each threshold.
# Input   : pooled provenance; threshold_plan.tsv; all_576_midgut_samples.txt
# Output  : Priya_postracon_thresholds_per_sample_576.tsv
# Run on  : Atlas9 (Python 3)
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================

import csv
from pathlib import Path
from collections import defaultdict

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")

PROV = RUN / "04_global_vsearch97/min2_complete576/pooled_racon_polished_576_provenance.tsv"
PLAN = RUN / "02_abundance_ratio/threshold_plan.tsv"
MAN  = RUN / "manifests/all_576_midgut_samples.txt"

OUT = RUN / "04_global_vsearch97/postracon_thresholds_576/Priya_postracon_thresholds_per_sample_576.tsv"

# ------------------------------------------------------------
# Read exact 576-sample order
# ------------------------------------------------------------
with open(MAN) as fh:
    samples = [x.strip() for x in fh if x.strip()]

if len(samples) != 576:
    raise RuntimeError("Expected 576 samples, found {}".format(len(samples)))

# ------------------------------------------------------------
# Read sample-specific percentage cutoffs
# ------------------------------------------------------------
cutoffs = {}

with open(PLAN, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        cutoffs[row["sample_id"]] = {
            "0p1pct": int(row["cutoff_0p1pct"]),
            "0p2pct": int(row["cutoff_0p2pct"]),
            "0p5pct": int(row["cutoff_0p5pct"]),
            "abs5": 5,
            "abs10": 10,
            "abs20": 20,
        }

# ------------------------------------------------------------
# Count retained FINAL Racon consensuses per sample
# ------------------------------------------------------------
conditions = ["0p1pct", "0p2pct", "0p5pct", "abs5", "abs10", "abs20"]

counts = {
    sample: {condition: 0 for condition in conditions}
    for sample in samples
}

racon_total = defaultdict(int)

with open(PROV, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    n_records = 0

    for row in reader:
        n_records += 1

        sample = row["sample_id"]
        support = int(row["supporting_reads"])

        if sample not in counts:
            raise RuntimeError(
                "Provenance sample not in final 576 manifest: {}".format(sample)
            )

        if sample not in cutoffs:
            raise RuntimeError(
                "Sample missing threshold plan: {}".format(sample)
            )

        racon_total[sample] += 1

        for condition in conditions:
            if support >= cutoffs[sample][condition]:
                counts[sample][condition] += 1

if n_records != 22528:
    raise RuntimeError(
        "Expected 22528 Racon consensuses, found {}".format(n_records)
    )

# ------------------------------------------------------------
# Write per-sample table
# ------------------------------------------------------------
with open(OUT, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")

    writer.writerow([
        "sample_id",
        "racon_polishing_only",
        "0.1%_filter",
        "0.2%_filter",
        "0.5%_filter",
        "5_reads",
        "10_reads",
        "20_reads"
    ])

    for sample in samples:
        writer.writerow([
            sample,
            racon_total[sample],
            counts[sample]["0p1pct"],
            counts[sample]["0p2pct"],
            counts[sample]["0p5pct"],
            counts[sample]["abs5"],
            counts[sample]["abs10"],
            counts[sample]["abs20"]
        ])

print("===== FINAL PER-SAMPLE POST-RACON THRESHOLD TABLE =====")
print("Samples:", len(samples))

print("Racon-only total:",
      sum(racon_total[s] for s in samples))

for c in conditions:
    print("{} total: {}".format(
        c,
        sum(counts[s][c] for s in samples)
    ))

print()
print("Output:")
print(OUT)
