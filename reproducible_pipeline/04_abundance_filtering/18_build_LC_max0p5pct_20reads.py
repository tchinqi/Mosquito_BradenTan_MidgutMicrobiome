#!/usr/bin/env python3
# =============================================================================
# Script 18 - 18_build_LC_max0p5pct_20reads.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Abundance-support filtering
# Purpose : Builds the final MC dataset: consensuses supported by >= max(0.5%
#           of the sample's reads, 20 reads).
# Input   : pooled FASTA + provenance; threshold_plan.tsv
# Output  : LC_max0p5pct_20reads/ (FASTA, kept IDs, per-sample table)
# Run on  : Atlas9 (Python 3)
# Note    : Historical file/folder names swap LC and MC: '0p1pct' /
#           'MC_0p1pct' = final LC (0.1%); 'LC_max0p5pct_20reads' = final MC
#           (max(0.5%, 20 reads)).
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================

import csv
from pathlib import Path
from collections import defaultdict

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")

PROV = RUN / "04_global_vsearch97/min2_complete576/pooled_racon_polished_576_provenance.tsv"
FASTA = RUN / "04_global_vsearch97/min2_complete576/pooled_racon_polished_576.fasta"
PLAN = RUN / "02_abundance_ratio/threshold_plan.tsv"
MAN = RUN / "manifests/all_576_midgut_samples.txt"

OUTDIR = RUN / "04_global_vsearch97/postracon_thresholds_576/LC_max0p5pct_20reads"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_FASTA = OUTDIR / "LC_max0p5pct_20reads_postracon_filtered.fasta"
OUT_IDS = OUTDIR / "LC_max0p5pct_20reads_keep_ids.txt"
OUT_SAMPLE = OUTDIR / "LC_max0p5pct_20reads_per_sample.tsv"

# Read final 576 samples
with open(MAN) as fh:
    samples = [x.strip() for x in fh if x.strip()]

if len(samples) != 576:
    raise RuntimeError("Expected 576 samples, found {}".format(len(samples)))

# Read 0.5% cutoff and build LC cutoff = max(0.5%, 20 reads)
lc_cutoff = {}
cutoff_source = {"20_reads": 0, "0.5_percent": 0, "equal": 0}
plan_rows = {}

with open(PLAN, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        sample = row["sample_id"]
        plan_rows[sample] = row

        c05 = int(row["cutoff_0p5pct"])
        lc = max(c05, 20)

        lc_cutoff[sample] = lc

        if c05 < 20:
            cutoff_source["20_reads"] += 1
        elif c05 > 20:
            cutoff_source["0.5_percent"] += 1
        else:
            cutoff_source["equal"] += 1

# Select final Racon consensuses passing LC
keep_ids = set()
racon_total = defaultdict(int)
lc_total = defaultdict(int)

with open(PROV, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    n = 0

    for row in reader:
        n += 1

        seq_id = row["sequence_id"]
        sample = row["sample_id"]
        support = int(row["supporting_reads"])

        if sample not in lc_cutoff:
            raise RuntimeError("Missing LC cutoff for {}".format(sample))

        racon_total[sample] += 1

        if support >= lc_cutoff[sample]:
            keep_ids.add(seq_id)
            lc_total[sample] += 1

if n != 22528:
    raise RuntimeError(
        "Expected 22528 final Racon consensuses, found {}".format(n)
    )

# Save retained IDs
with open(OUT_IDS, "w") as out:
    for seq_id in sorted(keep_ids):
        out.write(seq_id + "\n")

# Filter pooled FASTA
written = 0
found_ids = set()
current_id = None
current_lines = []

def flush_record(out):
    global written

    if current_id is not None and current_id in keep_ids:
        out.write("".join(current_lines))
        found_ids.add(current_id)
        written += 1

with open(FASTA) as inp, open(OUT_FASTA, "w") as out:
    for line in inp:
        if line.startswith(">"):
            flush_record(out)
            current_id = line[1:].strip().split()[0]
            current_lines = [line]
        else:
            current_lines.append(line)

    flush_record(out)

if found_ids != keep_ids:
    raise RuntimeError(
        "{} retained IDs missing from FASTA".format(
            len(keep_ids - found_ids)
        )
    )

# Write per-sample table
with open(OUT_SAMPLE, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "sample_id",
        "racon_polishing_only",
        "cutoff_0p5pct",
        "LC_cutoff_max_0p5pct_20reads",
        "LC_retained_consensuses"
    ])

    for sample in samples:
        writer.writerow([
            sample,
            racon_total[sample],
            int(plan_rows[sample]["cutoff_0p5pct"]),
            lc_cutoff[sample],
            lc_total[sample]
        ])

zero_samples = sum(1 for s in samples if lc_total[s] == 0)

print("===== FINAL LC FILTER COMPLETE =====")
print("LC definition: max(0.5%, 20 reads)")
print()
print("Starting Racon consensuses:", n)
print("LC retained consensuses:", len(keep_ids))
print("Samples represented:", 576 - zero_samples)
print("Zero-consensus samples:", zero_samples)
print()
print("20 reads was stricter:", cutoff_source["20_reads"])
print("0.5% was stricter:", cutoff_source["0.5_percent"])
print("Equal:", cutoff_source["equal"])
print()
print("Output FASTA:")
print(OUT_FASTA)
