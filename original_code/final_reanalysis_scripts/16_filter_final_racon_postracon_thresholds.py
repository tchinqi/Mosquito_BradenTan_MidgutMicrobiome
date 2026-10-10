#!/usr/bin/env python3

import csv
from pathlib import Path

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")

POOLED = RUN / "04_global_vsearch97/min2_complete576/pooled_racon_polished_576.fasta"
PROV   = RUN / "04_global_vsearch97/min2_complete576/pooled_racon_polished_576_provenance.tsv"
PLAN   = RUN / "02_abundance_ratio/threshold_plan.tsv"
OUT    = RUN / "04_global_vsearch97/postracon_thresholds_576"

OUT.mkdir(parents=True, exist_ok=True)

conditions = ["0p1pct", "0p2pct", "0p5pct", "abs5", "abs10", "abs20"]

# ------------------------------------------------------------
# 1. Read sample-specific thresholds
# ------------------------------------------------------------
cutoffs = {}

with open(PLAN, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        sample = row["sample_id"]

        cutoffs[sample] = {
            "0p1pct": int(row["cutoff_0p1pct"]),
            "0p2pct": int(row["cutoff_0p2pct"]),
            "0p5pct": int(row["cutoff_0p5pct"]),
            "abs5": 5,
            "abs10": int(row["cutoff_abs10"]),
            "abs20": int(row["cutoff_abs20"]),
        }

# ------------------------------------------------------------
# 2. Read exact final 22,528-consensus provenance
# ------------------------------------------------------------
records = []

with open(PROV, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        records.append({
            "sequence_id": row["sequence_id"],
            "sample_id": row["sample_id"],
            "supporting_reads": int(row["supporting_reads"]),
        })

if len(records) != 22528:
    raise RuntimeError(
        f"Expected 22528 final Racon consensuses, found {len(records)}"
    )

# ------------------------------------------------------------
# 3. Determine which consensus survives each threshold
# ------------------------------------------------------------
keep = {c: set() for c in conditions}
per_sample = {c: {} for c in conditions}

for rec in records:
    seqid = rec["sequence_id"]
    sample = rec["sample_id"]
    support = rec["supporting_reads"]

    if sample not in cutoffs:
        raise RuntimeError(f"Sample missing from threshold plan: {sample}")

    for condition in conditions:
        cutoff = cutoffs[sample][condition]

        if support >= cutoff:
            keep[condition].add(seqid)
            per_sample[condition][sample] = (
                per_sample[condition].get(sample, 0) + 1
            )

# ------------------------------------------------------------
# 4. Read pooled Racon FASTA
# ------------------------------------------------------------
seqs = {}
current = None
parts = []

with open(POOLED) as fh:
    for line in fh:
        line = line.rstrip()

        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(parts)

            current = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)

    if current is not None:
        seqs[current] = "".join(parts)

if len(seqs) != 22528:
    raise RuntimeError(
        f"Expected 22528 sequences in pooled FASTA, found {len(seqs)}"
    )

# Confirm every provenance sequence exists in FASTA
prov_ids = {r["sequence_id"] for r in records}
missing = prov_ids - set(seqs)

if missing:
    raise RuntimeError(
        f"{len(missing)} provenance sequence IDs are missing from pooled FASTA"
    )

# ------------------------------------------------------------
# 5. Write six post-Racon filtered FASTAs
# ------------------------------------------------------------
summary_rows = []

for condition in conditions:
    fasta_out = OUT / f"{condition}_postracon_filtered.fasta"
    ids_out   = OUT / f"{condition}_keep_ids.txt"

    ordered_ids = [
        rec["sequence_id"]
        for rec in records
        if rec["sequence_id"] in keep[condition]
    ]

    with open(ids_out, "w") as fh:
        for seqid in ordered_ids:
            fh.write(seqid + "\n")

    with open(fasta_out, "w") as fh:
        for seqid in ordered_ids:
            fh.write(f">{seqid}\n{seqs[seqid]}\n")

    represented = len(per_sample[condition])
    zero = 576 - represented

    summary_rows.append({
        "condition": condition,
        "retained_polished_consensuses": len(ordered_ids),
        "samples_with_ge1_consensus": represented,
        "samples_with_zero_consensuses": zero,
    })

# ------------------------------------------------------------
# 6. Write summary
# ------------------------------------------------------------
SUMMARY = OUT / "postracon_threshold_summary.tsv"

with open(SUMMARY, "w", newline="") as fh:
    fieldnames = [
        "condition",
        "retained_polished_consensuses",
        "samples_with_ge1_consensus",
        "samples_with_zero_consensuses",
    ]

    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

print("===== POST-RACON THRESHOLD FILTERING COMPLETE =====")
print(f"Starting Racon consensuses: {len(records)}")
print()
print(
    "condition\tretained_consensuses\t"
    "samples_with_consensus\tzero_samples"
)

for row in summary_rows:
    print(
        f"{row['condition']}\t"
        f"{row['retained_polished_consensuses']}\t"
        f"{row['samples_with_ge1_consensus']}\t"
        f"{row['samples_with_zero_consensuses']}"
    )
