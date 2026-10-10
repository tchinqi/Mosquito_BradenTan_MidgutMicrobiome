#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict
from statistics import median

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

FULL_MATRIX = (
    RUN /
    "07_analysis/MC_GG2_97_classified_OTU_table/"
    "MC_867_OTUs_x_576_samples_supporting_reads.tsv"
)

CLASSIFIED_MATRIX = (
    RUN /
    "07_analysis/MC_GG2_97_classified_OTU_table/"
    "MC_GG2_97_classified_576_samples_x_328_OTUs.tsv"
)

TAXONOMY = (
    RUN /
    "06_taxonomy/greengenes2_MC_LC_final/"
    "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv"
)

PER_SAMPLE = (
    RUN /
    "04_global_vsearch97/postracon_thresholds_576/"
    "LC_max0p5pct_20reads/"
    "LC_max0p5pct_20reads_per_sample.tsv"
)

METADATA = (
    RUN /
    "metadata/sample_metadata_current.tsv"
)

OUTDIR = (
    RUN /
    "07_analysis/MC_GG2_97_control_audit"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = OUTDIR / "MC_control_candidate_read_audit.tsv"


KNOWN_CONTROLS = [
    "sample_Neg_AG_P01_16S",
    "sample_Negative_AG_P02_16S",
    "sample_Negative_ZX_P01_ZX_P01_16S",
    "sample_Negative_ZX_P02_16S",
]

BT_CANDIDATES = [
    "sample_BT_CLX_096_16S_BT_P01_16S",
    "sample_BT_CLX_120_16S_BT_P02_16S",
]

AUDIT_SAMPLES = KNOWN_CONTROLS + BT_CANDIDATES


# ============================================================
# METADATA
# ============================================================

metadata = {}

with open(METADATA, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        metadata[row["Sample_ID"]] = row


# ============================================================
# MC PER-SAMPLE CONSENSUS INFORMATION
# ============================================================

per_sample = {}

with open(PER_SAMPLE, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        per_sample[row["sample_id"]] = row


# ============================================================
# GG97 TAXONOMY / GENUS
# ============================================================

taxonomy = {}
otu_genus = {}

with open(TAXONOMY, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        otu = row["Feature ID"].strip()
        taxon = row["Taxon"].strip()
        taxonomy[otu] = taxon

        genus = None

        if taxon.lower() != "unassigned":
            for part in taxon.split(";"):
                part = part.strip()
                if part.startswith("g__") and len(part) > 3:
                    genus = part[3:]

        if genus is None:
            genus = "Unclassified_at_genus"

        otu_genus[otu] = genus


# ============================================================
# FULL 867-OTU SUPPORTING-READ MATRIX
# ============================================================

full_counts = {
    sid: {}
    for sid in metadata
}

with open(FULL_MATRIX, newline="") as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)

    sample_ids = header[1:]

    for row in reader:
        otu = row[0]

        for sid, value in zip(sample_ids, row[1:]):
            full_counts[sid][otu] = int(value)


# ============================================================
# CLASSIFIED 328-OTU MATRIX
# ============================================================

classified_counts = {}

with open(CLASSIFIED_MATRIX, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    classified_otus = reader.fieldnames[1:]

    for row in reader:
        sid = row["sample_id"]
        classified_counts[sid] = {
            otu: int(row[otu])
            for otu in classified_otus
        }


# ============================================================
# SUMMARISE EACH SAMPLE
# ============================================================

def summarise_sample(sid):

    all_counts = full_counts[sid]
    class_counts = classified_counts[sid]

    total_reads = sum(all_counts.values())
    classified_reads = sum(class_counts.values())

    observed_otus = sum(
        x > 0 for x in all_counts.values()
    )

    classified_otus_observed = sum(
        x > 0 for x in class_counts.values()
    )

    unassigned_reads = (
        total_reads - classified_reads
    )

    classified_pct = (
        100 * classified_reads / total_reads
        if total_reads > 0
        else 0
    )

    genus_counts = defaultdict(int)

    for otu, count in class_counts.items():
        if count > 0:
            genus_counts[otu_genus[otu]] += count

    top = sorted(
        genus_counts.items(),
        key=lambda x: x[1],
        reverse=True
    )[:5]

    top_text = "; ".join(
        f"{g}={n}"
        for g, n in top
    )

    ps = per_sample[sid]

    return {
        "Sample_ID": sid,
        "Status_for_audit": (
            "Known negative control"
            if sid in KNOWN_CONTROLS
            else "BT candidate control"
        ),
        "Current_Origin": metadata[sid]["Origin"],
        "Source": metadata[sid]["Source"],
        "Plate": metadata[sid]["Plate"],
        "MC_retained_consensuses":
            int(ps["LC_retained_consensuses"]),
        "Total_supporting_reads_867_OTUs":
            total_reads,
        "Observed_OTUs_867":
            observed_otus,
        "GG97_classified_supporting_reads":
            classified_reads,
        "GG97_unassigned_supporting_reads":
            unassigned_reads,
        "GG97_classified_percent":
            classified_pct,
        "Observed_GG97_classified_OTUs":
            classified_otus_observed,
        "Top_5_classified_genera":
            top_text,
    }


summaries = [
    summarise_sample(sid)
    for sid in AUDIT_SAMPLES
]


# ============================================================
# REFERENCE DISTRIBUTIONS
# ============================================================

known_summary = [
    x for x in summaries
    if x["Status_for_audit"] == "Known negative control"
]

bt_biological_ids = [
    sid
    for sid, row in metadata.items()
    if (
        row["Source"] == "BT"
        and row["Origin"] == "Lab-reared"
        and sid not in BT_CANDIDATES
    )
]

bt_bio_summary = [
    summarise_sample(sid)
    for sid in bt_biological_ids
]


metrics = [
    "MC_retained_consensuses",
    "Total_supporting_reads_867_OTUs",
    "Observed_OTUs_867",
    "GG97_classified_supporting_reads",
    "GG97_classified_percent",
    "Observed_GG97_classified_OTUs",
]


print("=" * 100)
print("CONTROL / CANDIDATE SAMPLE AUDIT")
print("=" * 100)

for x in summaries:

    print()
    print(x["Sample_ID"])
    print("  Audit status:", x["Status_for_audit"])
    print("  Current metadata origin:", x["Current_Origin"])
    print("  Source / Plate:",
          x["Source"], x["Plate"])
    print("  MC retained consensuses:",
          x["MC_retained_consensuses"])
    print("  Total supporting reads:",
          x["Total_supporting_reads_867_OTUs"])
    print("  Observed MC OTUs:",
          x["Observed_OTUs_867"])
    print("  GG97 classified reads:",
          x["GG97_classified_supporting_reads"])
    print("  GG97 unassigned reads:",
          x["GG97_unassigned_supporting_reads"])
    print(
        "  GG97 classified %:",
        f"{x['GG97_classified_percent']:.2f}%"
    )
    print(
        "  Observed classified OTUs:",
        x["Observed_GG97_classified_OTUs"]
    )
    print(
        "  Top classified genera:",
        x["Top_5_classified_genera"]
    )


print()
print("=" * 100)
print("REFERENCE MEDIANS")
print("=" * 100)

print()
print(
    "Known negative controls (n =",
    len(known_summary), ")"
)

for m in metrics:
    print(
        f"  {m}:",
        median(x[m] for x in known_summary)
    )

print()
print(
    "Ordinary BT lab-reared samples excluding candidates (n =",
    len(bt_bio_summary), ")"
)

for m in metrics:
    print(
        f"  {m}:",
        median(x[m] for x in bt_bio_summary)
    )


# ============================================================
# WRITE AUDIT TABLE
# ============================================================

fields = list(summaries[0].keys())

with open(OUTFILE, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summaries)


print()
print("=" * 100)
print("OUTPUT")
print("=" * 100)

print(OUTFILE)

