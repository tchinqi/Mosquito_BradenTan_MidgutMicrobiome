#!/usr/bin/env python3
# =============================================================================
# Script 44 - 44_make_MC_GG97_classified_OTU_sample_table.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : OTU x sample table
# Purpose : Builds the supporting-read table of the 328 MC OTUs classified by
#           Greengenes2 at >= 97% across the 576 wells (reads per consensus
#           summed per OTU).
# Input   : MC per-sample table; otu97_clusters.uc; centroids; provenance; GG2
#           97% taxonomy
# Output  : 07_analysis/MC_GG2_97_classified_OTU_table/
# Run on  : Atlas9 (Python 3)
# Note    : Historical file/folder names swap LC and MC: '0p1pct' /
#           'MC_0p1pct' = final LC (0.1%); 'LC_max0p5pct_20reads' = final MC
#           (max(0.5%, 20 reads)).
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================

import csv
from pathlib import Path
from collections import defaultdict, Counter

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

MC = (
    RUN /
    "04_global_vsearch97/postracon_thresholds_576/"
    "LC_max0p5pct_20reads"
)

PROVENANCE = (
    RUN /
    "04_global_vsearch97/min2_complete576/"
    "pooled_racon_polished_576_provenance.tsv"
)

PER_SAMPLE = (
    MC /
    "LC_max0p5pct_20reads_per_sample.tsv"
)

UC = (
    MC /
    "vsearch97/otu97_clusters.uc"
)

CENTROIDS = (
    MC /
    "vsearch97/otu97_centroids_UPPERCASE.fasta"
)

# This old directory name corresponds to the CURRENT MC branch:
# max(0.5%, 20 reads), 867 OTUs.
# This taxonomy file is the GG2 >=97% identity classification.
TAXONOMY = (
    RUN /
    "06_taxonomy/greengenes2_MC_LC_final/"
    "LC_max0p5pct_20reads/"
    "exported_taxonomy/taxonomy.tsv"
)

OUTDIR = (
    RUN /
    "07_analysis/MC_GG2_97_classified_OTU_table"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. READ 576-SAMPLE ORDER + EXPECTED RETAINED CONSENSUS COUNTS
# ============================================================

samples = []
expected_consensus_count = {}

with open(PER_SAMPLE, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:

        sample = row["sample_id"].strip()

        samples.append(sample)

        expected_consensus_count[sample] = int(
            row["LC_retained_consensuses"]
        )


print("=" * 76)
print("SAMPLE TABLE")
print("=" * 76)
print("Samples:", len(samples))
print("Unique samples:", len(set(samples)))

if len(samples) != 576:
    raise RuntimeError(
        f"Expected 576 samples, found {len(samples)}"
    )

if len(set(samples)) != 576:
    raise RuntimeError("Duplicate sample IDs detected")

print("576-sample validation: PASS")


# ============================================================
# 2. READ RACON PROVENANCE
# ============================================================

provenance = {}

with open(PROVENANCE, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:

        sequence_id = row["sequence_id"].strip()

        if sequence_id in provenance:
            raise RuntimeError(
                f"Duplicate provenance ID: {sequence_id}"
            )

        provenance[sequence_id] = {
            "sample": row["sample_id"].strip(),
            "supporting_reads": int(row["supporting_reads"]),
        }


print()
print("=" * 76)
print("PROVENANCE")
print("=" * 76)
print("Racon sequences in provenance:", len(provenance))


# ============================================================
# 3. READ CENTROID IDs
# ============================================================

centroid_ids = []

with open(CENTROIDS) as fh:

    for line in fh:

        if line.startswith(">"):
            centroid_ids.append(
                line[1:].strip().split()[0]
            )


print()
print("=" * 76)
print("MC OTU CENTROIDS")
print("=" * 76)
print("Centroid OTUs:", len(centroid_ids))

if len(centroid_ids) != 867:
    raise RuntimeError(
        f"Expected 867 MC OTUs, found {len(centroid_ids)}"
    )

expected_otu_ids = [
    f"OTU_{i}"
    for i in range(1, 868)
]

if centroid_ids != expected_otu_ids:
    raise RuntimeError(
        "Centroid IDs are not exactly OTU_1 ... OTU_867 "
        "in expected order"
    )

print("OTU_1 ... OTU_867 validation: PASS")


# ============================================================
# 4. READ VSEARCH .UC:
#    RACON SEQUENCE -> 97% CLUSTER -> OTU ID
#
# Cluster numbers are zero-based:
# cluster 0 -> OTU_1
# cluster 1 -> OTU_2
# ...
# ============================================================

racon_to_otu = {}
seed_clusters = set()

with open(UC) as fh:

    for line in fh:

        line = line.rstrip("\n")

        if not line:
            continue

        parts = line.split("\t")

        record_type = parts[0]

        if record_type not in {"S", "H"}:
            continue

        cluster_number = int(parts[1])
        sequence_id = parts[8].strip()

        otu_id = f"OTU_{cluster_number + 1}"

        if sequence_id in racon_to_otu:
            raise RuntimeError(
                f"Racon sequence appears multiple times in UC: "
                f"{sequence_id}"
            )

        racon_to_otu[sequence_id] = otu_id

        if record_type == "S":
            seed_clusters.add(cluster_number)


print()
print("=" * 76)
print("VSEARCH97 CLUSTERS")
print("=" * 76)
print("Racon sequences in UC:", len(racon_to_otu))
print("Seed clusters:", len(seed_clusters))
print(
    "Cluster range:",
    min(seed_clusters),
    "to",
    max(seed_clusters),
)

if len(seed_clusters) != 867:
    raise RuntimeError(
        f"Expected 867 seed clusters, found "
        f"{len(seed_clusters)}"
    )

if seed_clusters != set(range(867)):
    raise RuntimeError(
        "VSEARCH cluster IDs are not exactly 0...866"
    )

if set(racon_to_otu.values()) != set(centroid_ids):
    raise RuntimeError(
        "UC-derived OTU IDs do not match centroid OTU IDs"
    )

print("UC -> OTU numbering validation: PASS")


# ============================================================
# 5. VERIFY EVERY RETAINED RACON HAS PROVENANCE
# ============================================================

missing_provenance = [
    x
    for x in racon_to_otu
    if x not in provenance
]

print()
print("Racon sequences missing provenance:",
      len(missing_provenance))

if missing_provenance:
    raise RuntimeError(
        "Some MC-retained Racon sequences lack provenance"
    )


# ============================================================
# 6. VERIFY PER-SAMPLE RETAINED CONSENSUS COUNTS
# ============================================================

observed_consensus_count = Counter()

for sequence_id in racon_to_otu:

    sample = provenance[sequence_id]["sample"]

    if sample not in expected_consensus_count:
        raise RuntimeError(
            f"Unexpected sample in provenance: {sample}"
        )

    observed_consensus_count[sample] += 1


mismatches = []

for sample in samples:

    expected = expected_consensus_count[sample]
    observed = observed_consensus_count[sample]

    if expected != observed:
        mismatches.append(
            (sample, expected, observed)
        )


print()
print("=" * 76)
print("PER-SAMPLE CONSENSUS VALIDATION")
print("=" * 76)

print(
    "Expected retained consensuses:",
    sum(expected_consensus_count.values())
)

print(
    "Observed UC consensuses:",
    sum(observed_consensus_count.values())
)

print("Samples with mismatched counts:", len(mismatches))

if mismatches:

    for x in mismatches[:20]:
        print("MISMATCH:", x)

    raise RuntimeError(
        "Per-sample retained consensus counts do not match UC"
    )

print("Per-sample retained consensus counts: PASS")


# ============================================================
# 7. BUILD FULL 867 OTU x 576 SAMPLE SUPPORTING-READ MATRIX
# ============================================================

matrix = {
    otu: defaultdict(int)
    for otu in centroid_ids
}

total_support = 0

for sequence_id, otu in racon_to_otu.items():

    info = provenance[sequence_id]

    sample = info["sample"]
    reads = info["supporting_reads"]

    matrix[otu][sample] += reads
    total_support += reads


print()
print("=" * 76)
print("FULL MC ABUNDANCE MATRIX")
print("=" * 76)
print("OTUs:", len(matrix))
print("Samples:", len(samples))
print("Total supporting-read abundance:", total_support)

zero_otu_rows = [
    otu
    for otu in centroid_ids
    if sum(matrix[otu].values()) == 0
]

if zero_otu_rows:
    raise RuntimeError(
        f"{len(zero_otu_rows)} OTUs have zero abundance"
    )

print("All 867 OTUs have >0 abundance: PASS")


# ============================================================
# 8. READ GG2 >=97% TAXONOMY
# ============================================================

taxonomy = {}

with open(TAXONOMY, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:

        feature_id = row["Feature ID"].strip()
        taxon = (row["Taxon"] or "").strip()
        consensus = (row["Consensus"] or "").strip()

        taxonomy[feature_id] = {
            "taxon": taxon,
            "consensus": consensus,
        }


print()
print("=" * 76)
print("GREENGENES2 >=97% TAXONOMY")
print("=" * 76)
print("Taxonomy OTUs:", len(taxonomy))

if set(taxonomy) != set(centroid_ids):
    missing_tax = set(centroid_ids) - set(taxonomy)
    extra_tax = set(taxonomy) - set(centroid_ids)

    print("Missing taxonomy IDs:", len(missing_tax))
    print("Extra taxonomy IDs:", len(extra_tax))

    raise RuntimeError(
        "GG2 taxonomy IDs do not exactly match MC OTUs"
    )


unassigned = [
    otu
    for otu in centroid_ids
    if taxonomy[otu]["taxon"].lower() == "unassigned"
]

classified = [
    otu
    for otu in centroid_ids
    if taxonomy[otu]["taxon"].lower() != "unassigned"
]


print("Unassigned OTUs:", len(unassigned))
print("Classified OTUs:", len(classified))

if len(unassigned) != 539:
    raise RuntimeError(
        f"Expected 539 Unassigned OTUs, "
        f"found {len(unassigned)}"
    )

if len(classified) != 328:
    raise RuntimeError(
        f"Expected 328 classified OTUs, "
        f"found {len(classified)}"
    )

print("539 Unassigned + 328 classified = 867: PASS")


# ============================================================
# 9. WRITE FULL MC 867 OTU x 576 SAMPLE TABLE
# ============================================================

full_file = (
    OUTDIR /
    "MC_867_OTUs_x_576_samples_supporting_reads.tsv"
)

with open(full_file, "w", newline="") as fh:

    w = csv.writer(fh, delimiter="\t")

    w.writerow(
        ["Feature_ID"] + samples
    )

    for otu in centroid_ids:

        w.writerow(
            [otu] +
            [
                matrix[otu].get(sample, 0)
                for sample in samples
            ]
        )


# ============================================================
# 10. WRITE PRIYA'S FINAL:
#     328 CLASSIFIED OTUs x 576 SAMPLES
# ============================================================

classified_file = (
    OUTDIR /
    "MC_GG2_97_classified_328_OTUs_x_576_samples.tsv"
)

with open(classified_file, "w", newline="") as fh:

    w = csv.writer(fh, delimiter="\t")

    w.writerow(
        ["Feature_ID"] + samples
    )

    for otu in classified:

        w.writerow(
            [otu] +
            [
                matrix[otu].get(sample, 0)
                for sample in samples
            ]
        )


# ============================================================
# 11. WRITE TRANSPOSE:
#     576 SAMPLES x 328 CLASSIFIED OTUs
# ============================================================

transpose_file = (
    OUTDIR /
    "MC_GG2_97_classified_576_samples_x_328_OTUs.tsv"
)

with open(transpose_file, "w", newline="") as fh:

    w = csv.writer(fh, delimiter="\t")

    w.writerow(
        ["sample_id"] + classified
    )

    for sample in samples:

        w.writerow(
            [sample] +
            [
                matrix[otu].get(sample, 0)
                for otu in classified
            ]
        )


# ============================================================
# 12. WRITE TAXONOMY MAP FOR THE 328 CLASSIFIED OTUs
# ============================================================

taxonomy_file = (
    OUTDIR /
    "MC_GG2_97_classified_328_OTU_taxonomy.tsv"
)

with open(taxonomy_file, "w", newline="") as fh:

    w = csv.writer(fh, delimiter="\t")

    w.writerow([
        "Feature_ID",
        "GG2_97_Taxon",
        "Consensus",
        "Total_supporting_reads",
        "Prevalence_samples",
    ])

    for otu in classified:

        total = sum(
            matrix[otu].get(sample, 0)
            for sample in samples
        )

        prevalence = sum(
            matrix[otu].get(sample, 0) > 0
            for sample in samples
        )

        w.writerow([
            otu,
            taxonomy[otu]["taxon"],
            taxonomy[otu]["consensus"],
            total,
            prevalence,
        ])


# ============================================================
# 13. FINAL VALIDATION
# ============================================================

classified_support = sum(
    matrix[otu].get(sample, 0)
    for otu in classified
    for sample in samples
)

samples_with_classified_reads = sum(
    any(
        matrix[otu].get(sample, 0) > 0
        for otu in classified
    )
    for sample in samples
)


print()
print("=" * 76)
print("FINAL OUTPUT")
print("=" * 76)

print("Full matrix:")
print("  867 OTUs x 576 samples")

print()
print("GG2 >=97%-classified matrix:")
print("  328 OTUs x 576 samples")

print()
print(
    "Total supporting-read abundance in "
    "classified matrix:",
    classified_support
)

print(
    "Samples containing >=1 GG97-classified OTU:",
    samples_with_classified_reads,
    "/ 576"
)

print()
print("OUTPUT FILES")

print()
print(full_file)

print()
print(classified_file)

print()
print(transpose_file)

print()
print(taxonomy_file)

print()
print("=" * 76)
print("MC GG2 >=97% CLASSIFIED OTU x SAMPLE TABLE COMPLETE")
print("=" * 76)

