#!/usr/bin/env python3
# =============================================================================
# Script 51 - 51_summarise_MC_GG97_wild_vs_lab_genera_final.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : R inputs
# Purpose : Mean genus relative abundance, wild vs lab (cross-check used by R
#           script 01).
# Input   : biological plotting table
# Output  : MC_GG97_wild_vs_lab_genus_mean_relative_abundance_FINAL.tsv
# Run on  : Atlas9 (Python 3)
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================

import csv
from pathlib import Path

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

INFILE = (
    RUN /
    "07_analysis/MC_GG2_97_plotting_final/"
    "MC_GG97_biological_570_samples_genus_relative_abundance.tsv"
)

OUTDIR = (
    RUN /
    "07_analysis/MC_GG2_97_plotting_final/"
    "wild_vs_lab_genus"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

SUMMARY = (
    OUTDIR /
    "MC_GG97_wild_vs_lab_genus_mean_relative_abundance_FINAL.tsv"
)

metadata_cols = [
    "Sample_ID",
    "Specimen_Code",
    "Source",
    "Plate",
    "Origin",
    "Location_Code",
    "Location_Name",
    "Location_Status",
    "Include_in_Location_Analysis",
    "Exclusion_or_Review_Reason",
    "Source_Sheet",
    "Source_Well",
    "Location_Crosscheck",
    "Mapping_Source",
]

with open(INFILE, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)

genus_cols = [
    x for x in fields
    if x not in metadata_cols
]

print("=" * 76)
print("FINAL INPUT")
print("=" * 76)

print("Samples:", len(rows))
print("Genus labels:", len(genus_cols))

if len(rows) != 570:
    raise RuntimeError(
        f"Expected 570 biological samples, found {len(rows)}"
    )

if len(genus_cols) != 68:
    raise RuntimeError(
        f"Expected 68 genus labels, found {len(genus_cols)}"
    )


# ------------------------------------------------------------
# Validate relative abundance sums
# ------------------------------------------------------------

bad = []

for r in rows:

    total = sum(
        float(r[g])
        for g in genus_cols
    )

    if abs(total - 1.0) > 1e-8:
        bad.append(
            (r["Sample_ID"], total)
        )

print(
    "Samples not summing to 1:",
    len(bad)
)

if bad:
    raise RuntimeError(
        "Relative abundance validation failed"
    )

print("All 570 samples sum to 1.0: PASS")


# ------------------------------------------------------------
# Wild / Lab
# ------------------------------------------------------------

wild = [
    r for r in rows
    if r["Origin"] == "Wild"
]

lab = [
    r for r in rows
    if r["Origin"] == "Lab-reared"
]

print()
print("=" * 76)
print("GROUP COUNTS")
print("=" * 76)

print("Wild:", len(wild))
print("Lab-reared:", len(lab))

if len(wild) != 380:
    raise RuntimeError("Expected 380 Wild")

if len(lab) != 190:
    raise RuntimeError("Expected 190 Lab-reared")

print("380 Wild + 190 Lab-reared = 570: PASS")


# ------------------------------------------------------------
# Means + prevalence
# ------------------------------------------------------------

results = []

for genus in genus_cols:

    overall_values = [
        float(r[genus])
        for r in rows
    ]

    wild_values = [
        float(r[genus])
        for r in wild
    ]

    lab_values = [
        float(r[genus])
        for r in lab
    ]

    results.append({
        "Genus": genus,

        "Overall_mean_relative_abundance":
            sum(overall_values) / len(overall_values),

        "Wild_mean_relative_abundance":
            sum(wild_values) / len(wild_values),

        "Lab_mean_relative_abundance":
            sum(lab_values) / len(lab_values),

        "Overall_prevalence":
            sum(x > 0 for x in overall_values)
            / len(overall_values),

        "Wild_prevalence":
            sum(x > 0 for x in wild_values)
            / len(wild_values),

        "Lab_prevalence":
            sum(x > 0 for x in lab_values)
            / len(lab_values),
    })


results.sort(
    key=lambda x:
    x["Overall_mean_relative_abundance"],
    reverse=True
)


with open(SUMMARY, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(results[0].keys()),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(results)


# ------------------------------------------------------------
# Print top 15
# ------------------------------------------------------------

print()
print("=" * 76)
print("TOP 15 GENERA — FINAL 570-SAMPLE DATASET")
print("=" * 76)

print(
    f"{'Genus':35s}"
    f"{'Overall %':>12s}"
    f"{'Wild %':>12s}"
    f"{'Lab %':>12s}"
    f"{'Prev %':>12s}"
)

for r in results[:15]:

    print(
        f"{r['Genus'][:35]:35s}"
        f"{100*r['Overall_mean_relative_abundance']:12.2f}"
        f"{100*r['Wild_mean_relative_abundance']:12.2f}"
        f"{100*r['Lab_mean_relative_abundance']:12.2f}"
        f"{100*r['Overall_prevalence']:12.1f}"
    )


print()
print("=" * 76)
print("CHECKSUMS")
print("=" * 76)

overall_sum = sum(
    r["Overall_mean_relative_abundance"]
    for r in results
)

wild_sum = sum(
    r["Wild_mean_relative_abundance"]
    for r in results
)

lab_sum = sum(
    r["Lab_mean_relative_abundance"]
    for r in results
)

print("Overall:", overall_sum)
print("Wild:", wild_sum)
print("Lab:", lab_sum)

if abs(overall_sum - 1) > 1e-8:
    raise RuntimeError("Overall means do not sum to 1")

if abs(wild_sum - 1) > 1e-8:
    raise RuntimeError("Wild means do not sum to 1")

if abs(lab_sum - 1) > 1e-8:
    raise RuntimeError("Lab means do not sum to 1")

print("All composition sums = 1: PASS")

print()
print("FINAL SUMMARY:")
print(SUMMARY)

