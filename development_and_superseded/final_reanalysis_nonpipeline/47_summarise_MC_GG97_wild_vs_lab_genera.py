#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

INFILE = (
    RUN /
    "07_analysis/MC_GG2_97_plotting/"
    "MC_GG97_biological_572_samples_genus_relative_abundance.tsv"
)

OUTDIR = (
    RUN /
    "07_analysis/MC_GG2_97_plotting/"
    "wild_vs_lab_genus"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

SUMMARY = (
    OUTDIR /
    "MC_GG97_wild_vs_lab_genus_mean_relative_abundance.tsv"
)


# ============================================================
# READ TABLE
# ============================================================

with open(INFILE, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames
    rows = list(reader)

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

genus_cols = [
    x for x in fields
    if x not in metadata_cols
]


print("=" * 76)
print("INPUT")
print("=" * 76)

print("Samples:", len(rows))
print("Genus labels:", len(genus_cols))

if len(rows) != 572:
    raise RuntimeError(
        f"Expected 572 biological samples, found {len(rows)}"
    )

if len(genus_cols) != 68:
    raise RuntimeError(
        f"Expected 68 genus labels, found {len(genus_cols)}"
    )


# ============================================================
# VALIDATE RELATIVE ABUNDANCE
# ============================================================

bad_samples = []

for r in rows:

    total = sum(
        float(r[g])
        for g in genus_cols
    )

    if abs(total - 1.0) > 1e-8:
        bad_samples.append(
            (r["Sample_ID"], total)
        )

print()
print("Samples whose genus relative abundances do not sum to 1:",
      len(bad_samples))

if bad_samples:
    print(bad_samples[:10])
    raise RuntimeError(
        "Relative abundance validation failed"
    )

print("All 572 samples sum to 1.0: PASS")


# ============================================================
# GROUP BY ORIGIN
# ============================================================

groups = {
    "Wild": [
        r for r in rows
        if r["Origin"] == "Wild"
    ],

    "Lab-reared": [
        r for r in rows
        if r["Origin"] == "Lab-reared"
    ],
}

print()
print("=" * 76)
print("GROUP COUNTS")
print("=" * 76)

for group, group_rows in groups.items():
    print(group, len(group_rows))

if len(groups["Wild"]) != 380:
    raise RuntimeError("Wild count is not 380")

if len(groups["Lab-reared"]) != 192:
    raise RuntimeError("Lab-reared count is not 192")

print("Wild/Lab sample counts: PASS")


# ============================================================
# CALCULATE MEAN RELATIVE ABUNDANCE + PREVALENCE
# ============================================================

results = []

for genus in genus_cols:

    all_values = [
        float(r[genus])
        for r in rows
    ]

    wild_values = [
        float(r[genus])
        for r in groups["Wild"]
    ]

    lab_values = [
        float(r[genus])
        for r in groups["Lab-reared"]
    ]

    overall_mean = (
        sum(all_values) /
        len(all_values)
    )

    wild_mean = (
        sum(wild_values) /
        len(wild_values)
    )

    lab_mean = (
        sum(lab_values) /
        len(lab_values)
    )

    overall_prev = sum(
        x > 0
        for x in all_values
    ) / len(all_values)

    wild_prev = sum(
        x > 0
        for x in wild_values
    ) / len(wild_values)

    lab_prev = sum(
        x > 0
        for x in lab_values
    ) / len(lab_values)

    results.append({
        "Genus": genus,
        "Overall_mean_relative_abundance": overall_mean,
        "Wild_mean_relative_abundance": wild_mean,
        "Lab_mean_relative_abundance": lab_mean,
        "Overall_prevalence": overall_prev,
        "Wild_prevalence": wild_prev,
        "Lab_prevalence": lab_prev,
    })


results.sort(
    key=lambda x:
    x["Overall_mean_relative_abundance"],
    reverse=True
)


# ============================================================
# WRITE SUMMARY
# ============================================================

with open(SUMMARY, "w", newline="") as fh:

    fields_out = [
        "Genus",
        "Overall_mean_relative_abundance",
        "Wild_mean_relative_abundance",
        "Lab_mean_relative_abundance",
        "Overall_prevalence",
        "Wild_prevalence",
        "Lab_prevalence",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields_out,
        delimiter="\t"
    )

    writer.writeheader()

    for r in results:
        writer.writerow(r)


# ============================================================
# PRINT TOP 15
# ============================================================

print()
print("=" * 76)
print("TOP 15 GENERA BY OVERALL MEAN RELATIVE ABUNDANCE")
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

print("Overall mean composition sum:", overall_sum)
print("Wild mean composition sum:", wild_sum)
print("Lab mean composition sum:", lab_sum)

if abs(overall_sum - 1) > 1e-8:
    raise RuntimeError("Overall genus means do not sum to 1")

if abs(wild_sum - 1) > 1e-8:
    raise RuntimeError("Wild genus means do not sum to 1")

if abs(lab_sum - 1) > 1e-8:
    raise RuntimeError("Lab genus means do not sum to 1")

print("All composition sums = 1: PASS")


print()
print("=" * 76)
print("OUTPUT")
print("=" * 76)

print(SUMMARY)

print()
print("WILD VS LAB GENUS SUMMARY COMPLETE")

