#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

READY = (
    RUN /
    "07_analysis/MC_GG2_97_analysis_ready/"
    "MC_GG2_97_classified_576_samples_x_328_OTUs_with_metadata.tsv"
)

TAXONOMY = (
    RUN /
    "07_analysis/MC_GG2_97_classified_OTU_table/"
    "MC_GG2_97_classified_328_OTU_taxonomy.tsv"
)

OUTDIR = (
    RUN /
    "07_analysis/MC_GG2_97_plotting_final"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. READ TAXONOMY AND EXTRACT GENUS
# ============================================================

otu_to_taxon = {}
otu_to_genus = {}

with open(TAXONOMY, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:

        otu = row["Feature_ID"].strip()
        taxon = row["GG2_97_Taxon"].strip()

        otu_to_taxon[otu] = taxon

        genus = None

        for part in taxon.split(";"):

            part = part.strip()

            if part.startswith("g__"):

                value = part[3:].strip()

                if value:
                    genus = value

        if genus is None:
            genus = "Unclassified_at_genus"

        otu_to_genus[otu] = genus


print("=" * 76)
print("TAXONOMY")
print("=" * 76)

print("Classified OTUs:", len(otu_to_genus))
print("Unique genus labels:", len(set(otu_to_genus.values())))

if len(otu_to_genus) != 328:
    raise RuntimeError(
        f"Expected 328 classified OTUs, found {len(otu_to_genus)}"
    )

print("328-OTU taxonomy validation: PASS")


# ============================================================
# 2. READ ANALYSIS-READY TABLE
# ============================================================

with open(READY, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    fields = reader.fieldnames

    otu_cols = [
        x for x in fields
        if x.startswith("OTU_")
    ]

    metadata_cols = [
        x for x in fields
        if not x.startswith("OTU_")
    ]

    rows = list(reader)


print()
print("=" * 76)
print("INPUT TABLE")
print("=" * 76)

print("Samples:", len(rows))
print("Metadata columns:", len(metadata_cols))
print("OTU columns:", len(otu_cols))

if len(rows) != 576:
    raise RuntimeError(
        f"Expected 576 samples, found {len(rows)}"
    )

if len(otu_cols) != 328:
    raise RuntimeError(
        f"Expected 328 OTUs, found {len(otu_cols)}"
    )

if set(otu_cols) != set(otu_to_genus):
    raise RuntimeError(
        "OTU columns do not exactly match taxonomy OTUs"
    )

print("576 x 328 validation: PASS")


# ============================================================
# 3. DEFINE FINAL ANALYSIS SUBSETS
# ============================================================

biological = [
    r for r in rows
    if r["Origin"] in {"Wild", "Lab-reared"}
]

wild = [
    r for r in rows
    if r["Origin"] == "Wild"
]

lab = [
    r for r in rows
    if r["Origin"] == "Lab-reared"
]

controls = [
    r for r in rows
    if r["Origin"] == "Control"
]

location = [
    r for r in rows
    if (
        r["Origin"] == "Wild"
        and r["Include_in_Location_Analysis"] == "YES"
    )
]

wild_review = [
    r for r in rows
    if (
        r["Origin"] == "Wild"
        and r["Include_in_Location_Analysis"] != "YES"
    )
]


print()
print("=" * 76)
print("FINAL PLOTTING SUBSETS")
print("=" * 76)

print("Biological Wild + Lab-reared:", len(biological))
print("Wild overall:", len(wild))
print("Lab-reared:", len(lab))
print("Negative controls:", len(controls))
print("Wild mapped for location analysis:", len(location))
print("Wild excluded/review for location only:", len(wild_review))

if len(biological) != 570:
    raise RuntimeError(
        f"Expected 570 biological samples, found {len(biological)}"
    )

if len(wild) != 380:
    raise RuntimeError(
        f"Expected 380 Wild samples, found {len(wild)}"
    )

if len(lab) != 190:
    raise RuntimeError(
        f"Expected 190 Lab-reared samples, found {len(lab)}"
    )

if len(controls) != 6:
    raise RuntimeError(
        f"Expected 6 controls, found {len(controls)}"
    )

if len(location) != 345:
    raise RuntimeError(
        f"Expected 345 mapped Wild samples, found {len(location)}"
    )

if len(wild_review) != 35:
    raise RuntimeError(
        f"Expected 35 Wild review samples, found {len(wild_review)}"
    )

print("All final subset counts: PASS")


# ============================================================
# 4. WRITE OTU TABLES
# ============================================================

def write_otu_table(filename, subset):

    path = OUTDIR / filename

    with open(path, "w", newline="") as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=metadata_cols + otu_cols,
            delimiter="\t"
        )

        writer.writeheader()

        for r in subset:
            writer.writerow(r)

    return path


files = []

files.append(
    write_otu_table(
        "MC_GG97_biological_570_samples_x_328_OTUs.tsv",
        biological
    )
)

files.append(
    write_otu_table(
        "MC_GG97_wild_380_samples_x_328_OTUs.tsv",
        wild
    )
)

files.append(
    write_otu_table(
        "MC_GG97_lab_190_samples_x_328_OTUs.tsv",
        lab
    )
)

files.append(
    write_otu_table(
        "MC_GG97_location_345_samples_x_328_OTUs.tsv",
        location
    )
)

files.append(
    write_otu_table(
        "MC_GG97_controls_6_samples_x_328_OTUs.tsv",
        controls
    )
)


# ============================================================
# 5. GENUS TABLES FOR 570 BIOLOGICAL SAMPLES
# ============================================================

genera = sorted(set(otu_to_genus.values()))

GENUS_COUNTS = (
    OUTDIR /
    "MC_GG97_biological_570_samples_genus_counts.tsv"
)

GENUS_REL = (
    OUTDIR /
    "MC_GG97_biological_570_samples_genus_relative_abundance.tsv"
)


with open(GENUS_COUNTS, "w", newline="") as counts_fh, \
     open(GENUS_REL, "w", newline="") as rel_fh:

    counts_writer = csv.writer(
        counts_fh,
        delimiter="\t"
    )

    rel_writer = csv.writer(
        rel_fh,
        delimiter="\t"
    )

    header = metadata_cols + genera

    counts_writer.writerow(header)
    rel_writer.writerow(header)

    for r in biological:

        genus_counts = defaultdict(int)

        for otu in otu_cols:

            genus = otu_to_genus[otu]
            genus_counts[genus] += int(r[otu])

        total = sum(genus_counts.values())

        if total <= 0:
            raise RuntimeError(
                f"Zero abundance biological sample: "
                f"{r['Sample_ID']}"
            )

        meta_values = [
            r[x]
            for x in metadata_cols
        ]

        count_values = [
            genus_counts[g]
            for g in genera
        ]

        rel_values = [
            genus_counts[g] / total
            for g in genera
        ]

        counts_writer.writerow(
            meta_values + count_values
        )

        rel_writer.writerow(
            meta_values + rel_values
        )


# ============================================================
# 6. OTU -> GENUS MAP
# ============================================================

GENUS_MAP = (
    OUTDIR /
    "MC_GG97_328_OTU_to_genus.tsv"
)

with open(GENUS_MAP, "w", newline="") as fh:

    w = csv.writer(
        fh,
        delimiter="\t"
    )

    w.writerow([
        "Feature_ID",
        "GG2_97_Taxon",
        "Genus_for_plotting",
    ])

    for otu in otu_cols:

        w.writerow([
            otu,
            otu_to_taxon[otu],
            otu_to_genus[otu],
        ])


# ============================================================
# 7. LOCATION CHECK
# ============================================================

location_counts = defaultdict(int)

for r in location:
    location_counts[r["Location_Name"]] += 1


print()
print("=" * 76)
print("LOCATION COUNTS")
print("=" * 76)

for location_name, n in sorted(location_counts.items()):
    print(f"{location_name}: {n}")

if sum(location_counts.values()) != 345:
    raise RuntimeError(
        "Location totals do not sum to 345"
    )


# ============================================================
# 8. CONTROL CHECK
# ============================================================

print()
print("=" * 76)
print("SIX NEGATIVE CONTROLS")
print("=" * 76)

for r in controls:
    print(
        r["Sample_ID"],
        "|",
        r["Source"],
        r["Plate"],
        r["Source_Well"]
    )


# ============================================================
# 9. FINAL
# ============================================================

print()
print("=" * 76)
print("FINAL PLOT-READY FILES CREATED")
print("=" * 76)

for f in files:
    print(f)

print(GENUS_COUNTS)
print(GENUS_REL)
print(GENUS_MAP)

print()
print("Final biological dataset:")
print("  570 samples = 380 Wild + 190 Lab-reared")

print()
print("Location dataset:")
print("  345 mapped Wild samples")

print()
print("Negative controls:")
print("  6, retained separately and excluded from biological plots")

print()
print("=" * 76)
print("FINAL MC GG2 >=97% PLOTTING PREPARATION COMPLETE")
print("=" * 76)

