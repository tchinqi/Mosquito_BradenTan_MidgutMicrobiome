#!/usr/bin/env python3
# =============================================================================
# Script 45 - 45_join_MC_GG97_table_to_metadata.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Metadata
# Purpose : Joins the OTU table to sample metadata (origin, plate, well,
#           collection site).
# Input   : OTU table; sample metadata
# Output  : joined OTU table + sample_metadata_current.tsv
# Run on  : Atlas9 (Python 3)
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================

import csv
from pathlib import Path
from collections import Counter

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

INDIR = (
    RUN /
    "07_analysis/MC_GG2_97_classified_OTU_table"
)

OTU_TABLE = (
    INDIR /
    "MC_GG2_97_classified_576_samples_x_328_OTUs.tsv"
)

METADATA = (
    RUN /
    "metadata/sample_metadata_current.tsv"
)

OUTDIR = (
    RUN /
    "07_analysis/MC_GG2_97_analysis_ready"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

OUTFILE = (
    OUTDIR /
    "MC_GG2_97_classified_576_samples_x_328_OTUs_with_metadata.tsv"
)

METADATA_OUT = (
    OUTDIR /
    "MC_GG2_97_576_sample_metadata.tsv"
)


# ============================================================
# 1. READ METADATA
# ============================================================

metadata = {}
metadata_order = []

with open(METADATA, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    metadata_fields = reader.fieldnames

    if "Sample_ID" not in metadata_fields:
        raise RuntimeError(
            "Metadata does not contain Sample_ID column"
        )

    for row in reader:

        sid = row["Sample_ID"].strip()

        if sid in metadata:
            raise RuntimeError(
                f"Duplicate metadata Sample_ID: {sid}"
            )

        metadata[sid] = row
        metadata_order.append(sid)


print("=" * 76)
print("METADATA")
print("=" * 76)

print("Metadata rows:", len(metadata))
print("Unique Sample_IDs:", len(set(metadata)))

if len(metadata) != 576:
    raise RuntimeError(
        f"Expected 576 metadata samples, found {len(metadata)}"
    )

print("576 unique metadata samples: PASS")


# ============================================================
# 2. READ MC GG2-97 CLASSIFIED ABUNDANCE TABLE
# ============================================================

otu_rows = {}
otu_sample_order = []

with open(OTU_TABLE, newline="") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    fields = reader.fieldnames

    if not fields:
        raise RuntimeError("OTU table has no header")

    if fields[0] != "sample_id":
        raise RuntimeError(
            f"Expected first OTU-table column to be sample_id, "
            f"found {fields[0]}"
        )

    otu_columns = fields[1:]

    for row in reader:

        sid = row["sample_id"].strip()

        if sid in otu_rows:
            raise RuntimeError(
                f"Duplicate sample_id in OTU table: {sid}"
            )

        # Validate every abundance value.
        values = {}

        for otu in otu_columns:

            raw = row[otu].strip()

            try:
                value = int(raw)
            except Exception:
                raise RuntimeError(
                    f"Non-integer abundance: "
                    f"sample={sid}, OTU={otu}, value={raw!r}"
                )

            if value < 0:
                raise RuntimeError(
                    f"Negative abundance: "
                    f"sample={sid}, OTU={otu}, value={value}"
                )

            values[otu] = value

        otu_rows[sid] = values
        otu_sample_order.append(sid)


print()
print("=" * 76)
print("MC GG2 >=97% CLASSIFIED OTU TABLE")
print("=" * 76)

print("Samples:", len(otu_rows))
print("OTU columns:", len(otu_columns))
print("Unique OTU columns:", len(set(otu_columns)))

if len(otu_rows) != 576:
    raise RuntimeError(
        f"Expected 576 OTU-table samples, found {len(otu_rows)}"
    )

if len(otu_columns) != 328:
    raise RuntimeError(
        f"Expected 328 classified OTUs, found {len(otu_columns)}"
    )

if len(set(otu_columns)) != 328:
    raise RuntimeError(
        "Duplicate OTU columns detected"
    )

if not all(x.startswith("OTU_") for x in otu_columns):
    raise RuntimeError(
        "Unexpected non-OTU column among abundance columns"
    )

print("576 samples x 328 OTUs: PASS")
print("All abundance values are non-negative integers: PASS")


# ============================================================
# 3. EXACT SAMPLE-ID MATCH
# ============================================================

metadata_ids = set(metadata)
otu_ids = set(otu_rows)

missing_from_metadata = sorted(
    otu_ids - metadata_ids
)

missing_from_otu = sorted(
    metadata_ids - otu_ids
)


print()
print("=" * 76)
print("SAMPLE-ID MATCH")
print("=" * 76)

print(
    "OTU-table samples missing from metadata:",
    len(missing_from_metadata)
)

print(
    "Metadata samples missing from OTU table:",
    len(missing_from_otu)
)

if missing_from_metadata:
    print("Examples:", missing_from_metadata[:20])

if missing_from_otu:
    print("Examples:", missing_from_otu[:20])

if missing_from_metadata or missing_from_otu:
    raise RuntimeError(
        "Metadata and OTU-table Sample_ID sets do not match"
    )

print("Exact 576-sample ID match: PASS")


# ============================================================
# 4. ABUNDANCE VALIDATION
# ============================================================

sample_totals = {}

for sid in otu_sample_order:

    sample_totals[sid] = sum(
        otu_rows[sid].values()
    )

total_abundance = sum(sample_totals.values())

zero_samples = [
    sid
    for sid, total in sample_totals.items()
    if total == 0
]


print()
print("=" * 76)
print("ABUNDANCE VALIDATION")
print("=" * 76)

print(
    "Total supporting-read abundance:",
    total_abundance
)

print(
    "Samples with zero GG97-classified abundance:",
    len(zero_samples)
)

if total_abundance != 1432996:
    raise RuntimeError(
        f"Expected total classified abundance 1,432,996, "
        f"found {total_abundance}"
    )

if zero_samples:
    raise RuntimeError(
        "At least one sample has zero classified abundance"
    )

print(
    "Total abundance = 1,432,996: PASS"
)

print(
    "All 576 samples contain >=1 classified OTU: PASS"
)


# ============================================================
# 5. METADATA SUMMARY
# ============================================================

print()
print("=" * 76)
print("METADATA SUMMARY")
print("=" * 76)

for field in [
    "Source",
    "Plate",
    "Origin",
    "Location_Status",
    "Include_in_Location_Analysis",
]:

    if field not in metadata_fields:
        continue

    counts = Counter(
        (
            metadata[sid][field].strip()
            if metadata[sid][field].strip()
            else "<blank>"
        )
        for sid in otu_sample_order
    )

    print()
    print(field)

    for value, n in sorted(
        counts.items(),
        key=lambda x: (-x[1], x[0])
    ):
        print(f"  {value}: {n}")


# ============================================================
# 6. WRITE METADATA IN EXACT OTU-TABLE SAMPLE ORDER
# ============================================================

with open(METADATA_OUT, "w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=metadata_fields,
        delimiter="\t",
        extrasaction="ignore"
    )

    writer.writeheader()

    for sid in otu_sample_order:
        writer.writerow(metadata[sid])


# ============================================================
# 7. WRITE JOINED ANALYSIS-READY TABLE
#
# Metadata columns first, then 328 OTUs.
# ============================================================

joined_fields = (
    metadata_fields +
    otu_columns
)

with open(OUTFILE, "w", newline="") as fh:

    writer = csv.writer(
        fh,
        delimiter="\t"
    )

    writer.writerow(joined_fields)

    for sid in otu_sample_order:

        meta = metadata[sid]

        writer.writerow(
            [
                meta[field]
                for field in metadata_fields
            ]
            +
            [
                otu_rows[sid][otu]
                for otu in otu_columns
            ]
        )


# ============================================================
# 8. FINAL OUTPUT CHECK
# ============================================================

with open(OUTFILE, newline="") as fh:

    reader = csv.reader(
        fh,
        delimiter="\t"
    )

    header = next(reader)
    rows = list(reader)


print()
print("=" * 76)
print("FINAL ANALYSIS-READY TABLE")
print("=" * 76)

print("Rows / samples:", len(rows))
print("Total columns:", len(header))
print("Metadata columns:", len(metadata_fields))
print("OTU columns:", len(otu_columns))

if len(rows) != 576:
    raise RuntimeError(
        "Final joined table does not contain 576 samples"
    )

if len(header) != (
    len(metadata_fields) + 328
):
    raise RuntimeError(
        "Final joined table column count incorrect"
    )

print("Joined table dimensions: PASS")

print()
print("OUTPUT FILES")

print()
print("Analysis metadata:")
print(METADATA_OUT)

print()
print("Metadata + 328 classified OTUs:")
print(OUTFILE)

print()
print("=" * 76)
print("MC GG2 >=97% METADATA JOIN COMPLETE")
print("=" * 76)

