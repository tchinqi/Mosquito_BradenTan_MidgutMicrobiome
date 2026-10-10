#!/usr/bin/env python3

import csv
import shutil
from pathlib import Path
from collections import Counter

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

META = RUN / "metadata/sample_metadata_current.tsv"

BACKUP = (
    RUN /
    "metadata/sample_metadata_before_6control_fix.tsv"
)

TMP = (
    RUN /
    "metadata/sample_metadata_current.tmp.tsv"
)

TARGETS = {
    "sample_BT_CLX_096_16S_BT_P01_16S": "P01",
    "sample_BT_CLX_120_16S_BT_P02_16S": "P02",
}


# ============================================================
# BACKUP
# ============================================================

if not BACKUP.exists():
    shutil.copy2(META, BACKUP)

    print("Backup created:")
    print(BACKUP)

else:
    print("Backup already exists:")
    print(BACKUP)


# ============================================================
# READ + CORRECT METADATA
# ============================================================

rows = []
found = set()

with open(META, newline="", encoding="utf-8") as fh:

    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames

    for row in reader:

        sid = row["Sample_ID"].strip()

        if sid in TARGETS:

            found.add(sid)

            if row["Source"] != "BT":
                raise RuntimeError(
                    f"{sid}: expected Source=BT, "
                    f"found {row['Source']}"
                )

            if row["Plate"] != TARGETS[sid]:
                raise RuntimeError(
                    f"{sid}: expected Plate={TARGETS[sid]}, "
                    f"found {row['Plate']}"
                )

            row["Origin"] = "Control"

            row["Location_Code"] = ""
            row["Location_Name"] = ""

            row["Location_Status"] = (
                "Not applicable — negative control"
            )

            row["Include_in_Location_Analysis"] = "NO"

            row["Exclusion_or_Review_Reason"] = (
                "Negative control; final well of BT plate."
            )

            # User confirmed final physical well of each
            # 96-well BT plate was the negative control.
            row["Source_Well"] = "H12"

        rows.append(row)


if found != set(TARGETS):
    raise RuntimeError(
        "Did not find both BT control samples. "
        f"Found: {sorted(found)}"
    )


# ============================================================
# WRITE CORRECTED METADATA
# ============================================================

with open(TMP, "w", newline="", encoding="utf-8") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

TMP.replace(META)


# ============================================================
# VALIDATE
# ============================================================

origin_counts = Counter(
    row["Origin"]
    for row in rows
)

controls = [
    row
    for row in rows
    if row["Origin"] == "Control"
]


print()
print("=" * 72)
print("CORRECTED SAMPLE COUNTS")
print("=" * 72)

print("Wild:", origin_counts["Wild"])
print("Lab-reared:", origin_counts["Lab-reared"])
print("Control:", origin_counts["Control"])
print("Total:", len(rows))

print()
print("===== SIX CONTROLS =====")

for row in controls:

    print(
        row["Sample_ID"],
        "| Source =", row["Source"],
        "| Plate =", row["Plate"],
        "| Well =", row["Source_Well"]
    )


if len(rows) != 576:
    raise RuntimeError(
        f"Expected 576 total samples, found {len(rows)}"
    )

if origin_counts["Wild"] != 380:
    raise RuntimeError(
        f"Expected 380 Wild, found "
        f"{origin_counts['Wild']}"
    )

if origin_counts["Lab-reared"] != 190:
    raise RuntimeError(
        f"Expected 190 Lab-reared, found "
        f"{origin_counts['Lab-reared']}"
    )

if origin_counts["Control"] != 6:
    raise RuntimeError(
        f"Expected 6 Controls, found "
        f"{origin_counts['Control']}"
    )

print()
print(
    "380 Wild + 190 Lab-reared + "
    "6 Controls = 576: PASS"
)

print()
print("Canonical metadata updated:")
print(META)

