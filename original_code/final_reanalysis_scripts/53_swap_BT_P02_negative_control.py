#!/usr/bin/env python3
"""
53_swap_BT_P02_negative_control.py

Correction to script 49. According to the demultiplexing sheet
(01_demux_metadata/demux_metadata_audit.csv), BT specimen numbers follow
well order, and the final well (H12) of each BT plate was the negative
control:
    BT_P01 H12 = BT_CLX_096  (already a control - unchanged)
    BT_P02 H12 = BT_ALB_192  (was wrongly treated as a biological sample)
BT_CLX_120 is in well H3 of BT_P02 and is a laboratory Culex sample that
script 49 wrongly assumed was the final well.

This script:
  1. restores BT_CLX_120 to its pre-script-49 metadata (Lab sample, well H3)
  2. marks BT_ALB_192 as the BT_P02 negative control (well H12)
Afterwards re-run scripts 45, 50, 51 and 52 (52 updated with the new control).
"""

import csv
import shutil
from pathlib import Path
from collections import Counter

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")

META        = RUN / "metadata/sample_metadata_current.tsv"
PRE49       = RUN / "metadata/sample_metadata_before_6control_fix.tsv"   # made by script 49
BACKUP      = RUN / "metadata/sample_metadata_before_BTP02_swap.tsv"
TMP         = RUN / "metadata/sample_metadata_current.tmp.tsv"

RESTORE = "sample_BT_CLX_120_16S_BT_P02_16S"   # back to biological (Lab Culex)
NEWCTRL = "sample_BT_ALB_192_16S_BT_P02_16S"   # becomes BT_P02 negative control

# ---- backup -----------------------------------------------------------------
if not BACKUP.exists():
    shutil.copy2(META, BACKUP)
    print("Backup created:", BACKUP)
else:
    print("Backup already exists:", BACKUP)

if not PRE49.exists():
    raise SystemExit(f"Cannot find pre-script-49 backup: {PRE49}")

with open(PRE49, newline="", encoding="utf-8") as fh:
    pre49 = {r["Sample_ID"].strip(): r for r in csv.DictReader(fh, delimiter="\t")}

if RESTORE not in pre49:
    raise SystemExit(f"{RESTORE} not found in {PRE49}")
if pre49[RESTORE]["Origin"].strip().lower().startswith("control"):
    raise SystemExit(f"{RESTORE} is already a control in the pre-49 backup - stopping")

# ---- correct ----------------------------------------------------------------
rows, found = [], set()
with open(META, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    fields = reader.fieldnames
    for row in reader:
        sid = row["Sample_ID"].strip()

        if sid == RESTORE:
            found.add(sid)
            row = {k: pre49[RESTORE].get(k, "") for k in fields}
            if "Source_Well" in fields:
                row["Source_Well"] = "H3"
            if "Exclusion_or_Review_Reason" in fields:
                row["Exclusion_or_Review_Reason"] = (
                    "Restored by script 53: laboratory Culex sample in well H3 "
                    "(previously misassigned as BT_P02 negative control)."
                )

        elif sid == NEWCTRL:
            found.add(sid)
            if row.get("Source", "BT") != "BT" or row.get("Plate", "P02") != "P02":
                raise SystemExit(f"{sid}: unexpected Source/Plate {row.get('Source')}/{row.get('Plate')}")
            row["Origin"] = "Control"
            for k, v in [("Location_Code", ""), ("Location_Name", ""),
                         ("Location_Status", "Not applicable — negative control"),
                         ("Include_in_Location_Analysis", "NO"),
                         ("Exclusion_or_Review_Reason",
                          "Negative control; final well (H12) of BT_P02 per demultiplexing sheet."),
                         ("Source_Well", "H12")]:
                if k in fields:
                    row[k] = v

        rows.append(row)

if found != {RESTORE, NEWCTRL}:
    raise SystemExit(f"Did not find both samples. Found: {sorted(found)}")

with open(TMP, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
TMP.replace(META)

# ---- validate ---------------------------------------------------------------
def group(o):
    o = o.strip().lower()
    return "Wild" if o.startswith("wild") else "Lab" if o.startswith("lab") else "Control"

counts = Counter(group(r["Origin"]) for r in rows)
controls = sorted(r["Sample_ID"] for r in rows if group(r["Origin"]) == "Control")
print("\nOrigin counts:", dict(counts))
print("Controls:")
for c in controls:
    print("  ", c)

ok = counts == Counter({"Wild": 380, "Lab": 190, "Control": 6}) and NEWCTRL in controls and RESTORE not in controls
print("\nVALIDATION:", "PASS" if ok else "FAIL - check above")
if not ok:
    raise SystemExit(1)
print("\nNext: re-run scripts 45, 50, 51, then 52.")

