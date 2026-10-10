import pandas as pd
import re
from pathlib import Path

BASE = Path("/hpctmp/e0958567/Braden_16S_20260819")

otu_table = BASE / "13_global_otu97/otu97_biological_min500.tsv"
workbook = Path(
    "/hpctmp/e0958567/UG FYP Tan Chin Qi Braden/"
    "Copy_of_Braden_Tan_Demultiplexing_Sheet.xlsx"
)

out = BASE / "14_diversity/sample_metadata_549.tsv"

# ------------------------------------------------------------
# Get the exact 549 samples retained after QC
# ------------------------------------------------------------

with open(otu_table) as fh:
    header = fh.readline().rstrip("\n").split("\t")

samples = header[1:]

# ------------------------------------------------------------
# Read the six midgut plates from the workbook
# ------------------------------------------------------------

sheet_to_plate = {
    "AG_Demultiplexing_16S_P01": ("AG", "AG_P01"),
    "AG_Demultiplexing_16S_P02": ("AG", "AG_P02"),
    "ZX_Demultiplexing_16S_P01": ("ZX", "ZX_P01"),
    "ZX_Demultiplexing_16S_P02": ("ZX", "ZX_P02"),
    "BT_Demultiplexing_16S_P01": ("BT", "BT_P01"),
    "BT_Demultiplexing_16S_P02": ("BT", "BT_P02"),
}

records = {}

for sheet, (source, plate) in sheet_to_plate.items():

    df = pd.read_excel(workbook, sheet_name=sheet)

    # standardise column names
    df.columns = [
        str(c).strip().replace("\n", " ")
        for c in df.columns
    ]

    # locate specimen column
    specimen_col = None
    for c in df.columns:
        if c.lower().replace("_", " ") in [
            "specimen code number",
            "specimen-code-number"
        ]:
            specimen_col = c
            break

    if specimen_col is None:
        raise ValueError(
            f"Could not find specimen-code column in {sheet}"
        )

    for _, row in df.iterrows():

        specimen = row.get(specimen_col)

        if pd.isna(specimen):
            continue

        specimen = str(specimen).strip()

        if not specimen:
            continue

        tissue = str(row.get("Tissue Group", "")).strip()

        well = str(row.get("Well no.", "")).strip()

        # Store lookup using the specimen name exactly as it appears
        records[(source, plate, specimen)] = {
            "SpecimenCode": specimen,
            "Source": source,
            "Plate": plate,
            "Well": well,
            "Tissue": tissue
        }

# ------------------------------------------------------------
# Parse sample IDs
# ------------------------------------------------------------

rows = []
unmatched = []

for sample in samples:

    source = ""
    plate = ""
    specimen = ""

    # AG
    m = re.match(r"sample_(.+)_AG_P(01|02)_16S$", sample)
    if m:
        specimen = m.group(1)
        source = "AG"
        plate = f"AG_P{m.group(2)}"

    # ZX
    if not source:
        m = re.match(r"sample_(.+)_ZX_P(01|02)_16S$", sample)
        if m:
            specimen = m.group(1)
            source = "ZX"
            plate = f"ZX_P{m.group(2)}"

    # BT
    if not source:
        m = re.match(
            r"sample_(BT_(?:ALB|CLX)_\d+)_16S_BT_P(01|02)_16S$",
            sample
        )
        if m:
            specimen = m.group(1) + "_16S"
            source = "BT"
            plate = f"BT_P{m.group(2)}"

    key = (source, plate, specimen)

    info = records.get(key)

    if info is None:
        unmatched.append((sample, specimen, source, plate))

        rows.append({
            "#SampleID": sample,
            "SpecimenCode": specimen,
            "Source": source,
            "Plate": plate,
            "Well": "Unknown",
            "Tissue": "Midgut",
            "TaxonCode": "Unknown",
            "Origin": "Unknown",
            "Location": "Unknown",
        })

        continue

    # Only use ALB / CLX when it is explicitly encoded
    taxon = "Unknown"

    if specimen.startswith("BT_ALB_"):
        taxon = "ALB"

    elif specimen.startswith("BT_CLX_"):
        taxon = "CLX"

    rows.append({
        "#SampleID": sample,
        "SpecimenCode": info["SpecimenCode"],
        "Source": info["Source"],
        "Plate": info["Plate"],
        "Well": info["Well"],
        "Tissue": info["Tissue"],
        "TaxonCode": taxon,
        "Origin": "Unknown",
        "Location": "Unknown",
    })

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

meta = pd.DataFrame(rows)

meta.to_csv(out, sep="\t", index=False)

print("===== METADATA BUILD =====")
print("Samples in OTU table:", len(samples))
print("Metadata rows:", len(meta))
print("Matched to workbook:", len(samples) - len(unmatched))
print("Unmatched:", len(unmatched))

print()
print("Source counts:")
print(meta["Source"].value_counts(dropna=False).to_string())

print()
print("Plate counts:")
print(meta["Plate"].value_counts(dropna=False).to_string())

print()
print("TaxonCode currently known:")
print(meta["TaxonCode"].value_counts(dropna=False).to_string())

if unmatched:
    print()
    print("===== UNMATCHED =====")
    for x in unmatched[:30]:
        print(*x, sep="\t")

print()
print("Output:", out)
