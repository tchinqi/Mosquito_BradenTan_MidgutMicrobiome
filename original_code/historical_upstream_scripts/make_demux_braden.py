from pathlib import Path
import pandas as pd
import re
import sys

# ============================================================
# Braden FYP - 16S Minibar demultiplexing sheet preparation
# ============================================================

ROOT = Path("/hpctmp/e0958567/Braden_16S_20260819")

EXCEL = ROOT / "00_source" / "Copy_of_Braden_Tan_Demultiplexing_Sheet.xlsx"
OUTDIR = ROOT / "01_demux_metadata"

OUTDIR.mkdir(parents=True, exist_ok=True)

# 16S primers used in this run
F_PRIMER = "ACTCCTACGGGAGGCAGCAG"
R_PRIMER = "GGACTACHVGGGTWTCTAAT"

EXPECTED_SAMPLES = 768

# Only the eight 16S plates from this sequencing run
SHEETS = {
    "AG_Demultiplexing_16S_P01": "AG_P01_16S",
    "AG_Demultiplexing_16S_P02": "AG_P02_16S",
    "ZX_Demultiplexing_16S_P01": "ZX_P01_16S",
    "ZX_Demultiplexing_16S_P02": "ZX_P02_16S",
    "BT_Demultiplexing_16S_P01": "BT_P01_16S",
    "BT_Demultiplexing_16S_P02": "BT_P02_16S",
    "AG_Reproductive_Organs_16S_P01": "AG_Repro_P01_16S",
    "AG_Reproductive_Organs_16s_P02": "AG_Repro_P02_16S",
}


def find_column(df, possibilities):
    for col in possibilities:
        if col in df.columns:
            return col
    return None


def clean_name(x):
    """Make sample names safe for downstream filenames."""
    x = str(x).strip()
    x = re.sub(r"[\s.:]+", "_", x)
    return x


# ------------------------------------------------------------
# Check workbook
# ------------------------------------------------------------

if not EXCEL.exists():
    print(f"ERROR: Cannot find workbook: {EXCEL}")
    sys.exit(1)

xls = pd.ExcelFile(EXCEL)

missing = [s for s in SHEETS if s not in xls.sheet_names]

if missing:
    print("ERROR: Missing expected sheets:")
    for s in missing:
        print("  ", s)
    sys.exit(1)


# ------------------------------------------------------------
# Read all plates
# ------------------------------------------------------------

records = []
excluded = []

for sheet, plate_id in SHEETS.items():

    print(f"Reading {sheet}")

    df = pd.read_excel(
        EXCEL,
        sheet_name=sheet,
        dtype=str,
        keep_default_na=False
    )

    df.columns = [str(c).strip() for c in df.columns]

    specimen_col = find_column(
        df,
        [
            "Specimen_code_number",
            "Specimen-code-number",
            "Specimen code number"
        ]
    )

    well_col = find_column(
        df,
        [
            "Well no.",
            "Well no",
            "Well"
        ]
    )

    if specimen_col is None:
        print(f"ERROR: Cannot find specimen column in {sheet}")
        print(df.columns.tolist())
        sys.exit(1)

    if "F_PRIMER_TAG" not in df.columns or "R_PRIMER_TAG" not in df.columns:
        print(f"ERROR: Primer-tag columns missing from {sheet}")
        sys.exit(1)

    for idx, row in df.iterrows():

        specimen = str(row.get(specimen_col, "")).strip()

        ftag = str(row.get("F_PRIMER_TAG", "")).strip().upper()
        rtag = str(row.get("R_PRIMER_TAG", "")).strip().upper()

        well = ""
        if well_col is not None:
            well = str(row.get(well_col, "")).strip()

        # Ignore blank / incomplete rows
        if not specimen or not ftag or not rtag:

            excluded.append({
                "Sheet": sheet,
                "Excel_row": idx + 2,
                "Specimen": specimen,
                "F_PRIMER_TAG": ftag,
                "R_PRIMER_TAG": rtag,
                "Reason": "Incomplete specimen or primer-tag information"
            })

            continue

        # F_PRIMER_TAG should be:
        # 9-bp forward barcode + forward primer
        if not ftag.endswith(F_PRIMER):

            print("\nERROR: Forward primer suffix mismatch")
            print("Sheet:", sheet)
            print("Excel row:", idx + 2)
            print("Specimen:", specimen)
            print("F_PRIMER_TAG:", ftag)

            sys.exit(1)

        # R_PRIMER_TAG should be:
        # 9-bp reverse barcode + reverse primer
        if not rtag.endswith(R_PRIMER):

            print("\nERROR: Reverse primer suffix mismatch")
            print("Sheet:", sheet)
            print("Excel row:", idx + 2)
            print("Specimen:", specimen)
            print("R_PRIMER_TAG:", rtag)

            sys.exit(1)

        # Separate barcode from primer
        fw_index = ftag[:-len(F_PRIMER)]
        rv_index = rtag[:-len(R_PRIMER)]

        if len(fw_index) != 9 or len(rv_index) != 9:

            print("\nERROR: Barcode is not 9 bp")
            print("Sheet:", sheet)
            print("Specimen:", specimen)
            print("Forward barcode:", fw_index, len(fw_index))
            print("Reverse barcode:", rv_index, len(rv_index))

            sys.exit(1)

        sample_id = f"{clean_name(specimen)}_{plate_id}"

        records.append({
            "SampleID": sample_id,
            "Specimen": specimen,
            "Plate_ID": plate_id,
            "Sheet": sheet,
            "Excel_row": idx + 2,
            "Well": well,

            "FwIndex": fw_index,
            "FwPrimer": F_PRIMER,
            "RvIndex": rv_index,
            "RvPrimer": R_PRIMER,

            "Original_F_PRIMER_TAG": ftag,
            "Original_R_PRIMER_TAG": rtag
        })


data = pd.DataFrame(records)

print()
print("Usable rows before duplicate handling:", len(data))


# ------------------------------------------------------------
# Handle duplicated SampleIDs
# ------------------------------------------------------------

dup_mask = data.duplicated("SampleID", keep=False)

if dup_mask.any():

    print("Duplicate SampleIDs detected:")
    print(
        data.loc[
            dup_mask,
            ["SampleID", "Sheet", "Well", "Excel_row"]
        ].to_string(index=False)
    )

    print()
    print("Adding well number to duplicated SampleIDs...")

    for idx in data[dup_mask].index:

        well = str(data.loc[idx, "Well"]).strip()

        if not well:
            print("ERROR: Duplicate SampleID has no well information.")
            print(data.loc[idx])
            sys.exit(1)

        data.loc[idx, "SampleID"] = (
            data.loc[idx, "SampleID"]
            + "_"
            + clean_name(well)
        )


# Check duplicates again
remaining_dupes = data[data.duplicated("SampleID", keep=False)]

if not remaining_dupes.empty:

    remaining_dupes.to_csv(
        OUTDIR / "duplicate_sampleIDs.csv",
        index=False
    )

    print("ERROR: Duplicate SampleIDs remain.")
    sys.exit(1)


# ------------------------------------------------------------
# Check barcode-pair clashes
# ------------------------------------------------------------

clashes = data[
    data.duplicated(
        ["FwIndex", "RvIndex"],
        keep=False
    )
]

if not clashes.empty:

    clashes.to_csv(
        OUTDIR / "barcode_clashes.csv",
        index=False
    )

    print("ERROR: Duplicate forward/reverse barcode combinations found.")
    sys.exit(1)


# ------------------------------------------------------------
# Expected sample number
# ------------------------------------------------------------

if len(data) != EXPECTED_SAMPLES:

    print(
        f"ERROR: Expected {EXPECTED_SAMPLES} usable samples "
        f"but found {len(data)}"
    )

    sys.exit(1)


# ------------------------------------------------------------
# Write Minibar input
# ------------------------------------------------------------

minibar = data[
    [
        "SampleID",
        "FwIndex",
        "FwPrimer",
        "RvIndex",
        "RvPrimer"
    ]
].copy()

demux_file = OUTDIR / "Demux_sheet_minibar.txt"

minibar.to_csv(
    demux_file,
    sep="\t",
    index=False
)


# Full metadata audit trail
data.to_csv(
    OUTDIR / "demux_metadata_audit.csv",
    index=False
)

pd.DataFrame(excluded).to_csv(
    OUTDIR / "excluded_rows.csv",
    index=False
)


# ------------------------------------------------------------
# Final QC
# ------------------------------------------------------------

print()
print("========================================")
print("DEMULTIPLEXING SHEET QC")
print("========================================")

print("Usable samples:", len(data))
print("Unique SampleIDs:", data["SampleID"].nunique())

print(
    "Unique barcode pairs:",
    data[["FwIndex", "RvIndex"]]
    .drop_duplicates()
    .shape[0]
)

print(
    "Forward barcode lengths:",
    sorted(data["FwIndex"].str.len().unique())
)

print(
    "Reverse barcode lengths:",
    sorted(data["RvIndex"].str.len().unique())
)

print("Excluded/incomplete rows:", len(excluded))

print()
print("Samples per plate:")

print(
    data["Plate_ID"]
    .value_counts()
    .sort_index()
    .to_string()
)

print()
print("PASS: Demultiplexing sheet generated successfully.")
print("Output:", demux_file)
