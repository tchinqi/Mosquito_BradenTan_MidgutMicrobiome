#!/usr/bin/env python3

import csv
from pathlib import Path

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

# ------------------------------------------------------------
# FINAL CORRECTED DATASET DEFINITIONS
#
# LC = 0.1%                 = 7,642 OTUs
# MC = max(0.5%, 20 reads)  =   867 OTUs
# ------------------------------------------------------------

DATASETS = {

    "LC": {
        "expected_raw": 7642,
        "expected_unassigned": 4870,

        "taxonomy": (
            RUN /
            "06_taxonomy/greengenes2_corrected_names/"
            "LC_0p1pct/default80/exported_taxonomy/taxonomy.tsv"
        ),

        "fasta": (
            RUN /
            "04_global_vsearch97/postracon_thresholds_576/"
            "vsearch97/0p1pct/"
            "otu97_centroids_UPPERCASE.fasta"
        ),
    },

    "MC": {
        "expected_raw": 867,
        "expected_unassigned": 186,

        "taxonomy": (
            RUN /
            "06_taxonomy/greengenes2_corrected_names/"
            "MC_max0p5pct_20reads/default80/"
            "exported_taxonomy/taxonomy.tsv"
        ),

        "fasta": (
            RUN /
            "04_global_vsearch97/postracon_thresholds_576/"
            "LC_max0p5pct_20reads/vsearch97/"
            "otu97_centroids_UPPERCASE.fasta"
        ),
    },
}

OUTDIR = (
    RUN /
    "06_taxonomy/kraken2_final/"
    "GG2_default80_unassigned"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# FASTA READER
# ------------------------------------------------------------

def read_fasta(path):

    records = {}
    current_id = None
    seq_parts = []

    with open(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if current_id is not None:
                    records[current_id] = "".join(seq_parts)

                # Feature ID = first item after ">"
                current_id = line[1:].split()[0]
                seq_parts = []

            else:
                seq_parts.append(line)

        if current_id is not None:
            records[current_id] = "".join(seq_parts)

    return records


# ------------------------------------------------------------
# PROCESS LC + MC
# ------------------------------------------------------------

for dataset, cfg in DATASETS.items():

    taxonomy_path = cfg["taxonomy"]
    fasta_path = cfg["fasta"]

    print()
    print("=" * 70)
    print(dataset)
    print("=" * 70)

    print("Taxonomy:")
    print(taxonomy_path)

    print("FASTA:")
    print(fasta_path)

    # --------------------------------------------------------
    # Read GG2 taxonomy
    # --------------------------------------------------------

    all_taxonomy_ids = set()
    unassigned_ids = []

    with open(taxonomy_path, newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            feature_id = row["Feature ID"].strip()
            taxon = (row["Taxon"] or "").strip()

            all_taxonomy_ids.add(feature_id)

            # Priya requested completely unassigned OTUs only.
            # Do NOT include OTUs assigned above family,
            # or genus/species-level "unclassified" OTUs.
            if taxon.lower() == "unassigned":
                unassigned_ids.append(feature_id)

    # --------------------------------------------------------
    # Validate taxonomy counts
    # --------------------------------------------------------

    print()
    print("GG2 default/80% taxonomy OTUs:",
          len(all_taxonomy_ids))

    print("Completely Unassigned OTUs:",
          len(unassigned_ids))

    if len(all_taxonomy_ids) != cfg["expected_raw"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected_raw']} total taxonomy "
            f"records but found {len(all_taxonomy_ids)}"
        )

    if len(unassigned_ids) != cfg["expected_unassigned"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected_unassigned']} "
            f"Unassigned OTUs but found {len(unassigned_ids)}"
        )

    if len(unassigned_ids) != len(set(unassigned_ids)):
        raise RuntimeError(
            f"{dataset}: duplicate Feature IDs detected "
            f"in Unassigned list"
        )

    # --------------------------------------------------------
    # Read complete OTU FASTA
    # --------------------------------------------------------

    fasta_records = read_fasta(fasta_path)

    print("OTUs in source FASTA:",
          len(fasta_records))

    if len(fasta_records) != cfg["expected_raw"]:
        raise RuntimeError(
            f"{dataset}: source FASTA contains "
            f"{len(fasta_records)} records; "
            f"expected {cfg['expected_raw']}"
        )

    # --------------------------------------------------------
    # Verify every GG2-unassigned ID exists in source FASTA
    # --------------------------------------------------------

    missing = [
        feature_id
        for feature_id in unassigned_ids
        if feature_id not in fasta_records
    ]

    print("GG2-unassigned IDs missing from FASTA:",
          len(missing))

    if missing:
        print("First missing IDs:", missing[:20])

        raise RuntimeError(
            f"{dataset}: cannot continue because "
            f"{len(missing)} GG2-unassigned OTUs "
            f"are absent from source FASTA"
        )

    # --------------------------------------------------------
    # Write ID list
    # --------------------------------------------------------

    ids_out = (
        OUTDIR /
        f"{dataset}_GG2_default80_unassigned_OTU_ids.txt"
    )

    with open(ids_out, "w") as fh:

        for feature_id in unassigned_ids:
            fh.write(feature_id + "\n")

    # --------------------------------------------------------
    # Write subset FASTA
    # --------------------------------------------------------

    fasta_out = (
        OUTDIR /
        f"{dataset}_GG2_default80_unassigned.fasta"
    )

    with open(fasta_out, "w") as fh:

        for feature_id in unassigned_ids:

            sequence = fasta_records[feature_id]

            fh.write(f">{feature_id}\n")

            # 80-character FASTA wrapping
            for i in range(0, len(sequence), 80):
                fh.write(sequence[i:i+80] + "\n")

    # --------------------------------------------------------
    # Re-read output as final validation
    # --------------------------------------------------------

    subset_records = read_fasta(fasta_out)

    output_ids = set(subset_records)

    expected_ids = set(unassigned_ids)

    if output_ids != expected_ids:
        raise RuntimeError(
            f"{dataset}: output FASTA IDs do not exactly match "
            f"GG2 Unassigned IDs"
        )

    if len(subset_records) != cfg["expected_unassigned"]:
        raise RuntimeError(
            f"{dataset}: output FASTA contains "
            f"{len(subset_records)} records; "
            f"expected {cfg['expected_unassigned']}"
        )

    print()
    print("FINAL VALIDATION")
    print(
        f"Expected Unassigned OTUs: "
        f"{cfg['expected_unassigned']}"
    )
    print(
        f"Extracted FASTA OTUs: "
        f"{len(subset_records)}"
    )
    print(
        "Exact Feature-ID match:",
        "PASS"
    )

    print()
    print("ID list:")
    print(ids_out)

    print("FASTA for Kraken2:")
    print(fasta_out)


print()
print("=" * 70)
print("EXTRACTION COMPLETE")
print("=" * 70)

print()
print("Expected final FASTAs:")
print(
    "LC = 4,870 GG2-default/80%-Unassigned OTUs"
)
print(
    "MC = 186 GG2-default/80%-Unassigned OTUs"
)

print()
print("Output directory:")
print(OUTDIR)

