#!/usr/bin/env python3
# =============================================================================
# Script 41 - 41_extract_final_GG80_rescue_FASTAs.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Kraken2 cross-check
# Purpose : Extracts OTUs unassigned or assigned only to d__Bacteria by
#           Greengenes2 at 80% (LC 4,884; MC 188).
# Input   : default80 taxonomy.tsv; OTU centroid FASTAs
# Output  : FINAL_GG80_unassigned_domain_rescue/<LC|MC>_GG80_Unassigned_plus_d
#           omain_only.fasta
# Run on  : Atlas9 (Python 3)
# Note    : Historical file/folder names swap LC and MC: '0p1pct' /
#           'MC_0p1pct' = final LC (0.1%); 'LC_max0p5pct_20reads' = final MC
#           (max(0.5%, 20 reads)).
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================

import csv
from pathlib import Path

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

DATASETS = {
    "LC": {
        "expected_raw": 7642,
        "expected_unassigned": 4870,
        "expected_domain_only": 14,
        "expected_final": 4884,

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
        "expected_domain_only": 2,
        "expected_final": 188,

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
    "FINAL_GG80_unassigned_domain_rescue"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


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

                current_id = line[1:].split()[0]
                seq_parts = []

            else:
                seq_parts.append(line)

        if current_id is not None:
            records[current_id] = "".join(seq_parts)

    return records


for dataset, cfg in DATASETS.items():

    print()
    print("=" * 70)
    print(dataset)
    print("=" * 70)

    unassigned_ids = []
    domain_only_ids = []
    all_ids = []

    # --------------------------------------------------------
    # Read GG2 80% taxonomy
    # --------------------------------------------------------

    with open(cfg["taxonomy"], newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            feature_id = row["Feature ID"].strip()
            taxon = (row["Taxon"] or "").strip()

            all_ids.append(feature_id)

            if taxon.lower() == "unassigned":
                unassigned_ids.append(feature_id)

            elif taxon == "d__Bacteria":
                domain_only_ids.append(feature_id)

    rescue_ids = unassigned_ids + domain_only_ids

    # --------------------------------------------------------
    # Validate GG2 counts
    # --------------------------------------------------------

    print("Raw GG2 OTUs:", len(all_ids))
    print("Unassigned:", len(unassigned_ids))
    print("Domain-only d__Bacteria:", len(domain_only_ids))
    print("Final rescue OTUs:", len(rescue_ids))

    if len(all_ids) != cfg["expected_raw"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected_raw']} raw OTUs, "
            f"found {len(all_ids)}"
        )

    if len(unassigned_ids) != cfg["expected_unassigned"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected_unassigned']} "
            f"Unassigned OTUs, found {len(unassigned_ids)}"
        )

    if len(domain_only_ids) != cfg["expected_domain_only"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected_domain_only']} "
            f"domain-only OTUs, found {len(domain_only_ids)}"
        )

    if len(rescue_ids) != cfg["expected_final"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected_final']} "
            f"final rescue OTUs, found {len(rescue_ids)}"
        )

    if len(rescue_ids) != len(set(rescue_ids)):
        raise RuntimeError(
            f"{dataset}: duplicate rescue Feature IDs detected"
        )

    # --------------------------------------------------------
    # Read source OTU FASTA
    # --------------------------------------------------------

    fasta = read_fasta(cfg["fasta"])

    if len(fasta) != cfg["expected_raw"]:
        raise RuntimeError(
            f"{dataset}: source FASTA expected "
            f"{cfg['expected_raw']} records, found {len(fasta)}"
        )

    missing = [x for x in rescue_ids if x not in fasta]

    print("Rescue IDs missing from FASTA:", len(missing))

    if missing:
        raise RuntimeError(
            f"{dataset}: {len(missing)} rescue OTUs missing "
            f"from source FASTA"
        )

    # --------------------------------------------------------
    # Write audit ID table
    # --------------------------------------------------------

    idfile = OUTDIR / f"{dataset}_GG80_rescue_OTU_ids.tsv"

    unassigned_set = set(unassigned_ids)

    with open(idfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Feature_ID",
            "GG2_80_category",
        ])

        for feature_id in rescue_ids:

            category = (
                "Unassigned"
                if feature_id in unassigned_set
                else "Domain_only_d__Bacteria"
            )

            w.writerow([
                feature_id,
                category,
            ])

    # --------------------------------------------------------
    # Write final rescue FASTA
    # --------------------------------------------------------

    outfile = (
        OUTDIR /
        f"{dataset}_GG80_Unassigned_plus_domain_only.fasta"
    )

    with open(outfile, "w") as fh:

        for feature_id in rescue_ids:

            seq = fasta[feature_id]

            fh.write(f">{feature_id}\n")

            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    # Re-read output
    check = read_fasta(outfile)

    if len(check) != cfg["expected_final"]:
        raise RuntimeError(
            f"{dataset}: final FASTA expected "
            f"{cfg['expected_final']} sequences, "
            f"found {len(check)}"
        )

    if set(check) != set(rescue_ids):
        raise RuntimeError(
            f"{dataset}: final FASTA Feature IDs "
            f"do not exactly match rescue IDs"
        )

    print()
    print("FINAL VALIDATION")
    print("Expected:", cfg["expected_final"])
    print("Extracted:", len(check))
    print("Exact Feature-ID match: PASS")

    print()
    print("FASTA:")
    print(outfile)

    print("Audit IDs:")
    print(idfile)


print()
print("=" * 70)
print("FINAL GG80 RESCUE FASTA EXTRACTION COMPLETE")
print("=" * 70)

print()
print("Expected:")
print("LC = 4,884 OTUs")
print("MC = 188 OTUs")

print()
print("Output directory:")
print(OUTDIR)

