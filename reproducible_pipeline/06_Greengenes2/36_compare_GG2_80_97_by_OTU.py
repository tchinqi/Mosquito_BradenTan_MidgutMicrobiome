#!/usr/bin/env python3
# =============================================================================
# Script 36 - 36_compare_GG2_80_97_by_OTU.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Greengenes2 comparison
# Purpose : Pairs each OTU's 80% and 97% assignments (unchanged / reassigned /
#           assigned->unassigned / unassigned->assigned); Table S2B.
# Input   : taxonomy.tsv at 80% and 97% (LC and MC)
# Output  : GG2_80_97_OTU_pairing/Greengenes2_<LC|MC>_80_vs_97_per_OTU.tsv
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

DATASETS = {
    "LC": {
        "expected": 7642,
        "tax80": (
            RUN /
            "06_taxonomy/greengenes2_corrected_names/"
            "LC_0p1pct/default80/exported_taxonomy/taxonomy.tsv"
        ),
        "tax97": (
            RUN /
            "06_taxonomy/greengenes2_MC_LC_final/"
            "MC_0p1pct/exported_taxonomy/taxonomy.tsv"
        ),
    },

    "MC": {
        "expected": 867,
        "tax80": (
            RUN /
            "06_taxonomy/greengenes2_corrected_names/"
            "MC_max0p5pct_20reads/default80/"
            "exported_taxonomy/taxonomy.tsv"
        ),
        "tax97": (
            RUN /
            "06_taxonomy/greengenes2_MC_LC_final/"
            "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv"
        ),
    },
}

OUTDIR = (
    RUN /
    "06_taxonomy/greengenes2_corrected_names/"
    "GG2_80_97_OTU_pairing"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


def read_taxonomy(path):

    records = {}

    with open(path, newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            feature_id = row["Feature ID"].strip()
            taxon = (row["Taxon"] or "").strip()

            records[feature_id] = taxon

    return records


def get_rank(taxonomy, prefix):

    if not taxonomy or taxonomy.lower() == "unassigned":
        return None

    for part in taxonomy.split(";"):

        part = part.strip()

        if part.startswith(prefix):

            value = part[len(prefix):].strip()

            if value:
                return value

    return None


def status(taxonomy):

    if not taxonomy or taxonomy.lower() == "unassigned":
        return "Unassigned"

    return "Assigned"


for dataset, cfg in DATASETS.items():

    print()
    print("=" * 70)
    print(dataset)
    print("=" * 70)

    tax80 = read_taxonomy(cfg["tax80"])
    tax97 = read_taxonomy(cfg["tax97"])

    ids80 = set(tax80)
    ids97 = set(tax97)

    print("80% OTUs:", len(ids80))
    print("97% OTUs:", len(ids97))

    if len(ids80) != cfg["expected"]:
        raise RuntimeError(
            f"{dataset}: 80% expected {cfg['expected']} "
            f"but found {len(ids80)}"
        )

    if len(ids97) != cfg["expected"]:
        raise RuntimeError(
            f"{dataset}: 97% expected {cfg['expected']} "
            f"but found {len(ids97)}"
        )

    if ids80 != ids97:

        print("Only in 80%:", len(ids80 - ids97))
        print("Only in 97%:", len(ids97 - ids80))

        raise RuntimeError(
            f"{dataset}: Feature ID sets differ"
        )

    print("Exact OTU-ID set match: PASS")

    counters = Counter()
    rows = []

    for feature_id in sorted(ids80):

        t80 = tax80[feature_id]
        t97 = tax97[feature_id]

        s80 = status(t80)
        s97 = status(t97)

        f80 = get_rank(t80, "f__")
        g80 = get_rank(t80, "g__")
        sp80 = get_rank(t80, "s__")

        f97 = get_rank(t97, "f__")
        g97 = get_rank(t97, "g__")
        sp97 = get_rank(t97, "s__")

        # Overall assignment transition
        if s80 == "Assigned" and s97 == "Assigned":

            if t80 == t97:
                transition = "ASSIGNED_SAME_FULL_TAXONOMY"
            else:
                transition = "ASSIGNED_TAXONOMY_CHANGED"

        elif s80 == "Assigned" and s97 == "Unassigned":
            transition = "ASSIGNED80_TO_UNASSIGNED97"

        elif s80 == "Unassigned" and s97 == "Assigned":
            transition = "UNASSIGNED80_TO_ASSIGNED97"

        else:
            transition = "UNASSIGNED_BOTH"

        counters[transition] += 1

        # Rank-by-rank comparison
        family_same = (
            f80 is not None
            and f97 is not None
            and f80 == f97
        )

        genus_same = (
            g80 is not None
            and g97 is not None
            and g80 == g97
        )

        species_same = (
            sp80 is not None
            and sp97 is not None
            and sp80 == sp97
        )

        rows.append([
            feature_id,
            s80,
            s97,
            transition,
            t80,
            t97,
            f80 or "",
            f97 or "",
            "YES" if family_same else "NO",
            g80 or "",
            g97 or "",
            "YES" if genus_same else "NO",
            sp80 or "",
            sp97 or "",
            "YES" if species_same else "NO",
        ])

    outfile = (
        OUTDIR /
        f"Greengenes2_{dataset}_80_vs_97_per_OTU.tsv"
    )

    with open(outfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Feature_ID",
            "Status_80",
            "Status_97",
            "Transition",
            "Taxonomy_80",
            "Taxonomy_97",
            "Family_80",
            "Family_97",
            "Same_family",
            "Genus_80",
            "Genus_97",
            "Same_genus",
            "Species_80",
            "Species_97",
            "Same_species",
        ])

        w.writerows(rows)

    print()
    print("OTU-level transitions:")

    for key in [
        "ASSIGNED_SAME_FULL_TAXONOMY",
        "ASSIGNED_TAXONOMY_CHANGED",
        "ASSIGNED80_TO_UNASSIGNED97",
        "UNASSIGNED80_TO_ASSIGNED97",
        "UNASSIGNED_BOTH",
    ]:
        print(f"  {key}: {counters[key]}")

    total = sum(counters.values())

    print()
    print(
        "Transition accounting:",
        total,
        "/",
        cfg["expected"],
        "PASS" if total == cfg["expected"] else "FAIL"
    )

    # Specifically inspect Acidomonas methanolica issue
    acidomonas = [
        r for r in rows
        if (
            "Acidomonas methanolica" in r[4]
            or "Acidomonas methanolica" in r[5]
        )
    ]

    print()
    print(
        "OTUs involving Acidomonas methanolica:",
        len(acidomonas)
    )

    for r in acidomonas[:20]:

        print()
        print("Feature ID:", r[0])
        print("80%:", r[4])
        print("97%:", r[5])
        print("Transition:", r[3])

    print()
    print("Audit file:")
    print(outfile)


print()
print("=" * 70)
print("GG2 80% vs 97% OTU PAIRING COMPLETE")
print("=" * 70)

print("Output directory:")
print(OUTDIR)

