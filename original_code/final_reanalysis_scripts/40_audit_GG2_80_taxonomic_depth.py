#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

FILES = {
    "LC": (
        RUN /
        "06_taxonomy/greengenes2_corrected_names/"
        "LC_0p1pct/default80/exported_taxonomy/taxonomy.tsv"
    ),
    "MC": (
        RUN /
        "06_taxonomy/greengenes2_corrected_names/"
        "MC_max0p5pct_20reads/default80/"
        "exported_taxonomy/taxonomy.tsv"
    ),
}


def deepest_rank(taxon):

    if not taxon or taxon.strip().lower() == "unassigned":
        return "Unassigned"

    ranks = [
        ("s__", "Species"),
        ("g__", "Genus"),
        ("f__", "Family"),
        ("o__", "Order"),
        ("c__", "Class"),
        ("p__", "Phylum"),
        ("d__", "Domain"),
    ]

    parts = [x.strip() for x in taxon.split(";")]

    for prefix, label in ranks:

        for part in parts:

            if part.startswith(prefix):
                value = part[len(prefix):].strip()

                if value:
                    return label

    return "Other / unresolved"


for dataset, path in FILES.items():

    counts = Counter()
    domain_only_strings = Counter()

    with open(path, newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            taxon = (row["Taxon"] or "").strip()
            depth = deepest_rank(taxon)

            counts[depth] += 1

            if depth == "Domain":
                domain_only_strings[taxon] += 1

    print()
    print("=" * 70)
    print(dataset)
    print("=" * 70)

    for category in [
        "Unassigned",
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Other / unresolved",
    ]:
        print(f"{category:20s}: {counts[category]}")

    print()
    print("Domain-only taxonomy strings:")

    if domain_only_strings:
        for taxon, n in domain_only_strings.items():
            print(f"  {n:6d}  {taxon}")
    else:
        print("  NONE")

    print()
    print("Total:", sum(counts.values()))

