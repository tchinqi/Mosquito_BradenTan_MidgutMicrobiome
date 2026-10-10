#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")
GG = RUN / "06_taxonomy/greengenes2_MC_LC_final"

SETS = {
    "MC": {
        "raw": 7642,
        "tax": GG / "MC_0p1pct/exported_taxonomy/taxonomy.tsv",
        "out": GG / "MC_0p1pct"
    },
    "LC": {
        "raw": 867,
        "tax": GG / "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv",
        "out": GG / "LC_max0p5pct_20reads"
    }
}

def extract_rank(taxonomy, prefix):
    for part in taxonomy.split(";"):
        part = part.strip()
        if part.startswith(prefix):
            value = part[len(prefix):].strip()
            if value:
                return value
    return ""

for label, info in SETS.items():

    genus_counts = Counter()
    species_counts = Counter()
    assigned_rows = []

    with open(info["tax"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            tax = (row["Taxon"] or "").strip()

            if not tax or tax.lower() == "unassigned":
                continue

            genus = extract_rank(tax, "g__")
            species = extract_rank(tax, "s__")

            if genus:
                genus_counts[genus] += 1

            if species:
                species_counts[species] += 1

            assigned_rows.append([
                row["Feature ID"],
                tax,
                genus if genus else "Unresolved",
                species if species else "Unresolved"
            ])

    assigned_n = len(assigned_rows)

    # Full assigned taxonomy table
    with open(info["out"] / "assigned_OTU_taxonomy.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Feature_ID", "Full_Taxonomy", "Genus", "Species"])
        w.writerows(assigned_rows)

    # Genus summary
    with open(info["out"] / "genus_counts.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "Genus",
            "OTU_count",
            "Percent_of_assigned_OTUs",
            "Percent_of_raw_OTUs"
        ])

        for genus, n in genus_counts.most_common():
            w.writerow([
                genus,
                n,
                "{:.2f}".format(100*n/assigned_n),
                "{:.2f}".format(100*n/info["raw"])
            ])

    # Species summary
    with open(info["out"] / "species_counts.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "Species",
            "OTU_count",
            "Percent_of_assigned_OTUs",
            "Percent_of_raw_OTUs"
        ])

        for species, n in species_counts.most_common():
            w.writerow([
                species,
                n,
                "{:.2f}".format(100*n/assigned_n),
                "{:.2f}".format(100*n/info["raw"])
            ])

    print()
    print("========== {} ==========".format(label))
    print("Assigned OTUs:", assigned_n)
    print("OTUs with genus-level assignment:", sum(genus_counts.values()))
    print("OTUs with species-level assignment:", sum(species_counts.values()))

    print()
    print("Top 20 genera:")
    for name, n in genus_counts.most_common(20):
        print("{:<45} {}".format(name, n))

    print()
    print("Top 20 species:")
    for name, n in species_counts.most_common(20):
        print("{:<55} {}".format(name, n))

