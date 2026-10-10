#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

# ------------------------------------------------------------
# FINAL CORRECTED NAMING
#
# LC = 0.1% = 7,642 OTUs
# MC = max(0.5%, 20 reads) = 867 OTUs
#
# NOTE:
# Old 97% directories were named before Priya corrected MC/LC labels.
# ------------------------------------------------------------

FILES = {
    "LC": {
        "raw": 7642,

        # corrected LC = 0.1%
        "80": (
            RUN /
            "06_taxonomy/greengenes2_corrected_names/"
            "LC_0p1pct/default80/exported_taxonomy/taxonomy.tsv"
        ),

        # old MC directory = corrected LC
        "97": (
            RUN /
            "06_taxonomy/greengenes2_MC_LC_final/"
            "MC_0p1pct/exported_taxonomy/taxonomy.tsv"
        ),
    },

    "MC": {
        "raw": 867,

        # corrected MC = max(0.5%,20 reads)
        "80": (
            RUN /
            "06_taxonomy/greengenes2_corrected_names/"
            "MC_max0p5pct_20reads/default80/"
            "exported_taxonomy/taxonomy.tsv"
        ),

        # old LC directory = corrected MC
        "97": (
            RUN /
            "06_taxonomy/greengenes2_MC_LC_final/"
            "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv"
        ),
    },
}

OUTDIR = (
    RUN /
    "06_taxonomy/greengenes2_corrected_names/"
    "Priya_hierarchical_tables"
)

OUTDIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def read_taxonomy(path):

    records = {}

    with open(path, newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            feature = row["Feature ID"].strip()
            taxon = (row["Taxon"] or "").strip()

            records[feature] = taxon

    return records


def get_rank(taxonomy, prefix):

    for part in taxonomy.split(";"):

        part = part.strip()

        if part.startswith(prefix):

            value = part[len(prefix):].strip()

            if value:
                return value

    return None


def parse_lineage(taxonomy):

    # Completely unassigned at this identity threshold
    if not taxonomy or taxonomy.lower() == "unassigned":
        return None

    family = get_rank(taxonomy, "f__")
    genus = get_rank(taxonomy, "g__")
    species = get_rank(taxonomy, "s__")

    # If no family assignment was reached, keep it outside
    # the hierarchical Family–Genus–Species table.
    if not family:
        return ("NO_FAMILY", None, None)

    # Fill unresolved lower ranks explicitly
    if not genus:
        genus = "unclassified genus"
        species = "unclassified species"

    elif not species:
        species = "unclassified species"

    return (family, genus, species)


def summarise(records):

    lineages = Counter()

    unassigned = 0
    unclassified_below_family = 0

    reached_family = 0
    reached_genus = 0
    reached_species = 0

    for taxonomy in records.values():

        lineage = parse_lineage(taxonomy)

        if lineage is None:
            unassigned += 1
            continue

        family, genus, species = lineage

        # Assigned somewhere above family, but not to family
        if family == "NO_FAMILY":
            unclassified_below_family += 1
            continue

        # At least family reached
        reached_family += 1

        # "reached genus" excludes placeholder
        if genus != "unclassified genus":
            reached_genus += 1

        # "reached species" excludes placeholder
        if species != "unclassified species":
            reached_species += 1

        lineages[(family, genus, species)] += 1

    return {
        "lineages": lineages,
        "unassigned": unassigned,
        "unclassified_below_family": unclassified_below_family,
        "family": reached_family,
        "genus": reached_genus,
        "species": reached_species,
    }


# ------------------------------------------------------------
# PROCESS LC AND MC
# ------------------------------------------------------------

for dataset in ["LC", "MC"]:

    raw_expected = FILES[dataset]["raw"]

    tax80 = read_taxonomy(FILES[dataset]["80"])
    tax97 = read_taxonomy(FILES[dataset]["97"])

    # Sanity checks
    if len(tax80) != raw_expected:
        raise RuntimeError(
            f"{dataset} 80%: expected {raw_expected} OTUs, "
            f"found {len(tax80)}"
        )

    if len(tax97) != raw_expected:
        raise RuntimeError(
            f"{dataset} 97%: expected {raw_expected} OTUs, "
            f"found {len(tax97)}"
        )

    if set(tax80) != set(tax97):
        raise RuntimeError(
            f"{dataset}: Feature IDs differ between 80% and 97% outputs"
        )

    s80 = summarise(tax80)
    s97 = summarise(tax97)

    all_lineages = (
        set(s80["lineages"]) |
        set(s97["lineages"])
    )

    # Sort:
    # Family alphabetically,
    # then real genera before "unclassified genus",
    # then real species before "unclassified species"
    def sort_key(x):

        family, genus, species = x

        genus_unclassified = genus == "unclassified genus"
        species_unclassified = species == "unclassified species"

        return (
            family,
            genus_unclassified,
            genus,
            species_unclassified,
            species,
        )

    all_lineages = sorted(all_lineages, key=sort_key)

    outfile = OUTDIR / f"Greengenes2_{dataset}_Priya_format.tsv"

    with open(outfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        # ----------------------------------------
        # SUMMARY
        # ----------------------------------------

        w.writerow([
            "Section",
            "Family",
            "Genus",
            "Species",
            "OTU Count (>=80% identity)",
            "OTU Count (>=97% identity)",
        ])

        w.writerow([
            "Summary",
            "Raw OTUs",
            "",
            "",
            raw_expected,
            "NA",
        ])

        w.writerow([
            "Summary",
            "Unassigned",
            "",
            "",
            s80["unassigned"],
            s97["unassigned"],
        ])

        w.writerow([
            "Summary",
            "Unclassified below family",
            "",
            "",
            s80["unclassified_below_family"],
            s97["unclassified_below_family"],
        ])

        w.writerow([
            "Summary",
            "Classified to family level",
            "",
            "",
            s80["family"],
            s97["family"],
        ])

        w.writerow([
            "Summary",
            "Classified to genus level",
            "",
            "",
            s80["genus"],
            s97["genus"],
        ])

        w.writerow([
            "Summary",
            "Classified to species level",
            "",
            "",
            s80["species"],
            s97["species"],
        ])

        # ----------------------------------------
        # HIERARCHICAL TAXONOMY
        # ----------------------------------------

        for family, genus, species in all_lineages:

            w.writerow([
                "Taxonomy",
                family,
                genus,
                species,
                s80["lineages"].get(
                    (family, genus, species), 0
                ),
                s97["lineages"].get(
                    (family, genus, species), 0
                ),
            ])

    # --------------------------------------------
    # PRINT VALIDATION
    # --------------------------------------------

    print()
    print("=" * 70)
    print(dataset)
    print("=" * 70)

    print("Raw OTUs:", raw_expected)

    print()
    print(">=80%:")
    print("  Unassigned:", s80["unassigned"])
    print(
        "  Unclassified below family:",
        s80["unclassified_below_family"]
    )
    print(
        "  Classified to family:",
        s80["family"]
    )
    print(
        "  Classified to genus:",
        s80["genus"]
    )
    print(
        "  Classified to species:",
        s80["species"]
    )

    print()
    print(">=97%:")
    print("  Unassigned:", s97["unassigned"])
    print(
        "  Unclassified below family:",
        s97["unclassified_below_family"]
    )
    print(
        "  Classified to family:",
        s97["family"]
    )
    print(
        "  Classified to genus:",
        s97["genus"]
    )
    print(
        "  Classified to species:",
        s97["species"]
    )

    # Family accounting should equal raw:
    check80 = (
        s80["unassigned"]
        + s80["unclassified_below_family"]
        + s80["family"]
    )

    check97 = (
        s97["unassigned"]
        + s97["unclassified_below_family"]
        + s97["family"]
    )

    print()
    print(
        "80% accounting:",
        check80,
        "/",
        raw_expected,
        "PASS" if check80 == raw_expected else "FAIL"
    )

    print(
        "97% accounting:",
        check97,
        "/",
        raw_expected,
        "PASS" if check97 == raw_expected else "FAIL"
    )

    print()
    print("Output:", outfile)


print()
print("DONE")
print("Output directory:")
print(OUTDIR)

