#!/usr/bin/env python3

import csv
import re
from pathlib import Path
from collections import Counter

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

DB = Path(
    "/hpctmp/e0958567/kraken2_db/"
    "k2_standard_08gb_20230605"
)

KROUT = RUN / "06_taxonomy/kraken2_final"

OUTDIR = (
    RUN /
    "06_taxonomy/kraken2_final/"
    "Priya_hierarchical_tables"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

FILES = {
    "LC": {
        "raw": 7642,
        "kraken": KROUT / "LC_full.kraken.tsv",
    },
    "MC": {
        "raw": 867,
        "kraken": KROUT / "MC_full.kraken.tsv",
    },
}

KTAX = DB / "ktaxonomy.tsv"


# ============================================================
# READ KRAKEN TAXONOMY TREE
# ============================================================

parent = {}
rank = {}
name = {}

with open(KTAX) as fh:

    for line in fh:

        # Format:
        # taxid | parent_taxid | rank | depth | name

        parts = [x.strip() for x in line.rstrip("\n").split("|")]

        if len(parts) < 5:
            continue

        taxid = parts[0]
        parent_taxid = parts[1]
        tax_rank = parts[2]
        tax_name = parts[4]

        parent[taxid] = parent_taxid
        rank[taxid] = tax_rank
        name[taxid] = tax_name


print("Taxonomy nodes loaded:", len(parent))


# ============================================================
# HELPERS
# ============================================================

TAXID_RE = re.compile(r"\(taxid\s+(\d+)\)")


def get_taxid(field):

    m = TAXID_RE.search(field)

    if not m:
        return None

    return m.group(1)


def get_lineage(taxid):

    """
    Walk from Kraken-assigned taxid towards root.

    Returns canonical:
    family, genus, species

    Kraken may classify directly at species, genus, family,
    order, etc. Canonical ancestors are obtained from the
    Kraken taxonomy tree.
    """

    family = None
    genus = None
    species = None

    current = str(taxid)
    seen = set()

    while current and current not in seen:

        seen.add(current)

        r = rank.get(current)
        n = name.get(current)

        # First encountered canonical rank is the
        # closest relevant ancestor.
        if r == "S" and species is None:
            species = n

        elif r == "G" and genus is None:
            genus = n

        elif r == "F" and family is None:
            family = n

        # Root
        if current == "1":
            break

        p = parent.get(current)

        if not p or p == current:
            break

        current = p

    return family, genus, species


# ============================================================
# PROCESS EACH DATASET
# ============================================================

for dataset in ["LC", "MC"]:

    infile = FILES[dataset]["kraken"]
    expected = FILES[dataset]["raw"]

    lineages = Counter()

    total = 0

    # Kraken "U" = completely unclassified
    unclassified = 0

    # Kraken classified somewhere in taxonomy,
    # but did not reach canonical family rank
    unclassified_below_family = 0

    reached_family = 0
    reached_genus = 0
    reached_species = 0

    # Audit output at OTU level
    audit_rows = []

    with open(infile) as fh:

        for line in fh:

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            total += 1

            status = parts[0]
            feature_id = parts[1]
            assignment = parts[2]

            if status == "U":

                unclassified += 1

                audit_rows.append([
                    feature_id,
                    status,
                    "0",
                    "unclassified",
                    "",
                    "",
                    "",
                ])

                continue

            taxid = get_taxid(assignment)

            if taxid is None:

                raise RuntimeError(
                    f"Could not extract taxid for "
                    f"{feature_id}: {assignment}"
                )

            family, genus, species = get_lineage(taxid)

            assigned_name = name.get(
                taxid,
                assignment
            )

            # Classified by Kraken but only above family
            if family is None:

                unclassified_below_family += 1

                audit_rows.append([
                    feature_id,
                    status,
                    taxid,
                    assigned_name,
                    "",
                    "",
                    "",
                ])

                continue

            reached_family += 1

            if genus is not None:
                reached_genus += 1
            else:
                genus = "unclassified genus"

            if species is not None:
                reached_species += 1
            else:
                species = "unclassified species"

            lineages[
                (family, genus, species)
            ] += 1

            audit_rows.append([
                feature_id,
                status,
                taxid,
                assigned_name,
                family,
                genus,
                species,
            ])

    # ========================================================
    # VALIDATION
    # ========================================================

    if total != expected:
        raise RuntimeError(
            f"{dataset}: expected {expected} OTUs, "
            f"found {total}"
        )

    accounting = (
        unclassified
        + unclassified_below_family
        + reached_family
    )

    if accounting != expected:
        raise RuntimeError(
            f"{dataset}: accounting failed "
            f"{accounting} != {expected}"
        )

    lineage_total = sum(lineages.values())

    if lineage_total != reached_family:
        raise RuntimeError(
            f"{dataset}: lineage accounting failed "
            f"{lineage_total} != {reached_family}"
        )

    if not (
        reached_species
        <= reached_genus
        <= reached_family
    ):
        raise RuntimeError(
            f"{dataset}: rank hierarchy check failed"
        )

    # ========================================================
    # OUTPUT PRIYA TABLE
    # ========================================================

    outfile = (
        OUTDIR /
        f"Kraken2_{dataset}_Priya_format.tsv"
    )

    with open(outfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Section",
            "Family",
            "Genus",
            "Species",
            "OTU Count",
        ])

        w.writerow([
            "Summary",
            "Raw OTUs",
            "",
            "",
            expected,
        ])

        w.writerow([
            "Summary",
            "Unclassified",
            "",
            "",
            unclassified,
        ])

        w.writerow([
            "Summary",
            "Unclassified below family",
            "",
            "",
            unclassified_below_family,
        ])

        w.writerow([
            "Summary",
            "Classified to family level",
            "",
            "",
            reached_family,
        ])

        w.writerow([
            "Summary",
            "Classified to genus level",
            "",
            "",
            reached_genus,
        ])

        w.writerow([
            "Summary",
            "Classified to species level",
            "",
            "",
            reached_species,
        ])

        # Alphabetical hierarchical ordering:
        # family -> genus -> species
        def sort_key(x):

            family, genus, species = x

            return (
                family.lower(),
                genus == "unclassified genus",
                genus.lower(),
                species == "unclassified species",
                species.lower(),
            )

        for family, genus, species in sorted(
            lineages,
            key=sort_key
        ):

            w.writerow([
                "Taxonomy",
                family,
                genus,
                species,
                lineages[
                    (family, genus, species)
                ],
            ])

    # ========================================================
    # OUTPUT PER-OTU AUDIT FILE
    # ========================================================

    auditfile = (
        OUTDIR /
        f"Kraken2_{dataset}_per_OTU_lineage.tsv"
    )

    with open(auditfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Feature_ID",
            "Kraken_status",
            "Kraken_taxid",
            "Kraken_assignment",
            "Family",
            "Genus",
            "Species",
        ])

        w.writerows(audit_rows)

    # ========================================================
    # REPORT
    # ========================================================

    print()
    print("=" * 70)
    print(dataset)
    print("=" * 70)

    print("Raw OTUs:", expected)
    print("Unclassified:", unclassified)

    print(
        "Unclassified below family:",
        unclassified_below_family
    )

    print(
        "Classified to family:",
        reached_family
    )

    print(
        "Classified to genus:",
        reached_genus
    )

    print(
        "Classified to species:",
        reached_species
    )

    print()
    print(
        "Accounting:",
        accounting,
        "/",
        expected,
        "PASS"
        if accounting == expected
        else "FAIL"
    )

    print(
        "Hierarchical lineage accounting:",
        lineage_total,
        "/",
        reached_family,
        "PASS"
        if lineage_total == reached_family
        else "FAIL"
    )

    print()
    print("Priya table:", outfile)
    print("Audit table:", auditfile)


print()
print("DONE")
print("Output directory:")
print(OUTDIR)

