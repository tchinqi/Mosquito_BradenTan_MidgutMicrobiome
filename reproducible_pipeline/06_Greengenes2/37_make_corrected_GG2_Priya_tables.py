#!/usr/bin/env python3
# =============================================================================
# Script 37 - 37_make_corrected_GG2_Priya_tables.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Greengenes2 comparison
# Purpose : Summary tables of the 80% vs 97% pairing for LC and MC.
# Input   : GG2_80_97_OTU_pairing/*.tsv
# Output  : GG2_80_97_Priya_corrected/
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

INDIR = (
    RUN /
    "06_taxonomy/greengenes2_corrected_names/"
    "GG2_80_97_OTU_pairing"
)

OUTDIR = (
    RUN /
    "06_taxonomy/greengenes2_corrected_names/"
    "GG2_80_97_Priya_corrected"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "LC": {
        "expected": 7642,
        "file": INDIR / "Greengenes2_LC_80_vs_97_per_OTU.tsv",
    },
    "MC": {
        "expected": 867,
        "file": INDIR / "Greengenes2_MC_80_vs_97_per_OTU.tsv",
    },
}


def clean(v):
    return (v or "").strip()


def summarise_rank(rows, suffix):

    unassigned = 0
    below_family = 0
    family = 0
    genus = 0
    species = 0

    for r in rows:

        status = clean(r[f"Status_{suffix}"])
        fam = clean(r[f"Family_{suffix}"])
        gen = clean(r[f"Genus_{suffix}"])
        sp = clean(r[f"Species_{suffix}"])

        if status == "Unassigned":
            unassigned += 1

        elif not fam:
            # Assigned somewhere above family
            below_family += 1

        else:
            family += 1

            if gen:
                genus += 1

            if sp:
                species += 1

    return {
        "unassigned": unassigned,
        "below_family": below_family,
        "family": family,
        "genus": genus,
        "species": species,
    }


for dataset, cfg in DATASETS.items():

    print()
    print("=" * 72)
    print(dataset)
    print("=" * 72)

    # --------------------------------------------------------
    # Read paired OTU-level file
    # --------------------------------------------------------

    with open(cfg["file"], newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if len(rows) != cfg["expected"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected']} OTUs, "
            f"found {len(rows)}"
        )

    ids = [clean(r["Feature_ID"]) for r in rows]

    if len(ids) != len(set(ids)):
        raise RuntimeError(
            f"{dataset}: duplicate Feature IDs detected"
        )

    # --------------------------------------------------------
    # Independent summary counts for 80% and 97%
    # --------------------------------------------------------

    s80 = summarise_rank(rows, "80")
    s97 = summarise_rank(rows, "97")

    # --------------------------------------------------------
    # Build taxonomy table anchored to the 80% lineage
    #
    # IMPORTANT:
    # The 97% count uses the SAME OTUs belonging to that
    # 80%-defined lineage.
    #
    # Therefore no row can have 97% > 80%.
    # --------------------------------------------------------

    count80 = Counter()
    count97_same_otus = Counter()

    changed_rows = []

    for r in rows:

        status80 = clean(r["Status_80"])
        status97 = clean(r["Status_97"])

        family80 = clean(r["Family_80"])
        genus80 = clean(r["Genus_80"])
        species80 = clean(r["Species_80"])

        # Preserve taxonomy-change audit
        if clean(r["Transition"]) == "ASSIGNED_TAXONOMY_CHANGED":
            changed_rows.append(r)

        # Priya taxonomy table starts from family.
        # Completely unassigned and assigned-above-family OTUs
        # are kept in the summary rather than taxonomy rows.
        if status80 == "Unassigned" or not family80:
            continue

        genus_display = (
            genus80 if genus80 else "unclassified genus"
        )

        species_display = (
            species80 if species80 else "unclassified species"
        )

        lineage80 = (
            family80,
            genus_display,
            species_display,
        )

        # Baseline count
        count80[lineage80] += 1

        # SAME OTU comparison:
        # count this OTU in the 97% column if it remains
        # assigned in the 97% GG2 run.
        #
        # Its 97% taxonomic label may become more specific
        # or otherwise change; that change is retained in
        # the audit file instead of moving the OTU into
        # another row.
        if status97 == "Assigned":
            count97_same_otus[lineage80] += 1

    # --------------------------------------------------------
    # Validate monotonicity
    # --------------------------------------------------------

    violations = []

    for lineage, n80 in count80.items():

        n97 = count97_same_otus.get(lineage, 0)

        if n97 > n80:
            violations.append(
                (lineage, n80, n97)
            )

    if violations:
        raise RuntimeError(
            f"{dataset}: found {len(violations)} rows "
            f"where 97% > 80%"
        )

    # --------------------------------------------------------
    # Sort taxonomy rows
    # --------------------------------------------------------

    def sort_key(lineage):

        family, genus, species = lineage

        return (
            family.lower(),
            genus == "unclassified genus",
            genus.lower(),
            species == "unclassified species",
            species.lower(),
        )

    lineages = sorted(
        count80.keys(),
        key=sort_key
    )

    # --------------------------------------------------------
    # Write corrected Priya table
    # --------------------------------------------------------

    outfile = (
        OUTDIR /
        f"Greengenes2_{dataset}_Priya_CORRECTED.tsv"
    )

    with open(outfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Section",
            "Family",
            "Genus",
            "Species",
            "OTU Count (>=80% identity)",
            "OTU Count (same OTUs after >=97% cutoff)",
        ])

        # Summary
        w.writerow([
            "Summary",
            "Raw OTUs",
            "",
            "",
            cfg["expected"],
            cfg["expected"],
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
            s80["below_family"],
            s97["below_family"],
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

        # Taxonomy rows
        for family, genus, species in lineages:

            lineage = (family, genus, species)

            w.writerow([
                "Taxonomy",
                family,
                genus,
                species,
                count80[lineage],
                count97_same_otus.get(lineage, 0),
            ])

    # --------------------------------------------------------
    # Write changed-taxonomy audit subset
    # --------------------------------------------------------

    changefile = (
        OUTDIR /
        f"Greengenes2_{dataset}_80_97_taxonomy_changes.tsv"
    )

    if changed_rows:

        with open(changefile, "w", newline="") as fh:

            fieldnames = list(changed_rows[0].keys())

            w = csv.DictWriter(
                fh,
                delimiter="\t",
                fieldnames=fieldnames
            )

            w.writeheader()
            w.writerows(changed_rows)

    # --------------------------------------------------------
    # Validation reporting
    # --------------------------------------------------------

    print("Raw OTUs:", cfg["expected"])

    print()
    print("80% summary:")
    print("  Unassigned:", s80["unassigned"])
    print("  Unclassified below family:", s80["below_family"])
    print("  Family:", s80["family"])
    print("  Genus:", s80["genus"])
    print("  Species:", s80["species"])

    print()
    print("97% summary:")
    print("  Unassigned:", s97["unassigned"])
    print("  Unclassified below family:", s97["below_family"])
    print("  Family:", s97["family"])
    print("  Genus:", s97["genus"])
    print("  Species:", s97["species"])

    print()
    print(
        "80% taxonomy-table OTUs:",
        sum(count80.values())
    )

    print(
        "Same baseline OTUs still assigned at 97%:",
        sum(count97_same_otus.values())
    )

    print(
        "Rows with 97% > 80%:",
        len(violations),
        "PASS" if len(violations) == 0 else "FAIL"
    )

    print(
        "Taxonomy-changed OTUs retained in audit:",
        len(changed_rows)
    )

    # --------------------------------------------------------
    # Check Acidomonas methanolica issue
    # --------------------------------------------------------

    print()
    print("Rows involving Acidomonas / Acetobacteraceae:")

    shown = 0

    for lineage in lineages:

        family, genus, species = lineage

        if (
            "Acidomonas" in genus
            or "Acidomonas" in species
            or (
                family == "Acetobacteraceae"
                and genus == "unclassified genus"
            )
        ):

            print(
                lineage,
                "80% =", count80[lineage],
                "97% =", count97_same_otus.get(lineage, 0)
            )

            shown += 1

            if shown >= 10:
                break

    print()
    print("Corrected table:")
    print(outfile)

    print("Taxonomy-change audit:")
    print(changefile)


print()
print("=" * 72)
print("CORRECTED GG2 PRIYA TABLES COMPLETE")
print("=" * 72)

print()
print("Output directory:")
print(OUTDIR)

