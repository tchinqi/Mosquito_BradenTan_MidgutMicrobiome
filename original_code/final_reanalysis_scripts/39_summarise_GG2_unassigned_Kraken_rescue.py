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

KROOT = RUN / "06_taxonomy/kraken2_final"

OUTDIR = (
    KROOT /
    "GG2_unassigned_rescue_summary"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "LC": {
        "expected80": 4870,
        "expected97": 7134,

        "kraken80": (
            KROOT /
            "GG2_default80_unassigned/"
            "LC_GG2_default80_unassigned.kraken.tsv"
        ),

        "kraken97": (
            KROOT /
            "LC_GG2_97_unassigned.kraken.tsv"
        ),
    },

    "MC": {
        "expected80": 186,
        "expected97": 539,

        "kraken80": (
            KROOT /
            "GG2_default80_unassigned/"
            "MC_GG2_default80_unassigned.kraken.tsv"
        ),

        "kraken97": (
            KROOT /
            "MC_GG2_97_unassigned.kraken.tsv"
        ),
    },
}

KTAX = DB / "ktaxonomy.tsv"

TAXID_RE = re.compile(r"\(taxid\s+(\d+)\)")


# ============================================================
# READ KRAKEN TAXONOMY TREE
# ============================================================

parent = {}
rank = {}
name = {}

with open(KTAX) as fh:

    for line in fh:

        parts = [
            x.strip()
            for x in line.rstrip("\n").split("|")
        ]

        if len(parts) < 5:
            continue

        taxid = parts[0]
        parent_taxid = parts[1]
        tax_rank = parts[2]
        tax_name = parts[4]

        parent[taxid] = parent_taxid
        rank[taxid] = tax_rank
        name[taxid] = tax_name

print("Kraken taxonomy nodes loaded:", len(parent))


# ============================================================
# HELPERS
# ============================================================

def get_taxid(assignment):

    m = TAXID_RE.search(assignment)

    if not m:
        return None

    return m.group(1)


def get_lineage(taxid):

    """
    Walk assigned Kraken taxid towards root and recover
    canonical Family, Genus and Species ancestors.
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

        if r == "S" and species is None:
            species = n

        elif r == "G" and genus is None:
            genus = n

        elif r == "F" and family is None:
            family = n

        if current == "1":
            break

        p = parent.get(current)

        if not p or p == current:
            break

        current = p

    return family, genus, species


def read_kraken(path):

    records = {}

    with open(path) as fh:

        for line in fh:

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            status = parts[0]
            feature_id = parts[1]
            assignment = parts[2]

            if feature_id in records:
                raise RuntimeError(
                    f"Duplicate Kraken Feature ID: {feature_id}"
                )

            taxid = None

            if status == "C":
                taxid = get_taxid(assignment)

                if taxid is None:
                    raise RuntimeError(
                        f"Could not parse taxid for "
                        f"{feature_id}: {assignment}"
                    )

            records[feature_id] = {
                "status": status,
                "assignment": assignment,
                "taxid": taxid,
            }

    return records


def summarise(records):

    lineages = Counter()

    total = len(records)
    classified = 0
    unclassified = 0

    above_family = 0
    reached_family = 0
    reached_genus = 0
    reached_species = 0

    otu_details = {}

    for feature_id, rec in records.items():

        if rec["status"] == "U":

            unclassified += 1

            otu_details[feature_id] = {
                "status": "U",
                "taxid": "0",
                "assignment": "unclassified",
                "family": "",
                "genus": "",
                "species": "",
            }

            continue

        classified += 1

        family, genus, species = get_lineage(
            rec["taxid"]
        )

        otu_details[feature_id] = {
            "status": "C",
            "taxid": rec["taxid"],
            "assignment": rec["assignment"],
            "family": family or "",
            "genus": genus or "",
            "species": species or "",
        }

        # Kraken classified it, but no canonical family
        # ancestor was reached.
        if family is None:

            above_family += 1
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

    if classified + unclassified != total:
        raise RuntimeError(
            "Classified/unclassified accounting failed"
        )

    if (
        above_family + reached_family
        != classified
    ):
        raise RuntimeError(
            "Family-level accounting failed"
        )

    if not (
        reached_species
        <= reached_genus
        <= reached_family
        <= classified
    ):
        raise RuntimeError(
            "Taxonomic-rank hierarchy failed"
        )

    return {
        "total": total,
        "classified": classified,
        "unclassified": unclassified,
        "above_family": above_family,
        "family": reached_family,
        "genus": reached_genus,
        "species": reached_species,
        "lineages": lineages,
        "details": otu_details,
    }


# ============================================================
# PROCESS LC + MC
# ============================================================

for dataset, cfg in DATASETS.items():

    print()
    print("=" * 74)
    print(dataset)
    print("=" * 74)

    rec80 = read_kraken(cfg["kraken80"])
    rec97 = read_kraken(cfg["kraken97"])

    # --------------------------------------------------------
    # Input validation
    # --------------------------------------------------------

    print(
        "GG >=80%-unassigned rescue OTUs:",
        len(rec80)
    )

    print(
        "GG >=97%-unassigned rescue OTUs:",
        len(rec97)
    )

    if len(rec80) != cfg["expected80"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected80']} "
            f"GG80-unassigned rescue OTUs, "
            f"found {len(rec80)}"
        )

    if len(rec97) != cfg["expected97"]:
        raise RuntimeError(
            f"{dataset}: expected {cfg['expected97']} "
            f"GG97-unassigned rescue OTUs, "
            f"found {len(rec97)}"
        )

    ids80 = set(rec80)
    ids97 = set(rec97)

    # Since no OTUs went:
    # Unassigned80 -> Assigned97,
    # the GG80-unassigned set should be contained
    # within the GG97-unassigned set.
    if not ids80.issubset(ids97):

        missing = ids80 - ids97

        raise RuntimeError(
            f"{dataset}: GG80-unassigned IDs are not "
            f"a subset of GG97-unassigned IDs. "
            f"Missing from 97 set: {len(missing)}"
        )

    print(
        "GG80-unassigned IDs subset of "
        "GG97-unassigned IDs: PASS"
    )

    print(
        "Additional OTUs becoming unassigned at GG97:",
        len(ids97 - ids80)
    )

    # --------------------------------------------------------
    # Shared-OTU Kraken consistency
    # --------------------------------------------------------

    inconsistent = []

    for feature_id in sorted(ids80):

        a = rec80[feature_id]
        b = rec97[feature_id]

        if (
            a["status"] != b["status"]
            or a["taxid"] != b["taxid"]
        ):
            inconsistent.append(feature_id)

    print(
        "Shared OTUs with different Kraken result:",
        len(inconsistent)
    )

    if inconsistent:

        print(
            "First inconsistent IDs:",
            inconsistent[:20]
        )

        raise RuntimeError(
            f"{dataset}: Kraken classifications differ "
            f"for shared OTUs"
        )

    print(
        "Shared-OTU Kraken classification consistency: PASS"
    )

    # --------------------------------------------------------
    # Summaries
    # --------------------------------------------------------

    s80 = summarise(rec80)
    s97 = summarise(rec97)

    print()
    print("GG >=80%-UNASSIGNED -> KRAKEN")
    print("  Input OTUs:", s80["total"])
    print("  Kraken classified:", s80["classified"])
    print("  Kraken unclassified:", s80["unclassified"])
    print(
        "  Classified but not reaching family:",
        s80["above_family"]
    )
    print(
        "  Classified to family level:",
        s80["family"]
    )
    print(
        "  Classified to genus level:",
        s80["genus"]
    )
    print(
        "  Classified to species level:",
        s80["species"]
    )

    print()
    print("GG >=97%-UNASSIGNED -> KRAKEN")
    print("  Input OTUs:", s97["total"])
    print("  Kraken classified:", s97["classified"])
    print("  Kraken unclassified:", s97["unclassified"])
    print(
        "  Classified but not reaching family:",
        s97["above_family"]
    )
    print(
        "  Classified to family level:",
        s97["family"]
    )
    print(
        "  Classified to genus level:",
        s97["genus"]
    )
    print(
        "  Classified to species level:",
        s97["species"]
    )

    # --------------------------------------------------------
    # Write Priya-style combined rescue table
    # --------------------------------------------------------

    all_lineages = (
        set(s80["lineages"])
        | set(s97["lineages"])
    )

    def sort_key(x):

        family, genus, species = x

        return (
            family.lower(),
            genus == "unclassified genus",
            genus.lower(),
            species == "unclassified species",
            species.lower(),
        )

    outfile = (
        OUTDIR /
        f"Kraken2_{dataset}_GG2_unassigned_rescue_Priya.tsv"
    )

    with open(outfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Section",
            "Family",
            "Genus",
            "Species",
            "GG2 >=80% Unassigned -> Kraken OTU Count",
            "GG2 >=97% Unassigned -> Kraken OTU Count",
        ])

        summary_rows = [
            (
                "Input GG2-unassigned OTUs",
                s80["total"],
                s97["total"],
            ),
            (
                "Kraken classified",
                s80["classified"],
                s97["classified"],
            ),
            (
                "Kraken unclassified",
                s80["unclassified"],
                s97["unclassified"],
            ),
            (
                "Kraken classified but not reaching family",
                s80["above_family"],
                s97["above_family"],
            ),
            (
                "Classified to family level",
                s80["family"],
                s97["family"],
            ),
            (
                "Classified to genus level",
                s80["genus"],
                s97["genus"],
            ),
            (
                "Classified to species level",
                s80["species"],
                s97["species"],
            ),
        ]

        for label, n80, n97 in summary_rows:

            w.writerow([
                "Summary",
                label,
                "",
                "",
                n80,
                n97,
            ])

        for lineage in sorted(
            all_lineages,
            key=sort_key
        ):

            family, genus, species = lineage

            w.writerow([
                "Taxonomy",
                family,
                genus,
                species,
                s80["lineages"].get(lineage, 0),
                s97["lineages"].get(lineage, 0),
            ])

    # --------------------------------------------------------
    # Per-OTU audit file
    # --------------------------------------------------------

    auditfile = (
        OUTDIR /
        f"Kraken2_{dataset}_GG2_unassigned_rescue_per_OTU.tsv"
    )

    with open(auditfile, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Feature_ID",
            "In_GG80_unassigned_set",
            "In_GG97_unassigned_set",
            "Kraken_status",
            "Kraken_taxid",
            "Kraken_assignment",
            "Family",
            "Genus",
            "Species",
        ])

        for feature_id in sorted(ids97):

            # If shared, rec80 and rec97 are already
            # validated as identical.
            rec = (
                s80["details"][feature_id]
                if feature_id in ids80
                else s97["details"][feature_id]
            )

            w.writerow([
                feature_id,
                "YES" if feature_id in ids80 else "NO",
                "YES",
                rec["status"],
                rec["taxid"],
                rec["assignment"],
                rec["family"],
                rec["genus"],
                rec["species"],
            ])

    # --------------------------------------------------------
    # Final validation
    # --------------------------------------------------------

    print()
    print("VALIDATION")

    print(
        "80 rescue accounting:",
        s80["classified"] + s80["unclassified"],
        "/",
        s80["total"],
        "PASS"
    )

    print(
        "97 rescue accounting:",
        s97["classified"] + s97["unclassified"],
        "/",
        s97["total"],
        "PASS"
    )

    print()
    print("Priya rescue table:")
    print(outfile)

    print("Per-OTU audit:")
    print(auditfile)


print()
print("=" * 74)
print("GG2-UNASSIGNED -> KRAKEN RESCUE SUMMARY COMPLETE")
print("=" * 74)

print()
print("Output directory:")
print(OUTDIR)

