#!/usr/bin/env python3
# =============================================================================
# Script 43 - 43_summarise_FINAL_GG80_Kraken_rescue.py
# LSM4288 FYP: mosquito midgut 16S rRNA (Nanopore, V3-V4)
# Author  : Braden Tan (Tan Chin Qi Braden), NUS; Puniamoorthy lab
# Step    : Kraken2 cross-check
# Purpose : Summarises Kraken2 results by deepest rank and taxon; Fig. 1B,
#           Tables S2C-D.
# Input   : Kraken2 outputs; ktaxonomy.tsv
# Output  : FINAL_rescue_summary/
# Run on  : Atlas9 (Python 3)
# Release : Header added for the code release (10 Oct 2026); comments only,
#           code unchanged. Unedited original: original_code/.
# =============================================================================

import csv
import re
from pathlib import Path
from collections import Counter

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

KROOT = RUN / "06_taxonomy/kraken2_final"

RESCUE_ROOT = (
    KROOT /
    "FINAL_GG80_unassigned_domain_rescue"
)

KRAKEN_DIR = (
    RESCUE_ROOT /
    "kraken_default_vs_conf097"
)

OUTDIR = (
    RESCUE_ROOT /
    "FINAL_rescue_summary"
)

OUTDIR.mkdir(parents=True, exist_ok=True)

DB = Path(
    "/hpctmp/e0958567/kraken2_db/"
    "k2_standard_08gb_20230605"
)

KTAX = DB / "ktaxonomy.tsv"

DATASETS = {
    "LC": {
        "expected": 4884,
        "expected_default_classified": 86,
        "expected_conf097_classified": 0,

        "ids": (
            RESCUE_ROOT /
            "LC_GG80_rescue_OTU_ids.tsv"
        ),

        "default": (
            KRAKEN_DIR /
            "LC_GG80_rescue_default.kraken.tsv"
        ),

        "conf097": (
            KRAKEN_DIR /
            "LC_GG80_rescue_conf097.kraken.tsv"
        ),
    },

    "MC": {
        "expected": 188,
        "expected_default_classified": 16,
        "expected_conf097_classified": 0,

        "ids": (
            RESCUE_ROOT /
            "MC_GG80_rescue_OTU_ids.tsv"
        ),

        "default": (
            KRAKEN_DIR /
            "MC_GG80_rescue_default.kraken.tsv"
        ),

        "conf097": (
            KRAKEN_DIR /
            "MC_GG80_rescue_conf097.kraken.tsv"
        ),
    },
}


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


def read_rescue_ids(path):

    result = {}

    with open(path, newline="") as fh:

        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:

            result[row["Feature_ID"].strip()] = (
                row["GG2_80_category"].strip()
            )

    return result


def read_kraken(path):

    results = {}

    with open(path) as fh:

        for line in fh:

            line = line.rstrip("\n")

            if not line:
                continue

            parts = line.split("\t")

            if len(parts) < 3:
                raise RuntimeError(
                    f"Unexpected Kraken row: {line}"
                )

            status = parts[0].strip()
            feature_id = parts[1].strip()
            assignment = parts[2].strip()

            if status not in {"C", "U"}:
                raise RuntimeError(
                    f"Unexpected Kraken status: {status}"
                )

            taxid = None
            family = None
            genus = None
            species = None

            if status == "C":

                taxid = get_taxid(assignment)

                if taxid is None:
                    raise RuntimeError(
                        f"Could not recover taxid from "
                        f"classified assignment: {assignment}"
                    )

                family, genus, species = get_lineage(taxid)

            results[feature_id] = {
                "status": status,
                "assignment": assignment,
                "taxid": taxid,
                "family": family,
                "genus": genus,
                "species": species,
            }

    return results


def summarise(records):

    total = len(records)

    classified = sum(
        x["status"] == "C"
        for x in records.values()
    )

    unclassified = sum(
        x["status"] == "U"
        for x in records.values()
    )

    family = sum(
        x["status"] == "C" and x["family"] is not None
        for x in records.values()
    )

    genus = sum(
        x["status"] == "C" and x["genus"] is not None
        for x in records.values()
    )

    species = sum(
        x["status"] == "C" and x["species"] is not None
        for x in records.values()
    )

    above_family = classified - family

    return {
        "Input rescue OTUs": total,
        "Kraken classified": classified,
        "Kraken unclassified": unclassified,
        "Kraken classified but not reaching family": above_family,
        "Classified to family level": family,
        "Classified to genus level": genus,
        "Classified to species level": species,
    }


def make_lineage_counter(records):

    counter = Counter()

    for r in records.values():

        if r["status"] != "C":
            continue

        if r["family"] is None:
            continue

        family = r["family"]

        genus = (
            r["genus"]
            if r["genus"] is not None
            else "unclassified genus"
        )

        species = (
            r["species"]
            if r["species"] is not None
            else "unclassified species"
        )

        counter[(family, genus, species)] += 1

    return counter


# ============================================================
# PROCESS LC + MC
# ============================================================

for dataset, cfg in DATASETS.items():

    print()
    print("=" * 74)
    print(dataset)
    print("=" * 74)

    rescue_ids = read_rescue_ids(cfg["ids"])

    default = read_kraken(cfg["default"])
    conf097 = read_kraken(cfg["conf097"])

    expected_ids = set(rescue_ids)
    default_ids = set(default)
    conf097_ids = set(conf097)

    print("Expected rescue OTUs:", len(expected_ids))
    print("Default Kraken OTUs:", len(default_ids))
    print("Conf0.97 Kraken OTUs:", len(conf097_ids))

    if len(expected_ids) != cfg["expected"]:
        raise RuntimeError(
            f"{dataset}: rescue ID count mismatch"
        )

    if default_ids != expected_ids:
        raise RuntimeError(
            f"{dataset}: default Kraken IDs do not "
            f"exactly match rescue FASTA IDs"
        )

    if conf097_ids != expected_ids:
        raise RuntimeError(
            f"{dataset}: conf0.97 Kraken IDs do not "
            f"exactly match rescue FASTA IDs"
        )

    print("Exact Feature-ID match across both Kraken runs: PASS")

    # --------------------------------------------------------
    # Summaries
    # --------------------------------------------------------

    s_default = summarise(default)
    s_097 = summarise(conf097)

    print()
    print("KRAKEN DEFAULT CONFIDENCE 0.0")

    for k, v in s_default.items():
        print(f"  {k}: {v}")

    print()
    print("KRAKEN CONFIDENCE 0.97")

    for k, v in s_097.items():
        print(f"  {k}: {v}")

    # --------------------------------------------------------
    # Expected result validation
    # --------------------------------------------------------

    if (
        s_default["Kraken classified"]
        != cfg["expected_default_classified"]
    ):
        raise RuntimeError(
            f"{dataset}: expected "
            f"{cfg['expected_default_classified']} "
            f"default classified OTUs, found "
            f"{s_default['Kraken classified']}"
        )

    if (
        s_097["Kraken classified"]
        != cfg["expected_conf097_classified"]
    ):
        raise RuntimeError(
            f"{dataset}: expected "
            f"{cfg['expected_conf097_classified']} "
            f"conf0.97 classified OTUs, found "
            f"{s_097['Kraken classified']}"
        )

    if (
        s_default["Kraken classified"]
        + s_default["Kraken unclassified"]
        != cfg["expected"]
    ):
        raise RuntimeError(
            f"{dataset}: default accounting FAIL"
        )

    if (
        s_097["Kraken classified"]
        + s_097["Kraken unclassified"]
        != cfg["expected"]
    ):
        raise RuntimeError(
            f"{dataset}: conf0.97 accounting FAIL"
        )

    print()
    print("VALIDATION")
    print(
        "Default rescue accounting:",
        s_default["Kraken classified"]
        + s_default["Kraken unclassified"],
        "/",
        cfg["expected"],
        "PASS"
    )

    print(
        "Conf0.97 rescue accounting:",
        s_097["Kraken classified"]
        + s_097["Kraken unclassified"],
        "/",
        cfg["expected"],
        "PASS"
    )

    # --------------------------------------------------------
    # Detailed taxonomy counters
    # --------------------------------------------------------

    c_default = make_lineage_counter(default)
    c_097 = make_lineage_counter(conf097)

    lineages = sorted(
        set(c_default) | set(c_097),
        key=lambda x: (
            x[0].lower(),
            x[1].lower(),
            x[2].lower(),
        )
    )

    # --------------------------------------------------------
    # Priya-format TSV
    # --------------------------------------------------------

    priya_file = (
        OUTDIR /
        f"Kraken2_{dataset}_FINAL_GG80_rescue_Priya.tsv"
    )

    with open(priya_file, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Section",
            "Family",
            "Genus",
            "Species",
            "Kraken default confidence 0.0 OTU Count",
            "Kraken confidence 0.97 OTU Count",
        ])

        for metric in [
            "Input rescue OTUs",
            "Kraken classified",
            "Kraken unclassified",
            "Kraken classified but not reaching family",
            "Classified to family level",
            "Classified to genus level",
            "Classified to species level",
        ]:

            w.writerow([
                "Summary",
                metric,
                "",
                "",
                s_default[metric],
                s_097[metric],
            ])

        for family, genus, species in lineages:

            w.writerow([
                "Taxonomy",
                family,
                genus,
                species,
                c_default[(family, genus, species)],
                c_097[(family, genus, species)],
            ])

    # --------------------------------------------------------
    # Per-OTU audit TSV
    # --------------------------------------------------------

    audit_file = (
        OUTDIR /
        f"Kraken2_{dataset}_FINAL_GG80_rescue_per_OTU.tsv"
    )

    with open(audit_file, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Feature_ID",
            "GG2_80_rescue_category",

            "Kraken_default_status",
            "Kraken_default_assignment",
            "Kraken_default_taxid",
            "Kraken_default_family",
            "Kraken_default_genus",
            "Kraken_default_species",

            "Kraken_conf097_status",
            "Kraken_conf097_assignment",
            "Kraken_conf097_taxid",
            "Kraken_conf097_family",
            "Kraken_conf097_genus",
            "Kraken_conf097_species",
        ])

        for feature_id in sorted(expected_ids):

            d = default[feature_id]
            c = conf097[feature_id]

            w.writerow([
                feature_id,
                rescue_ids[feature_id],

                d["status"],
                d["assignment"],
                d["taxid"] or "",
                d["family"] or "",
                d["genus"] or "",
                d["species"] or "",

                c["status"],
                c["assignment"],
                c["taxid"] or "",
                c["family"] or "",
                c["genus"] or "",
                c["species"] or "",
            ])

    print()
    print("Priya table:")
    print(priya_file)

    print()
    print("Per-OTU audit:")
    print(audit_file)


print()
print("=" * 74)
print("FINAL GG80 -> KRAKEN RESCUE SUMMARY COMPLETE")
print("=" * 74)

print()
print("Output directory:")
print(OUTDIR)

