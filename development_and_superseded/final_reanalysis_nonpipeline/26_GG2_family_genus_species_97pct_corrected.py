#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")
GG = RUN / "06_taxonomy/greengenes2_MC_LC_final"

# ------------------------------------------------------------
# CORRECTED PRIYA NAMING
#
# MC = max(0.5%, 20 reads) = 867 OTUs
#      stored in old LC directory
#
# LC = 0.1% = 7642 OTUs
#      stored in old MC directory
# ------------------------------------------------------------

SETS = {
    "MC": {
        "raw": 867,
        "tax": GG / "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv"
    },
    "LC": {
        "raw": 7642,
        "tax": GG / "MC_0p1pct/exported_taxonomy/taxonomy.tsv"
    }
}

RANKS = {
    "family": "f__",
    "genus": "g__",
    "species": "s__"
}


def extract_rank(taxonomy, prefix):
    for part in taxonomy.split(";"):
        part = part.strip()
        if part.startswith(prefix):
            value = part[len(prefix):].strip()
            if value:
                return value
    return ""


def summarise(path, rank_prefix, expected):
    counts = Counter()
    total = 0

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            total += 1
            tax = (row["Taxon"] or "").strip()

            # No qualifying >=97% Greengenes assignment
            if not tax or tax.lower() == "unassigned":
                counts["Unassigned"] += 1
                continue

            value = extract_rank(tax, rank_prefix)

            # Greengenes assigned the OTU, but not to this rank
            if not value:
                counts["Unclassified"] += 1
            else:
                counts[value] += 1

    if total != expected:
        raise RuntimeError(
            "Expected {} OTUs, found {} in {}".format(
                expected, total, path
            )
        )

    if sum(counts.values()) != total:
        raise RuntimeError(
            "Accounting error: {} != {}".format(
                sum(counts.values()), total
            )
        )

    return counts


for rank, prefix in RANKS.items():

    mc = summarise(
        SETS["MC"]["tax"],
        prefix,
        SETS["MC"]["raw"]
    )

    lc = summarise(
        SETS["LC"]["tax"],
        prefix,
        SETS["LC"]["raw"]
    )

    taxa = (set(mc) | set(lc)) - {
        "Unassigned",
        "Unclassified"
    }

    taxa = sorted(
        taxa,
        key=lambda x: (
            -mc.get(x, 0),
            -lc.get(x, 0),
            x
        )
    )

    out = GG / "GG2_{}_corrected_MC_LC_table.tsv".format(rank)

    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        writer.writerow([
            "{}_level_ID".format(rank.capitalize()),
            "MC",
            "MC_97pct_identity",
            "LC",
            "LC_97pct_identity"
        ])

        # Priya-style raw OTU row
        writer.writerow([
            "Raw OTU",
            SETS["MC"]["raw"],
            "NA",
            SETS["LC"]["raw"],
            "NA"
        ])

        # Broader/default GG2 columns left blank until 80% run finishes
        writer.writerow([
            "Unassigned",
            "",
            mc["Unassigned"],
            "",
            lc["Unassigned"]
        ])

        writer.writerow([
            "Unclassified",
            "",
            mc["Unclassified"],
            "",
            lc["Unclassified"]
        ])

        for taxon in taxa:
            writer.writerow([
                taxon,
                "",
                mc.get(taxon, 0),
                "",
                lc.get(taxon, 0)
            ])

    print()
    print("========== {} ==========".format(rank.upper()))
    print("MC raw OTUs:", SETS["MC"]["raw"])
    print("MC >=97% Unassigned:", mc["Unassigned"])
    print("MC >=97% Unclassified at {}: {}".format(
        rank, mc["Unclassified"]
    ))
    print("MC distinct {} labels: {}".format(
        rank, len([x for x in mc if x not in {
            "Unassigned", "Unclassified"
        }])
    ))

    print()
    print("LC raw OTUs:", SETS["LC"]["raw"])
    print("LC >=97% Unassigned:", lc["Unassigned"])
    print("LC >=97% Unclassified at {}: {}".format(
        rank, lc["Unclassified"]
    ))
    print("LC distinct {} labels: {}".format(
        rank, len([x for x in lc if x not in {
            "Unassigned", "Unclassified"
        }])
    ))

    print()
    print("Output:", out)

