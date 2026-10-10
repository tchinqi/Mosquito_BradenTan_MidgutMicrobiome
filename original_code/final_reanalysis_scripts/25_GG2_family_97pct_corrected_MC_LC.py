#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")
GG = RUN / "06_taxonomy/greengenes2_MC_LC_final"

# IMPORTANT: corrected Priya naming
#
# MC = max(0.5%, 20 reads) = old directory labelled LC
# LC = 0.1% = old directory labelled MC

SETS = {
    "MC": {
        "raw": 867,
        "tax": GG / "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv",
    },
    "LC": {
        "raw": 7642,
        "tax": GG / "MC_0p1pct/exported_taxonomy/taxonomy.tsv",
    },
}

OUT = GG / "GG2_family_97pct_corrected_MC_vs_LC.tsv"


def get_family(tax):
    for part in tax.split(";"):
        part = part.strip()
        if part.startswith("f__"):
            value = part[3:].strip()
            if value:
                return value
    return ""


results = {}

for label, info in SETS.items():

    counts = Counter()
    total = 0

    with open(info["tax"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            total += 1

            tax = (row["Taxon"] or "").strip()

            if not tax or tax.lower() == "unassigned":
                counts["Unassigned"] += 1
                continue

            family = get_family(tax)

            if family:
                counts[family] += 1
            else:
                counts["Unclassified"] += 1

    if total != info["raw"]:
        raise RuntimeError(
            "{} expected {} raw OTUs, found {}".format(
                label, info["raw"], total
            )
        )

    if sum(counts.values()) != total:
        raise RuntimeError(
            "{} accounting failed: {} != {}".format(
                label, sum(counts.values()), total
            )
        )

    results[label] = counts


families = (
    set(results["MC"].keys()) |
    set(results["LC"].keys())
)

families.discard("Unassigned")
families.discard("Unclassified")

families = sorted(
    families,
    key=lambda f: (
        -results["LC"].get(f, 0),
        -results["MC"].get(f, 0),
        f
    )
)


with open(OUT, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")

    writer.writerow([
        "Family_level_ID",
        "MC_97pct_OTUs",
        "LC_97pct_OTUs"
    ])

    writer.writerow([
        "Raw OTU",
        "NA",
        "NA"
    ])

    writer.writerow([
        "Unassigned",
        results["MC"]["Unassigned"],
        results["LC"]["Unassigned"]
    ])

    writer.writerow([
        "Unclassified",
        results["MC"]["Unclassified"],
        results["LC"]["Unclassified"]
    ])

    for family in families:
        writer.writerow([
            family,
            results["MC"].get(family, 0),
            results["LC"].get(family, 0)
        ])


print("===== GG2 FAMILY-LEVEL >=97% TABLE =====")
print("Corrected naming:")
print("MC = max(0.5%, 20 reads) = 867 raw OTUs")
print("LC = 0.1% = 7642 raw OTUs")
print()

for label in ["MC", "LC"]:
    c = results[label]

    family_assigned = sum(
        n for name, n in c.items()
        if name not in {"Unassigned", "Unclassified"}
    )

    print(label)
    print("  Raw OTUs:", SETS[label]["raw"])
    print("  Unassigned:", c["Unassigned"])
    print("  Unclassified at family:", c["Unclassified"])
    print("  Assigned to family:", family_assigned)
    print("  Accounting total:", sum(c.values()))
    print()

print("Distinct families across MC/LC:", len(families))
print("Output:", OUT)
