#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")

OLD97 = RUN / "06_taxonomy/greengenes2_MC_LC_final"
NEW80 = RUN / "06_taxonomy/greengenes2_corrected_names"

# FINAL CORRECTED NAMING
# MC = max(0.5%, 20 reads) = 867 OTUs
# LC = 0.1% = 7642 OTUs

FILES = {
    "MC": {
        "raw": 867,
        "80": NEW80 / "MC_max0p5pct_20reads/default80/exported_taxonomy/taxonomy.tsv",
        "97": OLD97 / "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv",
    },
    "LC": {
        "raw": 7642,
        "80": NEW80 / "LC_0p1pct/default80/exported_taxonomy/taxonomy.tsv",
        "97": OLD97 / "MC_0p1pct/exported_taxonomy/taxonomy.tsv",
    },
}

RANKS = {
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}

OUTDIR = RUN / "06_taxonomy/greengenes2_corrected_names/final_80_vs_97_tables"
OUTDIR.mkdir(parents=True, exist_ok=True)


def read_taxonomy(path):
    rows = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            feature = row["Feature ID"].strip()
            taxon = (row["Taxon"] or "").strip()
            rows[feature] = taxon
    return rows


def extract_rank(taxonomy, prefix):
    for part in taxonomy.split(";"):
        part = part.strip()
        if part.startswith(prefix):
            value = part[len(prefix):].strip()
            if value:
                return value
    return ""


def summarise(taxdict, prefix):
    counts = Counter()

    for tax in taxdict.values():

        if not tax or tax.lower() == "unassigned":
            counts["Unassigned"] += 1
            continue

        rank_value = extract_rank(tax, prefix)

        if rank_value:
            counts[rank_value] += 1
        else:
            counts["Unclassified"] += 1

    return counts


# Load and validate all four taxonomy outputs
tax = {}

for dataset in ["MC", "LC"]:
    tax[dataset] = {}

    for threshold in ["80", "97"]:
        path = FILES[dataset][threshold]
        tax[dataset][threshold] = read_taxonomy(path)

        observed = len(tax[dataset][threshold])
        expected = FILES[dataset]["raw"]

        if observed != expected:
            raise RuntimeError(
                f"{dataset} {threshold}%: expected {expected} OTUs, found {observed}"
            )

    ids80 = set(tax[dataset]["80"])
    ids97 = set(tax[dataset]["97"])

    if ids80 != ids97:
        raise RuntimeError(
            f"{dataset}: OTU IDs differ between 80% and 97% classifications"
        )

    print(
        f"{dataset}: PASS — {FILES[dataset]['raw']} identical OTU IDs "
        f"between 80% and 97% runs"
    )


for rank, prefix in RANKS.items():

    MC80 = summarise(tax["MC"]["80"], prefix)
    MC97 = summarise(tax["MC"]["97"], prefix)
    LC80 = summarise(tax["LC"]["80"], prefix)
    LC97 = summarise(tax["LC"]["97"], prefix)

    all_taxa = (
        set(MC80) | set(MC97) |
        set(LC80) | set(LC97)
    )

    all_taxa.discard("Unassigned")
    all_taxa.discard("Unclassified")

    # Largest observed counts first
    all_taxa = sorted(
        all_taxa,
        key=lambda x: (
            -max(
                MC80.get(x, 0),
                MC97.get(x, 0),
                LC80.get(x, 0),
                LC97.get(x, 0)
            ),
            x
        )
    )

    outfile = OUTDIR / f"GG2_{rank}_MC_LC_80_vs_97.tsv"

    with open(outfile, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        writer.writerow([
            f"{rank.capitalize()}-level identification",
            "MC",
            "MC (>=97% identity)",
            "LC",
            "LC (>=97% identity)",
        ])

        writer.writerow([
            "Raw OTU",
            FILES["MC"]["raw"],
            "NA",
            FILES["LC"]["raw"],
            "NA",
        ])

        writer.writerow([
            "Unassigned",
            MC80["Unassigned"],
            MC97["Unassigned"],
            LC80["Unassigned"],
            LC97["Unassigned"],
        ])

        writer.writerow([
            "Unclassified",
            MC80["Unclassified"],
            MC97["Unclassified"],
            LC80["Unclassified"],
            LC97["Unclassified"],
        ])

        for name in all_taxa:
            writer.writerow([
                name,
                MC80.get(name, 0),
                MC97.get(name, 0),
                LC80.get(name, 0),
                LC97.get(name, 0),
            ])

    print()
    print("=" * 65)
    print(rank.upper())
    print("=" * 65)

    for dataset, d80, d97 in [
        ("MC", MC80, MC97),
        ("LC", LC80, LC97),
    ]:
        raw = FILES[dataset]["raw"]

        assigned80 = raw - d80["Unassigned"] - d80["Unclassified"]
        assigned97 = raw - d97["Unassigned"] - d97["Unclassified"]

        print(f"{dataset}:")
        print(
            f"  80%  -> Unassigned {d80['Unassigned']}, "
            f"Unclassified {d80['Unclassified']}, "
            f"classified to {rank} {assigned80}"
        )
        print(
            f"  97%  -> Unassigned {d97['Unassigned']}, "
            f"Unclassified {d97['Unclassified']}, "
            f"classified to {rank} {assigned97}"
        )

    print(f"Output: {outfile}")

print()
print("ALL THREE TABLES COMPLETE")
print("Output directory:", OUTDIR)
