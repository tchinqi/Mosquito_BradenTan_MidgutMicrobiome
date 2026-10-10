#!/usr/bin/env python3

import csv
from pathlib import Path

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")
GG = RUN / "06_taxonomy/greengenes2_MC_LC_final"

MC_RAW = 7642
LC_RAW = 867

MC_GENUS = GG / "MC_0p1pct/all_genera_counts.tsv"
LC_GENUS = GG / "LC_max0p5pct_20reads/all_genera_counts.tsv"

MC_SPECIES = GG / "MC_0p1pct/all_species_counts.tsv"
LC_SPECIES = GG / "LC_max0p5pct_20reads/all_species_counts.tsv"

OUT_GENUS = GG / "GG2_genus_MC_vs_LC.tsv"
OUT_SPECIES = GG / "GG2_species_MC_vs_LC.tsv"


def read_counts(path, taxon_col):
    counts = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            counts[row[taxon_col]] = int(row["OTU_count"])

    return counts


def make_table(mc_file, lc_file, taxon_col, out_file):

    mc = read_counts(mc_file, taxon_col)
    lc = read_counts(lc_file, taxon_col)

    taxa = set(mc) | set(lc)

    # Sort primarily by MC count, then LC count, then taxon name
    taxa = sorted(
        taxa,
        key=lambda x: (-mc.get(x, 0), -lc.get(x, 0), x)
    )

    with open(out_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        writer.writerow([
            taxon_col,
            "MC_OTU_count",
            "MC_percent_of_raw_OTUs",
            "LC_OTU_count",
            "LC_percent_of_raw_OTUs"
        ])

        for taxon in taxa:
            mc_n = mc.get(taxon, 0)
            lc_n = lc.get(taxon, 0)

            writer.writerow([
                taxon,
                mc_n,
                "{:.2f}".format(100 * mc_n / MC_RAW),
                lc_n,
                "{:.2f}".format(100 * lc_n / LC_RAW)
            ])

    print()
    print("==== {} TABLE ====".format(taxon_col.upper()))
    print("Distinct taxa in MC:", len(mc))
    print("Distinct taxa in LC:", len(lc))
    print("Distinct taxa in union:", len(taxa))
    print("Output:", out_file)


make_table(
    MC_GENUS,
    LC_GENUS,
    "Genus",
    OUT_GENUS
)

make_table(
    MC_SPECIES,
    LC_SPECIES,
    "Species",
    OUT_SPECIES
)
