#!/usr/bin/env python3

import csv
from pathlib import Path

RUN = Path(
    "/hpctmp/e0958567/Braden_16S_20260819/"
    "20_FINAL_reanalysis_20260903"
)

GG = (
    RUN / "06_taxonomy/greengenes2_corrected_names/"
    "final_80_vs_97_tables"
)

KR = RUN / "06_taxonomy/kraken2_final"

OUT = RUN / "06_taxonomy/final_combined_taxonomy_tables"
OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# GREENGENES2
# ============================================================

GG_FILES = {
    "Family": GG / "GG2_family_MC_LC_80_vs_97.tsv",
    "Genus": GG / "GG2_genus_MC_LC_80_vs_97.tsv",
    "Species": GG / "GG2_species_MC_LC_80_vs_97.tsv",
}


def read_gg(path):
    with open(path, newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


gg_data = {
    rank: read_gg(path)
    for rank, path in GG_FILES.items()
}


def make_gg_table(dataset):

    # Source columns:
    # 0 taxon
    # 1 MC80
    # 2 MC97
    # 3 LC80
    # 4 LC97

    if dataset == "MC":
        c80, c97 = 1, 2
    else:
        c80, c97 = 3, 4

    outpath = OUT / f"Greengenes2_{dataset}_combined.tsv"

    with open(outpath, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Taxonomic_level",
            "Identification",
            f"{dataset}_>=80pct_identity",
            f"{dataset}_>=97pct_identity"
        ])

        # Use Family table for common input summary.
        fam = gg_data["Family"]

        raw = fam[1]
        unassigned = fam[2]

        w.writerow([
            "Summary", "Raw OTUs",
            raw[c80], raw[c97]
        ])

        w.writerow([
            "Summary", "Unassigned",
            unassigned[c80], unassigned[c97]
        ])

        for rank in ["Family", "Genus", "Species"]:

            rows = gg_data[rank]

            # Unclassified is row 4 in original table
            unclassified = rows[3]

            w.writerow([
                rank,
                f"Unclassified at {rank.lower()}",
                unclassified[c80],
                unclassified[c97]
            ])

            for row in rows[4:]:

                if not row or not row[0].strip():
                    continue

                w.writerow([
                    rank,
                    row[0],
                    row[c80],
                    row[c97]
                ])

    print("Created:", outpath)


make_gg_table("LC")
make_gg_table("MC")


# ============================================================
# KRAKEN2
# ============================================================

KR_FILES = {
    "LC": {
        "full_report": KR / "LC_full.report.tsv",
        "full_output": KR / "LC_full.kraken.tsv",
        "rescue_report": KR / "LC_GG2_97_unassigned.report.tsv",
        "rescue_output": KR / "LC_GG2_97_unassigned.kraken.tsv",
    },

    "MC": {
        "full_report": KR / "MC_full.report.tsv",
        "full_output": KR / "MC_full.kraken.tsv",
        "rescue_report": KR / "MC_GG2_97_unassigned.report.tsv",
        "rescue_output": KR / "MC_GG2_97_unassigned.kraken.tsv",
    },
}


def read_report(path):

    records = []

    with open(path) as fh:

        for line in fh:

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 6:
                continue

            records.append({
                "percent": float(parts[0]),
                "clade": int(parts[1]),
                "direct": int(parts[2]),
                "rank": parts[3].strip(),
                "taxid": parts[4].strip(),
                "name": parts[5].strip(),
            })

    return records


def output_counts(path):

    total = classified = unclassified = 0

    with open(path) as fh:

        for line in fh:

            total += 1

            if line.startswith("C\t"):
                classified += 1

            elif line.startswith("U\t"):
                unclassified += 1

    return total, classified, unclassified


def taxid_clade(records, taxid):

    for r in records:
        if r["taxid"] == str(taxid):
            return r["clade"]

    return 0


def rank_map(records, rank):

    return {
        r["taxid"]: (r["name"], r["clade"])
        for r in records
        if r["rank"] == rank
    }


def make_kraken_table(dataset):

    files = KR_FILES[dataset]

    full = read_report(files["full_report"])
    rescue = read_report(files["rescue_report"])

    ftotal, fclass, funclass = output_counts(
        files["full_output"]
    )

    rtotal, rclass, runclass = output_counts(
        files["rescue_output"]
    )

    # Broad domains by NCBI taxid
    broad = [
        ("Bacteria", 2),
        ("Archaea", 2157),
        ("Eukaryota", 2759),
        ("Viruses", 10239),
    ]

    full_broad = {
        name: taxid_clade(full, tid)
        for name, tid in broad
    }

    rescue_broad = {
        name: taxid_clade(rescue, tid)
        for name, tid in broad
    }

    full_other = (
        fclass
        - sum(full_broad.values())
    )

    rescue_other = (
        rclass
        - sum(rescue_broad.values())
    )

    outpath = OUT / f"Kraken2_{dataset}_combined.tsv"

    with open(outpath, "w", newline="") as fh:

        w = csv.writer(fh, delimiter="\t")

        w.writerow([
            "Taxonomic_level",
            "Identification",
            f"{dataset}_full_Kraken2",
            f"{dataset}_GG2_97pct_unassigned_Kraken2"
        ])

        # ---------------- SUMMARY ----------------

        w.writerow([
            "Summary", "Input OTUs",
            ftotal, rtotal
        ])

        w.writerow([
            "Summary", "Classified",
            fclass, rclass
        ])

        w.writerow([
            "Summary", "Unclassified",
            funclass, runclass
        ])

        for name, _ in broad:

            w.writerow([
                "Summary",
                name,
                full_broad[name],
                rescue_broad[name]
            ])

        w.writerow([
            "Summary",
            "Other classified",
            full_other,
            rescue_other
        ])

        # ---------------- TAXONOMIC RANKS ----------------

        for label, code in [
            ("Family", "F"),
            ("Genus", "G"),
            ("Species", "S")
        ]:

            fm = rank_map(full, code)
            rm = rank_map(rescue, code)

            full_rank_total = sum(
                count for _, count in fm.values()
            )

            rescue_rank_total = sum(
                count for _, count in rm.values()
            )

            w.writerow([
                label,
                f"Unclassified at {label.lower()}",
                fclass - full_rank_total,
                rclass - rescue_rank_total
            ])

            taxids = set(fm) | set(rm)

            ordered = sorted(
                taxids,
                key=lambda tid: (
                    -max(
                        fm.get(tid, ("", 0))[1],
                        rm.get(tid, ("", 0))[1]
                    ),
                    fm.get(
                        tid,
                        rm.get(tid, ("", 0))
                    )[0]
                )
            )

            for tid in ordered:

                if tid in fm:
                    name = fm[tid][0]
                else:
                    name = rm[tid][0]

                full_count = fm.get(
                    tid, (name, 0)
                )[1]

                rescue_count = rm.get(
                    tid, (name, 0)
                )[1]

                w.writerow([
                    label,
                    name,
                    full_count,
                    rescue_count
                ])

    print("Created:", outpath)


make_kraken_table("LC")
make_kraken_table("MC")

print()
print("ALL FINAL COMBINED TABLES CREATED")
print("Output directory:")
print(OUT)
