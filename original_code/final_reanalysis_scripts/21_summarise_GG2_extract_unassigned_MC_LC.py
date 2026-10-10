#!/usr/bin/env python3

import csv
from pathlib import Path

RUN = Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903")
GG = RUN / "06_taxonomy/greengenes2_MC_LC_final"

SETS = {
    "MC": {
        "label": "MC_0p1pct",
        "expected": 7642,
        "tax": GG / "MC_0p1pct/exported_taxonomy/taxonomy.tsv",
        "fasta": RUN / "04_global_vsearch97/postracon_thresholds_576/vsearch97/0p1pct/otu97_centroids_UPPERCASE.fasta",
    },
    "LC": {
        "label": "LC_max0p5pct_20reads",
        "expected": 867,
        "tax": GG / "LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv",
        "fasta": RUN / "04_global_vsearch97/postracon_thresholds_576/LC_max0p5pct_20reads/vsearch97/otu97_centroids_UPPERCASE.fasta",
    }
}

def read_fasta(path):
    records = {}
    current = None
    lines = []

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(lines)
                current = line[1:].strip().split()[0]
                lines = []
            else:
                lines.append(line)

        if current is not None:
            records[current] = "".join(lines)

    return records

for name, info in SETS.items():

    out = GG / info["label"]
    summary_file = out / "Priya_GG2_summary.tsv"
    unassigned_ids_file = out / "unassigned_OTU_ids.txt"
    unassigned_fasta = out / "unassigned_OTUs.fasta"

    with open(info["tax"], newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if len(rows) != info["expected"]:
        raise RuntimeError(
            "{} expected {} OTUs but found {}".format(
                name, info["expected"], len(rows)
            )
        )

    bacteria = []
    archaea = []
    other = []
    unassigned = []
    wolbachia = []

    for row in rows:
        fid = row["Feature ID"]
        tax = (row["Taxon"] or "").strip()
        low = tax.lower()

        if not tax or low == "unassigned":
            unassigned.append(fid)
        elif "d__bacteria" in low:
            bacteria.append(fid)
            if "wolbachia" in low:
                wolbachia.append(fid)
        elif "d__archaea" in low:
            archaea.append(fid)
        else:
            other.append(fid)

    total = len(rows)

    summary = [
        ["Raw OTUs", total, 100.0,
         "All VSEARCH97 OTUs entering Greengenes2"],
        ["Bacteria assigned at >=97% identity", len(bacteria),
         100.0 * len(bacteria) / total,
         "Qualifying Greengenes2 bacterial assignment"],
        ["Archaea assigned at >=97% identity", len(archaea),
         100.0 * len(archaea) / total,
         "Qualifying Greengenes2 archaeal assignment"],
        ["Other assigned at >=97% identity", len(other),
         100.0 * len(other) / total,
         "Assigned but not Bacteria/Archaea"],
        ["Unassigned", len(unassigned),
         100.0 * len(unassigned) / total,
         "No qualifying Greengenes2 assignment"],
        ["Wolbachia (subset of Bacteria)", len(wolbachia),
         100.0 * len(wolbachia) / total,
         "Subset only; do not add to Bacteria total"],
    ]

    with open(summary_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "Classification_category",
            "Number_of_OTUs",
            "Percent_of_raw_OTUs",
            "Notes"
        ])

        for row in summary:
            writer.writerow([
                row[0],
                row[1],
                "{:.2f}".format(row[2]),
                row[3]
            ])

    with open(unassigned_ids_file, "w") as fh:
        for fid in sorted(unassigned):
            fh.write(fid + "\n")

    fasta = read_fasta(info["fasta"])

    missing = [fid for fid in unassigned if fid not in fasta]

    if missing:
        raise RuntimeError(
            "{} unassigned OTU IDs missing from {} FASTA".format(
                len(missing), name
            )
        )

    with open(unassigned_fasta, "w") as fh:
        for fid in sorted(unassigned):
            fh.write(">{}\n".format(fid))
            fh.write(fasta[fid])

    print()
    print("========== {} ==========".format(name))
    print("Raw OTUs:", total)
    print("Bacteria >=97%:", len(bacteria),
          "({:.2f}%)".format(100.0 * len(bacteria) / total))
    print("Archaea >=97%:", len(archaea))
    print("Other assigned >=97%:", len(other))
    print("Unassigned:", len(unassigned),
          "({:.2f}%)".format(100.0 * len(unassigned) / total))
    print("Wolbachia:", len(wolbachia),
          "({:.2f}% of raw OTUs)".format(
              100.0 * len(wolbachia) / total
          ))
    print("Unassigned FASTA sequences:",
          sum(1 for line in open(unassigned_fasta)
              if line.startswith(">")))
    print("Summary:", summary_file)
    print("Unassigned FASTA:", unassigned_fasta)

