from pathlib import Path
import csv
from collections import Counter

datasets = {
    "MC_0p1pct": Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903/06_taxonomy/greengenes2_MC_LC_final/MC_0p1pct/exported_taxonomy/taxonomy.tsv"),
    "LC_max0p5pct_20reads": Path("/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903/06_taxonomy/greengenes2_MC_LC_final/LC_max0p5pct_20reads/exported_taxonomy/taxonomy.tsv"),
}

def parse_rank(taxon, prefix):
    parts = [x.strip() for x in taxon.split(";")]
    for p in parts:
        if p.startswith(prefix):
            val = p.split("__", 1)[1].strip()
            if val != "":
                return val
    return ""

for label, taxfile in datasets.items():
    outdir = taxfile.parent.parent
    assigned_rows = []
    genus_counter = Counter()
    species_counter = Counter()

    with open(taxfile, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            otu = row["Feature ID"]
            taxon = row["Taxon"]
            consensus = row["Consensus"]

            if taxon == "Unassigned":
                continue

            genus = parse_rank(taxon, "g__")
            species = parse_rank(taxon, "s__")

            assigned_rows.append({
                "Feature ID": otu,
                "Taxon": taxon,
                "Consensus": consensus,
                "Genus": genus,
                "Species": species
            })

            if genus:
                genus_counter[genus] += 1
            if species:
                species_counter[species] += 1

    assigned_path = outdir / "all_assigned_OTUs.tsv"
    genus_path = outdir / "all_genera_counts.tsv"
    species_path = outdir / "all_species_counts.tsv"

    with open(assigned_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Feature ID", "Taxon", "Consensus", "Genus", "Species"])
        for r in assigned_rows:
            w.writerow([r["Feature ID"], r["Taxon"], r["Consensus"], r["Genus"], r["Species"]])

    with open(genus_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Genus", "OTU_count"])
        for genus, n in sorted(genus_counter.items(), key=lambda x: (-x[1], x[0])):
            w.writerow([genus, n])

    with open(species_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Species", "OTU_count"])
        for species, n in sorted(species_counter.items(), key=lambda x: (-x[1], x[0])):
            w.writerow([species, n])

    print(f"===== {label} =====")
    print(f"Assigned OTUs: {len(assigned_rows)}")
    print(f"Genera: {len(genus_counter)}")
    print(f"Species: {len(species_counter)}")
    print(f"Output dir: {outdir}")
    print()

