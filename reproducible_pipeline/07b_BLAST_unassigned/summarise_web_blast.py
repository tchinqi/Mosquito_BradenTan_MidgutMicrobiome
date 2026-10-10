#!/usr/bin/env python3
"""
summarise_web_blast.py - summarise NCBI web BLAST results (Single-file JSON)
for the 188 MC OTUs left unassigned by Greengenes2 (>=80%).

For each OTU: best hit (highest bit score), % identity, query coverage,
taxonomic category from the full NCBI lineage (taxonomy4blast.sqlite3), and
the OTU's supporting reads / number of samples from the MC OTU table.

Usage (in 08_unassigned_blast/):
    python3 summarise_web_blast.py web_blast.json
Outputs:
    unassigned_188_top_hits_web.tsv   (one row per OTU)
    printed summary by category + top 15 OTUs by reads
"""
import json, sqlite3, sys, os, csv
from collections import OrderedDict, defaultdict

JSON_IN   = sys.argv[1] if len(sys.argv) > 1 else "web_blast.json"
FASTA     = "MC_GG80_Unassigned_plus_domain_only.fasta"
TAXDB     = "taxdb/taxonomy4blast.sqlite3"
OTU_TABLE = ("../07_analysis/MC_GG2_97_classified_OTU_table/"
             "MC_867_OTUs_x_576_samples_supporting_reads.tsv")
OUT       = "unassigned_188_top_hits_web.tsv"
CHIMERA_QCOV = 70.0        # best hit covers <70% of the OTU -> possible chimera

# ---- 1. OTU names in the query FASTA ---------------------------------------
def first_token(s):
    return s.strip().lstrip(">").split()[0].split(";")[0]

fasta_ids = [first_token(l) for l in open(FASTA) if l.startswith(">")]
print("OTUs in FASTA:", len(fasta_ids))

# ---- 2. Parse BLAST JSON ----------------------------------------------------
data = json.load(open(JSON_IN))
reports = data["BlastOutput2"]
if isinstance(reports, dict):
    reports = [reports]

best = {}
for r in reports:
    search = r["report"]["results"]["search"]
    qtitle = search.get("query_title") or search.get("query_id")
    qid = first_token(qtitle)
    qlen = search["query_len"]
    top = None
    for h in search.get("hits", []):
        hsps = h["hsps"]
        bits = max(x["bit_score"] for x in hsps)
        if top is None or bits > top[0]:
            top = (bits, h)
    if top is None:
        best[qid] = None
        continue
    bits, h = top
    hsps = sorted(h["hsps"], key=lambda x: -x["bit_score"])
    hsp = hsps[0]
    # query coverage = union of all HSP query ranges for this subject / query length
    ranges = sorted((min(x["query_from"], x["query_to"]), max(x["query_from"], x["query_to"]))
                    for x in h["hsps"])
    covered, cs, ce = 0, None, None
    for s, e in ranges:
        if cs is None or s > ce + 1:
            if cs is not None:
                covered += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    covered += ce - cs + 1
    d = h["description"][0]
    best[qid] = dict(
        accession=d.get("accession", d.get("id", "")),
        taxid=int(d.get("taxid", 0) or 0),
        sciname=d.get("sciname", ""),
        title=d.get("title", ""),
        pident=round(100.0 * hsp["identity"] / hsp["align_len"], 2),
        align_len=hsp["align_len"],
        qcov=round(100.0 * covered / qlen, 1),
        evalue=hsp["evalue"],
        bitscore=bits,
        qlen=qlen)

missing = [q for q in fasta_ids if q not in best]
extra = [q for q in best if q not in set(fasta_ids)]
print("OTUs in BLAST JSON:", len(best), "| missing from JSON:", len(missing),
      "| names not matching FASTA:", len(extra))
if extra[:5]:
    print("  e.g. unmatched names:", extra[:5])

# ---- 3. Lineage from taxonomy4blast.sqlite3 ---------------------------------
con = sqlite3.connect(TAXDB)
tab = None
for (name,) in con.execute("SELECT name FROM sqlite_master WHERE type='table'"):
    cols = [c[1].lower() for c in con.execute("PRAGMA table_info(%s)" % name)]
    if "taxid" in cols and "parent" in cols:
        tab = name
        break
if tab is None:
    sys.exit("Could not find a taxid/parent table in " + TAXDB)

def lineage(t):
    seen = []
    while t and t not in seen and len(seen) < 100:
        seen.append(t)
        row = con.execute("SELECT parent FROM %s WHERE taxid=?" % tab, (t,)).fetchone()
        if row is None or row[0] == t:
            break
        t = row[0]
    return set(seen)

CATS = OrderedDict([            # first match wins
    (9606,  "Human"),
    (7157,  "Mosquito (Culicidae)"),
    (6656,  "Other arthropod"),
    (7742,  "Other vertebrate"),
    (4751,  "Fungi"),
    (33090, "Plant/alga (Viridiplantae)"),
    (2759,  "Other eukaryote"),
    (2157,  "Archaea"),
    (10239, "Virus"),
    (2,     "Bacteria"),
])

def category(b):
    if b is None:
        return "No hit"
    lin = lineage(b["taxid"])
    if len(lin) <= 1:
        return "Taxid not in local DB"
    for tid, lab in CATS.items():
        if tid in lin:
            if lab == "Bacteria":
                return "Bacteria (>=97% identity)" if b["pident"] >= 97 else "Bacteria (<97% identity)"
            return lab
    return "Other/unclassified"

# ---- 4. Reads and samples per OTU -------------------------------------------
reads, nsamp = {}, {}
with open(OTU_TABLE) as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd)
    for row in rd:
        vals = [float(x) if x not in ("", "NA") else 0.0 for x in row[1:]]
        otu = first_token(row[0])
        reads[otu] = sum(vals)
        nsamp[otu] = sum(1 for v in vals if v > 0)
print("Samples (columns) in OTU table:", len(header) - 1)
nomatch = [q for q in fasta_ids if q not in reads]
if nomatch:
    print("WARNING:", len(nomatch), "query OTUs not found in OTU table, e.g.", nomatch[:5])

# ---- 5. Write per-OTU table ---------------------------------------------------
rows = []
for q in fasta_ids:
    b = best.get(q)
    cat = category(b)
    chim = "yes" if (b is not None and b["qcov"] < CHIMERA_QCOV) else ""
    rows.append(dict(OTU=q, category=cat, reads=reads.get(q, 0), samples=nsamp.get(q, 0),
                     possible_chimera=chim,
                     **({k: b[k] for k in ("sciname", "taxid", "pident", "qcov", "align_len",
                                            "qlen", "evalue", "bitscore", "accession", "title")}
                        if b else {})))
cols = ["OTU", "category", "reads", "samples", "sciname", "taxid", "pident", "qcov",
        "possible_chimera", "align_len", "qlen", "evalue", "bitscore", "accession", "title"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in sorted(rows, key=lambda r: -r["reads"]):
        w.writerow(r)
print("Wrote", OUT)

# ---- 6. Summary ---------------------------------------------------------------
tot_reads = sum(r["reads"] for r in rows)
agg = defaultdict(lambda: [0, 0.0, 0])
for r in rows:
    a = agg[r["category"]]
    a[0] += 1
    a[1] += r["reads"]
    a[2] += 1 if r["possible_chimera"] else 0
print("\nSummary of %d unassigned OTUs (%.0f supporting reads)" % (len(rows), tot_reads))
print("%-30s %6s %10s %8s %16s" % ("Category", "OTUs", "Reads", "% reads", "qcov<70% (chim?)"))
for cat, (n, rd_, ch) in sorted(agg.items(), key=lambda kv: -kv[1][1]):
    print("%-30s %6d %10.0f %7.1f%% %16d" % (cat, n, rd_, 100 * rd_ / tot_reads if tot_reads else 0, ch))

print("\nTop 15 OTUs by supporting reads")
print("%-14s %8s %7s %-28s %6s %6s  %s" % ("OTU", "reads", "samples", "category", "pid", "qcov", "best hit"))
for r in sorted(rows, key=lambda r: -r["reads"])[:15]:
    print("%-14s %8.0f %7d %-28s %6s %6s  %s" % (r["OTU"], r["reads"], r["samples"], r["category"][:28],
          r.get("pident", ""), r.get("qcov", ""), r.get("sciname", "")[:45]))
