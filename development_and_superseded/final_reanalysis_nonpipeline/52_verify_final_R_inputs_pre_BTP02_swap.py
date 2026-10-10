#!/usr/bin/env python3
"""
52_verify_final_R_inputs.py
LSM4288 FYP - pre-R sanity check of the final MC + GG2 >=97% plotting files.

Read-only: it never modifies any input file. It checks every file R will use
against the agreed source-of-truth numbers and prints PASS / FAIL / WARN.

Outputs (written to $RUN/07_analysis/MC_GG2_97_plotting_final/verification_52/):
  verification_report.txt          - everything printed to screen
  per_sample_classified_depth.tsv  - sample, origin, location, classified reads,
                                     observed OTUs (for the depth/rarefaction check)

Usage on Atlas9:
  RUN=/hpctmp/e0958567/Braden_16S_20260819/20_FINAL_reanalysis_20260903
  /hpctmp/e0958567/conda_envs/ngspeciesid-0.3.1/bin/python \
      $RUN/scripts/52_verify_final_R_inputs.py --run "$RUN"

Standard library only (no pandas needed).
"""
import argparse
import csv
import hashlib
import os
import re
import statistics
import sys

# ---------------------------------------------------------------------------
# Source-of-truth expectations
# ---------------------------------------------------------------------------
EXP = dict(
    total_samples=576, wild=380, lab=190, controls=6, biological=570,
    location_samples=345, otus=328, classified_reads=1_432_996,
)
EXP_LOCATIONS = {"Mandai Zoo": 194, "Kranji Marshes": 88,
                 "Bird Paradise": 39, "River Safari": 24}
EXP_CONTROLS = {
    "sample_BT_CLX_096_16S_BT_P01_16S",
    "sample_BT_CLX_120_16S_BT_P02_16S",
    "sample_Neg_AG_P01_16S",
    "sample_Negative_AG_P02_16S",
    "sample_Negative_ZX_P01_ZX_P01_16S",
    "sample_Negative_ZX_P02_16S",
}
# Descriptive means from script 51 (percent), used as a reproducibility check
EXP_GENUS_MEANS = {
    "Elizabethkingia": (25.62, 15.43, 45.99),
    "Wolbachia": (20.49, 25.96, 9.55),
}

# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------
LOG = []
COUNTS = {"PASS": 0, "FAIL": 0, "WARN": 0}


def say(msg=""):
    print(msg)
    LOG.append(msg)


def check(ok, label, detail=""):
    tag = "PASS" if ok else "FAIL"
    COUNTS[tag] += 1
    say(f"  [{tag}] {label}" + (f"  ({detail})" if detail else ""))
    return ok


def warn(label, detail=""):
    COUNTS["WARN"] += 1
    say(f"  [WARN] {label}" + (f"  ({detail})" if detail else ""))


def header(title):
    say("")
    say("=" * 78)
    say(title)
    say("=" * 78)


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty file: {path}")
    return rows[0], rows[1:]


def find_col(header_row, *patterns, exclude=()):
    """First column whose name matches any regex (case-insensitive)."""
    for i, name in enumerate(header_row):
        low = name.lower()
        if any(re.search(e, low) for e in exclude):
            continue
        if any(re.search(p, low) for p in patterns):
            return i
    return None


def to_num(x):
    x = x.strip()
    if x in ("", "NA", "nan", "NaN"):
        return 0.0
    return float(x)


def norm_origin(v):
    v = v.strip().lower()
    if v.startswith("wild"):
        return "Wild"
    if v.startswith("lab"):
        return "Lab"
    if "control" in v or v.startswith("neg"):
        return "Control"
    return v or "<blank>"


# ---------------------------------------------------------------------------
# Table loader: samples x OTUs (tolerates extra metadata columns)
# ---------------------------------------------------------------------------
def load_sample_table(path, feature_ids=None, feature_regex=r"^OTU_\d+$"):
    hdr, rows = read_tsv(path)
    if feature_ids:
        feat_idx = [i for i, h in enumerate(hdr) if h in feature_ids]
    else:
        feat_idx = [i for i, h in enumerate(hdr) if re.match(feature_regex, h)]
    data = {}
    dup = 0
    for r in rows:
        if not r or not r[0].strip():
            continue
        sid = r[0].strip()
        if sid in data:
            dup += 1
        data[sid] = [to_num(r[i]) if i < len(r) else 0.0 for i in feat_idx]
    feats = [hdr[i] for i in feat_idx]
    return hdr, feats, data, dup


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=os.environ.get("RUN"),
                    help="Path to 20_FINAL_reanalysis_20260903")
    args = ap.parse_args()
    if not args.run:
        sys.exit("Give --run or export RUN=...")
    RUN = args.run.rstrip("/")

    P = f"{RUN}/07_analysis/MC_GG2_97_plotting_final"
    files = {
        "metadata": f"{RUN}/metadata/sample_metadata_current.tsv",
        "bio": f"{P}/MC_GG97_biological_570_samples_x_328_OTUs.tsv",
        "wild": f"{P}/MC_GG97_wild_380_samples_x_328_OTUs.tsv",
        "lab": f"{P}/MC_GG97_lab_190_samples_x_328_OTUs.tsv",
        "loc": f"{P}/MC_GG97_location_345_samples_x_328_OTUs.tsv",
        "ctrl": f"{P}/MC_GG97_controls_6_samples_x_328_OTUs.tsv",
        "genus_counts": f"{P}/MC_GG97_biological_570_samples_genus_counts.tsv",
        "genus_ra": f"{P}/MC_GG97_biological_570_samples_genus_relative_abundance.tsv",
        "otu2genus": f"{P}/MC_GG97_328_OTU_to_genus.tsv",
        "wvl": f"{P}/wild_vs_lab_genus/MC_GG97_wild_vs_lab_genus_mean_relative_abundance_FINAL.tsv",
        "joined": f"{RUN}/07_analysis/MC_GG2_97_analysis_ready/"
                  "MC_GG2_97_classified_576_samples_x_328_OTUs_with_metadata.tsv",
    }
    outdir = f"{P}/verification_52"

    # -------------------------------------------------------------- 0. files
    header("0. FILES PRESENT")
    missing = []
    for key, path in files.items():
        if os.path.isfile(path):
            say(f"  [ OK ] {key:13s} {os.path.getsize(path):>12,d} B  md5 {md5(path)}")
        else:
            missing.append(key)
            check(False, f"{key} missing", path)
    if any(k in missing for k in ("metadata", "bio", "otu2genus")):
        say("\nCore files missing - stopping.")
        finish(outdir)
        return

    old = f"{RUN}/07_analysis/MC_GG2_97_plotting"
    if os.path.isdir(old):
        warn("superseded 4-control folder still exists - do NOT load it in R", old)

    # ----------------------------------------------------------- 1. metadata
    header("1. METADATA  (sample_metadata_current.tsv)")
    mh, mrows = read_tsv(files["metadata"])
    say(f"  columns: {', '.join(mh)}")
    c_id = 0
    c_origin = find_col(mh, r"^origin$", r"origin")
    c_incl = find_col(mh, r"include.*location")
    c_loc = find_col(mh, r"location", exclude=(r"status", r"include", r"code"))
    c_well = find_col(mh, r"well")
    if c_origin is None:
        check(False, "Origin column found")
        finish(outdir)
        return

    meta = {}
    for r in mrows:
        if r and r[0].strip():
            meta[r[0].strip()] = r
    check(len(mrows) == EXP["total_samples"] and len(meta) == EXP["total_samples"],
          f"{EXP['total_samples']} rows, all sample IDs unique",
          f"rows={len(mrows)}, unique={len(meta)}")

    origin = {s: norm_origin(r[c_origin]) for s, r in meta.items()}
    oc = {}
    for v in origin.values():
        oc[v] = oc.get(v, 0) + 1
    say(f"  Origin counts: {oc}")
    check(oc.get("Wild", 0) == EXP["wild"], f"Wild = {EXP['wild']}", f"got {oc.get('Wild', 0)}")
    check(oc.get("Lab", 0) == EXP["lab"], f"Lab = {EXP['lab']}", f"got {oc.get('Lab', 0)}")
    check(oc.get("Control", 0) == EXP["controls"], f"Controls = {EXP['controls']}",
          f"got {oc.get('Control', 0)}")

    ctrl_ids = {s for s, o in origin.items() if o == "Control"}
    check(ctrl_ids == EXP_CONTROLS, "control IDs are exactly the six H12 wells",
          f"extra={sorted(ctrl_ids - EXP_CONTROLS)} missing={sorted(EXP_CONTROLS - ctrl_ids)}")
    if c_well is not None:
        bad = [s for s in EXP_CONTROLS if s in meta and meta[s][c_well].strip().upper() != "H12"]
        check(not bad, "all controls have well = H12", f"not H12: {bad}")

    bio_ids = {s for s, o in origin.items() if o in ("Wild", "Lab")}
    wild_ids = {s for s, o in origin.items() if o == "Wild"}
    lab_ids = {s for s, o in origin.items() if o == "Lab"}

    # Lab IDs should be BT_ALB / BT_CLX; wild should be AG / ZX
    lab_bad = [s for s in lab_ids if not re.search(r"BT_(ALB|CLX)", s)]
    check(not lab_bad, "all Lab IDs are BT_ALB_* / BT_CLX_*", f"{len(lab_bad)} odd: {lab_bad[:5]}")
    wild_bad = [s for s in wild_ids if not re.search(r"_(AG|ZX)_P0[12]", s)]
    if wild_bad:
        warn("Wild IDs not matching AG/ZX plate pattern (check by eye)",
             f"{len(wild_bad)}: {wild_bad[:5]}")
    alb = sum(1 for s in lab_ids if "BT_ALB" in s)
    clx = sum(1 for s in lab_ids if "BT_CLX" in s)
    say(f"  Lab split: BT_ALB = {alb}, BT_CLX = {clx}")

    loc_ids = set()
    locmap = {}
    if c_incl is not None:
        loc_ids = {s for s, r in meta.items() if r[c_incl].strip().upper() in ("YES", "Y", "TRUE", "1")}
        check(len(loc_ids) == EXP["location_samples"],
              f"Include_in_Location_Analysis = YES for {EXP['location_samples']}",
              f"got {len(loc_ids)}")
        check(loc_ids <= wild_ids, "all location-included samples are Wild",
              f"non-wild: {sorted(loc_ids - wild_ids)[:5]}")
    else:
        warn("no Include_in_Location_Analysis column found")
    if c_loc is not None:
        say(f"  location column used: '{mh[c_loc]}'")
        lc = {}
        for s in loc_ids:
            v = meta[s][c_loc].strip()
            locmap[s] = v
            lc[v] = lc.get(v, 0) + 1
        say(f"  Location counts (included samples): {lc}")
        for name, n in EXP_LOCATIONS.items():
            got = sum(c for k, c in lc.items()
                      if k.lower() == name.lower()
                      or k.upper() == "".join(w[0] for w in name.split()))
            check(got == n, f"{name} = {n}", f"got {got}")
    else:
        warn("could not identify the location column - check location counts by hand")

    c_species = find_col(mh, r"species", r"taxon", r"host")
    if c_species is None:
        say("  (no species column yet - expected, waiting for Wild species metadata)")
    else:
        say(f"  species-like column present: '{mh[c_species]}'")

    # --------------------------------------------------------- 2. OTU->genus
    header("2. OTU -> GENUS MAP")
    gh, grows = read_tsv(files["otu2genus"])
    say(f"  columns: {', '.join(gh)}")
    c_gen = find_col(gh, r"genus")
    otu_ids = [r[0].strip() for r in grows if r and r[0].strip()]
    check(len(otu_ids) == EXP["otus"] and len(set(otu_ids)) == EXP["otus"],
          f"{EXP['otus']} unique OTUs", f"rows={len(otu_ids)}, unique={len(set(otu_ids))}")
    otu2gen = {}
    if c_gen is not None:
        otu2gen = {r[0].strip(): r[c_gen].strip() for r in grows if r and r[0].strip()}
        genera = sorted(set(otu2gen.values()))
        say(f"  distinct genus labels: {len(genera)}  (expected 68)")
        check(len(genera) == 68, "68 genus labels", f"got {len(genera)}")
        unc = [g for g in genera if "unclass" in g.lower()]
        say(f"  unclassified-type labels: {unc}")
        blank = sum(1 for g in otu2gen.values() if g in ("", "NA"))
        check(blank == 0, "no blank genus labels", f"{blank} blank")
    else:
        warn("no genus column found in OTU map")
    otu_set = set(otu_ids)

    # ---------------------------------------------------- 3. biological table
    header("3. BIOLOGICAL 570 x 328 TABLE")
    _, feats, bio, dup = load_sample_table(files["bio"], feature_ids=otu_set)
    check(dup == 0, "no duplicate sample rows", f"{dup} duplicates")
    check(len(bio) == EXP["biological"], f"{EXP['biological']} samples", f"got {len(bio)}")
    check(len(feats) == EXP["otus"], f"{EXP['otus']} OTU columns", f"got {len(feats)}")
    check(set(feats) == otu_set, "OTU columns match OTU->genus map exactly")
    check(set(bio) == bio_ids, "sample IDs = Wild + Lab in metadata exactly",
          f"extra={sorted(set(bio) - bio_ids)[:5]} missing={sorted(bio_ids - set(bio))[:5]}")
    check(not (set(bio) & EXP_CONTROLS), "no negative controls in biological table")

    allvals = [v for row in bio.values() for v in row]
    check(all(v >= 0 for v in allvals), "no negative counts")
    check(all(abs(v - round(v)) < 1e-9 for v in allvals), "all counts are integers")
    depth = {s: sum(row) for s, row in bio.items()}
    rich = {s: sum(1 for v in row if v > 0) for s, row in bio.items()}
    zero = [s for s, d in depth.items() if d == 0]
    check(not zero, "every biological sample has >=1 classified read", f"{len(zero)} empty")
    empty_otus = [feats[j] for j in range(len(feats)) if all(row[j] == 0 for row in bio.values())]
    if empty_otus:
        warn("OTUs with zero reads in all 570 biological samples (present only in controls?)",
             f"{len(empty_otus)}: {empty_otus[:8]}")
    bio_total = int(round(sum(depth.values())))
    say(f"  total classified reads in 570 biological samples: {bio_total:,}")

    # ------------------------------------------------------------ 4. controls
    header("4. CONTROLS TABLE + READ TOTAL")
    ctrl_total = None
    if "ctrl" not in missing:
        _, cf, ctrl, _ = load_sample_table(files["ctrl"], feature_ids=otu_set)
        check(set(ctrl) == EXP_CONTROLS, "controls table holds exactly the 6 controls",
              f"got {sorted(ctrl)}")
        ctrl_total = int(round(sum(sum(r) for r in ctrl.values())))
        for s in sorted(ctrl):
            say(f"    {s:40s} {int(sum(ctrl[s])):>7,d} reads, {sum(1 for v in ctrl[s] if v > 0)} OTUs")
        check(bio_total + ctrl_total == EXP["classified_reads"],
              f"biological + controls = {EXP['classified_reads']:,}",
              f"{bio_total:,} + {ctrl_total:,} = {bio_total + ctrl_total:,}")

    # ------------------------------------------------------- 5. subset tables
    header("5. SUBSET TABLES (must be exact slices of the 570 table)")
    for key, ids, n in (("wild", wild_ids, EXP["wild"]),
                        ("lab", lab_ids, EXP["lab"]),
                        ("loc", loc_ids, EXP["location_samples"])):
        if key in missing:
            continue
        _, sf, sub, _ = load_sample_table(files[key], feature_ids=otu_set)
        check(len(sub) == n, f"{key}: {n} samples", f"got {len(sub)}")
        if ids:
            check(set(sub) == ids, f"{key}: IDs match metadata",
                  f"extra={sorted(set(sub) - ids)[:3]} missing={sorted(ids - set(sub))[:3]}")
        order = {f: i for i, f in enumerate(feats)}
        mism = 0
        for s, row in sub.items():
            if s not in bio:
                mism += 1
                continue
            ref = bio[s]
            if any(abs(row[i] - ref[order[f]]) > 1e-9 for i, f in enumerate(sf)):
                mism += 1
        check(mism == 0, f"{key}: values identical to 570 table", f"{mism} samples differ")

    # --------------------------------------------------------- 6. genus tables
    header("6. GENUS COUNT / RELATIVE-ABUNDANCE TABLES")
    gc_ok = False
    if "genus_counts" not in missing and otu2gen:
        gch, gcrows = read_tsv(files["genus_counts"])
        gcols = [i for i, h in enumerate(gch) if h in set(otu2gen.values())]
        say(f"  genus columns recognised: {len(gcols)}")
        check(len(gcols) == len(set(otu2gen.values())), "genus columns = labels in OTU map",
              f"{len(gcols)} vs {len(set(otu2gen.values()))}")
        gcount = {r[0].strip(): {gch[i]: to_num(r[i]) for i in gcols} for r in gcrows if r and r[0].strip()}
        check(set(gcount) == bio_ids, "genus counts: 570 biological samples")
        # recompute from OTU table
        mism = 0
        for s, row in bio.items():
            agg = {}
            for f, v in zip(feats, row):
                g = otu2gen[f]
                agg[g] = agg.get(g, 0) + v
            got = gcount.get(s, {})
            if any(abs(agg.get(g, 0) - got.get(g, 0)) > 1e-6 for g in set(agg) | set(got)):
                mism += 1
        gc_ok = check(mism == 0, "genus counts = OTU table aggregated by genus", f"{mism} samples differ")

    if "genus_ra" not in missing:
        rh, rrows = read_tsv(files["genus_ra"])
        rcols = [i for i, h in enumerate(rh) if not otu2gen or h in set(otu2gen.values())]
        ra = {r[0].strip(): {rh[i]: to_num(r[i]) for i in rcols} for r in rrows if r and r[0].strip()}
        sums = [sum(v.values()) for v in ra.values()]
        scale = 100.0 if statistics.median(sums) > 1.5 else 1.0
        say(f"  relative abundance scale detected: {'percent (0-100)' if scale == 100 else 'fraction (0-1)'}")
        check(all(abs(x - scale) < 1e-3 * scale for x in sums), "every sample sums to 1 / 100",
              f"min={min(sums):.4f} max={max(sums):.4f}")
        check(set(ra) == bio_ids, "rel. abundance: 570 biological samples")

        # reproduce script-51 means
        say("  reproducing script-51 mean relative abundance (%): overall / wild / lab")
        for g, (eo, ew, el) in EXP_GENUS_MEANS.items():
            def mean_for(ids):
                vals = [ra[s].get(g, 0) * (100 / scale) for s in ids if s in ra]
                return sum(vals) / len(vals) if vals else float("nan")
            o, w, l = mean_for(bio_ids), mean_for(wild_ids), mean_for(lab_ids)
            ok = all(abs(a - b) < 0.02 for a, b in ((o, eo), (w, ew), (l, el)))
            check(ok, f"{g}: {o:.2f} / {w:.2f} / {l:.2f}", f"expected {eo} / {ew} / {el}")

    # -------------------------------------------------- 7. joined 576 table
    header("7. ANALYSIS-READY JOINED TABLE (576 x 342)")
    if "joined" not in missing:
        jh, jrows = read_tsv(files["joined"])
        say(f"  {len(jrows)} rows x {len(jh)} columns")
        jo = find_col(jh, r"^origin$", r"origin")
        if jo is not None:
            jc = {}
            for r in jrows:
                v = norm_origin(r[jo])
                jc[v] = jc.get(v, 0) + 1
            if jc.get("Control", 0) == 6:
                check(True, "joined table already has 6 controls", str(jc))
            else:
                warn("joined table pre-dates the 6-control fix - in R, take Origin/Location "
                     "from sample_metadata_current.tsv, NOT from this file", str(jc))

    # ------------------------------------------------- 8. depth (for R plan)
    header("8. PER-SAMPLE CLASSIFIED DEPTH (informs the rarefaction decision)")

    def summ(vals):
        v = sorted(vals)
        q = statistics.quantiles(v, n=4) if len(v) >= 4 else [v[0], statistics.median(v), v[-1]]
        return f"n={len(v)}  min={v[0]:,.0f}  Q1={q[0]:,.0f}  median={statistics.median(v):,.0f}  " \
               f"Q3={q[2]:,.0f}  max={v[-1]:,.0f}"

    say(f"  All : {summ(depth.values())}")
    say(f"  Wild: {summ([depth[s] for s in wild_ids])}")
    say(f"  Lab : {summ([depth[s] for s in lab_ids])}")
    for cut in (100, 250, 500, 1000, 2000):
        nw = sum(1 for s in wild_ids if depth[s] < cut)
        nl = sum(1 for s in lab_ids if depth[s] < cut)
        say(f"  samples < {cut:>5,d} classified reads: Wild {nw:>3d}, Lab {nl:>3d}, total {nw + nl:>3d}")
    say(f"  Observed OTUs  All : {summ(rich.values())}")
    say(f"  Observed OTUs  Wild: {summ([rich[s] for s in wild_ids])}")
    say(f"  Observed OTUs  Lab : {summ([rich[s] for s in lab_ids])}")

    # Spearman rho between log depth and richness (no scipy needed)
    def rank(x):
        order = sorted(range(len(x)), key=lambda i: x[i])
        r = [0.0] * len(x)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    ids = sorted(bio)
    rd, rr = rank([depth[s] for s in ids]), rank([rich[s] for s in ids])
    rho = statistics.correlation(rd, rr) if hasattr(statistics, "correlation") else float("nan")
    say(f"  Spearman rho (classified depth vs observed OTUs), 570 samples: {rho:.3f}")

    os.makedirs(outdir, exist_ok=True)
    with open(f"{outdir}/per_sample_classified_depth.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id", "origin", "location", "classified_reads", "observed_otus"])
        for s in ids:
            w.writerow([s, origin[s], locmap.get(s, "NA"), int(round(depth[s])), rich[s]])
    say(f"\n  wrote {outdir}/per_sample_classified_depth.tsv")

    finish(outdir)


def finish(outdir):
    header("SUMMARY")
    say(f"  PASS {COUNTS['PASS']}   FAIL {COUNTS['FAIL']}   WARN {COUNTS['WARN']}")
    say("  All PASS -> safe to start R from MC_GG2_97_plotting_final/."
        if COUNTS["FAIL"] == 0 else "  Fix FAIL lines before starting R.")
    try:
        os.makedirs(outdir, exist_ok=True)
        with open(f"{outdir}/verification_report.txt", "w") as f:
            f.write("\n".join(LOG) + "\n")
        print(f"  report saved: {outdir}/verification_report.txt")
    except OSError as e:
        print(f"  could not save report: {e}")


if __name__ == "__main__":
    main()
