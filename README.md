# LSM4288 FYP – Midgut bacterial communities of wild-caught and laboratory-reared mosquitoes in Singapore

**Author:** Braden Tan (Tan Chin Qi Braden), National University of Singapore  
**Supervisor:** Nalini Puniamoorthy (Puniamoorthy lab) · **Mentor:** Priya  
**Data:** 576 individually sequenced midguts (570 mosquitoes + 6 negative controls), 16S rRNA V3–V4, Oxford Nanopore MinION

This repository contains all code used for the LSM4288 Final Year Project report: the
bioinformatics pipeline run on the NUS **Atlas9** HPC cluster, and the statistics and
figures produced locally in **R**.

## Where each part was run

| Folder | Run on | Contents |
|---|---|---|
| `reproducible_pipeline/01_demultiplexing` … `09_metadata_QC` | Atlas9 | Final pipeline, from raw reads to the OTU table and R inputs ("Scripts 01–53" in the report) |
| `reproducible_pipeline/07b_BLAST_unassigned` | NCBI web BLAST; summarised on Atlas9 | BLAST of the 188 OTUs Greengenes2 could not classify (Table S2E) |
| `R_analysis/` | Local (R 4.5.2) | Statistics and all report figures ("R scripts 00–12, 11a") |
| `reconstructed_code/`, `recovered_dependencies/` | Atlas9 | Reconstructed steps (labelled RECONSTRUCTED) and the patched NGSpeciesID module |
| `original_code/` | Atlas9 | Unedited original scripts (provenance; do not edit) |
| `development_and_superseded/` | Atlas9 | Test, diagnostic and superseded scripts kept for transparency |
| `documentation/` | – | Software versions, recovery notes and the master provenance ledger |

Every pipeline folder and `R_analysis/` has its own `README.md` describing each script, and
every curated pipeline script starts with a header (step, purpose, inputs, outputs). These
headers were added for this release as comment lines only; the code itself is unchanged,
and the untouched originals are in `original_code/`.

## Workflow and report sections

| Step | Folder | Report |
|---|---|---|
| Demultiplexing, 350–550 bp filtering | 01, 02 | Methods 2.5 |
| NGSpeciesID + Racon consensus, pooling | 03 | 2.6 |
| Abundance filtering (LC 0.1%; MC max(0.5%, 20 reads)), VSEARCH 97% | 04, 05 | 2.7 |
| Greengenes2 classification (≥80% / ≥97%) | 06 | 2.8, 3.2 |
| Kraken2 and BLAST cross-checks | 07, 07b | 2.9, 3.2 |
| OTU table, metadata, negative controls | 08, 09 | 2.10 |
| Statistics and figures | `R_analysis/` | 2.11, 3.1–3.5 |

Raw reads, intermediate files and COI identification spreadsheets are not included.

## Original repository overview
Code repository for the LSM4288 Final Year Project analysing bacterial
communities in mosquito midgut samples using Nanopore 16S V3–V4 sequencing.

## Current repository status

**Final release (2026-10-10).** This repository contains the frozen and
provenance-audited upstream bioinformatics workflow (run on NUS Atlas9) and
the final R statistical analysis and figure scripts (`R_analysis/`, run
locally), including wild host-species assignment from COI barcodes.

## Final upstream dataset

The final midgut dataset contains 576 samples across six 96-well plates:

- AG_P01
- AG_P02
- ZX_P01
- ZX_P02
- BT_P01
- BT_P02

Final metadata composition after negative-control correction
(script 49, corrected by script 53 on 2026-10-07):

- 380 Wild
- 190 Lab-reared
- 6 negative controls
- 570 biological samples

AG_Repro reproductive-organ samples are excluded from the final midgut
analysis.

## Main workflow

1. Demultiplex reads using Minibar
2. Apply 350–550 bp read-length filtering
3. Generate within-sample consensus sequences using NGSpeciesID
4. Exclude singleton clusters before Racon polishing
5. Rescue samples containing individual Racon-polishing failures
6. Pool 22,528 polished consensus sequences
7. Apply post-Racon abundance-support thresholds
8. Cluster retained consensuses globally with VSEARCH at 97% identity
9. Classify OTUs using Greengenes2
10. Cross-check unresolved Greengenes2 assignments using Kraken2
10b. BLAST the 188 unresolved MC OTUs against NCBI core_nt (web BLAST; `reproducible_pipeline/07b_BLAST_unassigned/`)
11. Construct the final MC + Greengenes2 ≥97% abundance table
12. Join sample metadata and correct six negative controls
13. Generate and verify final analysis-ready inputs
14. Correct the BT_P02 negative control (script 53) and regenerate inputs (45, 50, 51, 52)
15. Assign wild host species from COI barcodes (`R_analysis/11a_build_host_species.R`, in R)
16. Statistical analyses and figures in R (`R_analysis/00`–`12`; see `R_analysis/README_R.md`)

## LC and MC terminology

Final terminology used in this repository:

- LC = 0.1% per-sample support threshold
- MC = max(0.5% of sample reads, 20 reads)

LC is the less conservative threshold and produced 7,642 OTUs.

MC is the more conservative threshold and produced 867 OTUs.

Historical scripts and directories occasionally used the opposite LC/MC
labels. Those historical names have been preserved for provenance and are
documented in `documentation/MASTER_RECOVERY_LEDGER.md`.

## Repository structure

`original_code/`
Complete preserved historical source code. Files here should not be edited.

`recovered_dependencies/`
Recovered external or locally modified dependencies, including the
NGSpeciesID Racon-failure continuation patch.

`reconstructed_code/`
Code reconstructed where the original interactive command was unavailable.
These files are explicitly labelled RECONSTRUCTED.

`reproducible_pipeline/`
Curated final upstream workflow.

`development_and_superseded/`
Historical development, diagnostic, intermediate and superseded scripts that
are preserved for provenance but are not part of the final execution path.

`R_analysis/`
Final R statistics and figure scripts (run locally, after Atlas9), including
host-species assignment (11a). See `R_analysis/README_R.md`.

`documentation/`
Recovery records, software versions, script accounting, provenance audits and
the master recovery ledger.

## Provenance categories

The repository distinguishes:

- ORIGINAL
- RECOVERED DEPENDENCY
- RECONSTRUCTED
- RECONSTRUCTED — BYTE-IDENTICAL
- SUPERSEDED
- CURRENT WORK

See:

`documentation/MASTER_RECOVERY_LEDGER.md`

for the authoritative provenance record.

## Important reconstruction note

Several historically interactive operations did not have standalone source
files. Where possible, these were reconstructed from preserved inputs, logs
and historical outputs.

The following reconstructions were validated against historical outputs:

- final Racon consensus pooling and provenance construction
- final LC VSEARCH97 clustering
- final MC VSEARCH97 clustering
- uppercase centroid FASTA conversion

These reconstructions reproduced their corresponding historical outputs
byte-for-byte.

## Software environment

Recorded Atlas software versions include:

- QIIME2 / q2cli 2026.7.0
- NGSpeciesID 0.3.1
- VSEARCH 2.31.0
- Kraken2 2.17.1
- SeqKit 2.13.0
- Racon 1.5.0
- SPOA 4.1.5

See:

`documentation/ATLAS_SOFTWARE_VERSIONS_20261006.txt`

## Historical code versus final workflow

The complete final-reanalysis snapshot contains 55 scripts.

Of these:

- 23 are directly promoted into the curated reproducible pipeline
- 32 are retained as historical non-pipeline scripts

Additional verified upstream originals and reconstructed scripts are also
included in the curated pipeline.

No historical source files were deleted during repository construction.

## Data

Large sequencing data, intermediate FASTA files, QIIME2 artifacts, Kraken2
databases and reconstructed validation outputs are not intended to be stored
in the Git repository.

The repository focuses on source code, provenance and reproducibility
documentation.
