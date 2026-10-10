# Demultiplexing + 350–550 bp length filtering — CLOSED
Date: 2026-10-06

Status:
ORIGINAL HISTORICAL CODE RECOVERED.

Recovered original source files:
- scripts/make_demux_braden.py
- scripts/02_run_minibar.pbs
- scripts/07_length_filter.pbs

Historical dated copies:
logs/scripts_snapshot_20260820_FINAL/

Current scripts are byte-identical to the 2026-08-20 FINAL snapshot.

SHA-256:

02_run_minibar.pbs
f55155bd999e59fd46a5e673ff8f1fe7b8699d04f4fb14d8aec5007d82a1f24e

make_demux_braden.py
2dd4badbf017562b3be8f53b8e35bdcd166c1db1d7754e9fe5c35e6a5cbe6f70

07_length_filter.pbs
d7562a40f368ae6202a03cf49882fafea913abb614485172ffbdb788e44cd380

Demultiplexing:
- Minibar workflow
- demultiplexing sheet generated from original metadata
- primers:
  338F = ACTCCTACGGGAGGCAGCAG
  806R = GGACTACHVGGGTWTCTAAT

Production length filter:
seqkit seq -g -m 350 -M 550

Input directory:
03_no_gaps_fastq

Output directory:
03b_length_filtered_fastq

Historical totals used downstream:
- 2,289,751 demultiplexed reads
- 1,795,316 reads passed the 350–550 bp filter
