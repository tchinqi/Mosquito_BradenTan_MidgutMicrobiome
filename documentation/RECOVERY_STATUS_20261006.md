# LSM4288 code/provenance recovery checkpoint
Date: 2026-10-06

## 1. Original exact script snapshot

Original saved-code snapshot:
LSM4288_EXACT_CODE_SNAPSHOT_20261006.tar.gz

Original snapshot contained:
- 55 saved scripts
- ANALYSIS_LOG.md
- script inventory
- SHA-256 inventory
- combined ALL_SCRIPTS_SOURCE.txt

Previously verified SHA-256:
72c4f7bf0db8b498a1d2e666e34d5d7e78ba426fa030d8000f931c52aceaa2c3

The original archive must remain unchanged.

---

## 2. External NGSpeciesID patch recovered

Patch:
software/ngspeciesid_min2_continue/

Scripts referencing this patch:
- 12_min2_racon_baseline_midgut.pbs
- 13_test_rescue_one_sample.pbs
- 14_rescue_8_min2_racon.pbs
- 15_rescue_0p1_0p2_thresholds.pbs

Patch contains 17 files, including 8 Python source files.

Important files:
- modules/consensus.py
- modules/consensus.py.original

The modified consensus.py:
1. wraps Racon polishing in try/except
2. checks that Racon produced a non-empty consensus
3. writes failed clusters to racon_failed_clusters.tsv
4. excludes only the failed cluster
5. continues processing remaining clusters
6. removes failed None entries before returning final consensuses

This modifies NGSpeciesID failure handling, not Racon itself.

Patch archive:
NGSPECIESID_MIN2_CONTINUE_PATCH_20261006.tar.gz

Previously verified SHA-256:
5952087ac928cd56b80d4bf1ce53c5425588b63e39c90b65cd59de443999717f

---

## 3. Final min2 pooled dataset

Final files:
04_global_vsearch97/min2_complete576/pooled_racon_polished_576.fasta
04_global_vsearch97/min2_complete576/pooled_racon_polished_576_provenance.tsv

Provenance columns, in actual file order:
- supporting_reads
- sequence_id
- cluster_id
- source
- sample_id

Validated:
- 22,528 provenance rows / polished consensuses
- 576 unique samples
- 22,181 ORIGINAL consensuses
- 347 RESCUE consensuses
- 568 ORIGINAL samples
- 8 RESCUE samples

Therefore:
568 normally completed + 8 rescued = 576/576 represented.

---

## 4. Eight rescue samples

sample_0820CW-3_AG_P01_16S
- 52 retained consensuses
- 4,609 supporting reads

sample_240926MZFF121_ZX_P01_16S
- 52 retained consensuses
- 3,571 supporting reads

sample_240926MZFF122_ZX_P01_16S
- 35 retained consensuses
- 2,295 supporting reads

sample_241206MZAF9_ZX_P02_16S
- 40 retained consensuses
- 1,270 supporting reads

sample_250220MZAB3_ZX_P02_16S
- 54 retained consensuses
- 2,443 supporting reads

sample_BT_ALB_003_16S_BT_P01_16S
- 37 retained consensuses
- 1,845 supporting reads

sample_BT_ALB_024_16S_BT_P01_16S
- 39 retained consensuses
- 1,998 supporting reads

sample_BT_ALB_080_16S_BT_P01_16S
- 38 retained consensuses
- 1,401 supporting reads

Total RESCUE retained supporting reads:
19,432

Rescue directory:
03_ngspeciesid_racon/min2_rescue/samples/

Baseline directory:
03_ngspeciesid_racon/min2_baseline/

---

## 5. Concrete rescue example recovered

Sample:
sample_0820CW-3_AG_P01_16S

The sample contains:
- final_clusters.tsv
- final_cluster_origins.tsv
- reads_to_consensus_<cluster>.fastq
- racon_cl_id_<cluster>/consensus.fasta for successful clusters
- racon_failed_clusters.tsv

Cluster 3072:
- 2 supporting reads
- entered Racon polishing
- Racon returned non-zero exit status 1
- no final racon_cl_id_3072/consensus.fasta was produced
- patched NGSpeciesID recorded the failure and continued

The final provenance contains 52 retained consensuses for this sample.

This provides direct evidence of the rescue mechanism.

---

## 6. VSEARCH97 provenance

VSEARCH version observed:
v2.31.0_linux_x86_64

Final unthresholded VSEARCH97 input:
04_global_vsearch97/min2_complete576/pooled_racon_polished_576.fasta

Parameters documented in ANALYSIS_LOG.md:
- --cluster_fast
- --id 0.97
- --iddef 2
- --strand both
- --relabel OTU_
- --threads 1

Result:
- 22,528 input polished consensuses
- 10,915 97% OTUs

The literal interactive shell command was NOT recovered from ~/.bash_history.

Therefore final VSEARCH command should later be labelled:
RECONSTRUCTED FROM CONTEMPORARY SCRIPTS + ANALYSIS LOG + FINAL OUTPUT FILES.

---

## 7. Historical MC/LC naming warning

Historical folder labels were reversed.

FINAL terminology:
- LC = 0.1% per-sample support threshold = 7,642 OTUs
- MC = max(0.5% of sample reads, 20 reads) = 867 OTUs

Historical:
- MC_0p1pct = current LC
- LC_max0p5pct_20reads = current MC

Never silently rename historical source files; document the correction in the reproducible pipeline.

---

## 8. Original code not yet recovered

The following historical operations are confirmed by outputs/logs but their literal original interactive code has not yet been recovered:

1. construction of pooled_racon_polished_576.fasta
2. construction of pooled_racon_polished_576_provenance.tsv
3. final VSEARCH97 interactive command
4. post-threshold LC/MC VSEARCH97 interactive commands
5. uppercase FASTA conversion producing otu97_centroids_UPPERCASE.fasta
6. initial metadata/location construction
7. original demultiplexing + 350-550 bp filtering workflow
8. Greengenes2 V3-V4 reference preparation
9. final R plotting/statistical scripts

These must not be presented as original code if reconstructed later.

---

## 9. Next recovery step

Next test:

For sample_0820CW-3_AG_P01_16S, verify that:

final_clusters.tsv
→ raw-read counts per cluster
→ successful racon_cl_id_<cluster>/consensus.fasta only
→ 52 retained clusters
→ 4,609 supporting reads
→ exact match to the 52 RESCUE rows in pooled provenance.

If exact, repeat algorithm across all 568 baseline + 8 rescue samples.

If the regenerated provenance exactly matches the existing 22,528-row table,
create a clearly-labelled RECONSTRUCTED pooling/provenance script.

Do not modify the historical source directories during reconstruction.
