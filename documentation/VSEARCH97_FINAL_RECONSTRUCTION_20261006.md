# Final LC/MC VSEARCH97 reconstruction — CLOSED
Date: 2026-10-06

VSEARCH version:
v2.31.0_linux_x86_64

Recovered command parameters:
- --cluster_fast
- --id 0.97
- --iddef 2
- --strand both
- --relabel OTU_
- --threads 1

FINAL terminology:
LC = 0.1% per-sample support threshold
- input consensuses: 18,751
- 97% OTUs: 7,642

MC = max(0.5% of sample reads, 20 reads)
- input consensuses: 7,375
- 97% OTUs: 867

Historical folder warning:
LC_max0p5pct_20reads is the current MC dataset despite its historical name.

Reconstructed outputs were tested against the historical outputs.

LC:
- otu97_centroids.fasta: EXACT BYTE-FOR-BYTE MATCH
- otu97_clusters.uc: EXACT BYTE-FOR-BYTE MATCH

MC:
- otu97_centroids.fasta: EXACT BYTE-FOR-BYTE MATCH
- otu97_clusters.uc: EXACT BYTE-FOR-BYTE MATCH

Status:
RECONSTRUCTED AND BYTE-IDENTICAL TO HISTORICAL OUTPUTS.
