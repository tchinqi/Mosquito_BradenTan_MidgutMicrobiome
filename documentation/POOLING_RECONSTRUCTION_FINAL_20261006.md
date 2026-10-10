# Final Racon pooling reconstruction — CLOSED
Date: 2026-10-06

The literal original interactive pooling code was not recovered.

However, the historical pooling algorithm was reconstructed from the
preserved NGSpeciesID/Racon sample-level outputs and validated against
all 22,528 final provenance records.

Recovered algorithm:

1. Iterate samples in manifests/all_576_midgut_samples.txt order.
2. ORIGINAL samples use:
   03_ngspeciesid_racon/min2_baseline/samples/<sample_id>/
3. RESCUE samples use:
   03_ngspeciesid_racon/min2_rescue/samples/<sample_id>/
4. Retain clusters containing:
   racon_cl_id_<cluster_id>/consensus.fasta
5. Sort retained cluster IDs lexicographically as strings.
6. supporting_reads equals the number of FASTQ records in:
   reads_to_consensus_<cluster_id>.fastq
7. Assign global sequence IDs sequentially:
   RACON_000001 through RACON_022528
8. Write one unwrapped nucleotide sequence per FASTA record.
9. Write provenance columns in this exact order:
   sequence_id
   sample_id
   cluster_id
   supporting_reads
   source

Final composition:
- 576 samples
- 568 ORIGINAL samples
- 8 RESCUE samples
- 22,181 ORIGINAL consensuses
- 347 RESCUE consensuses
- 22,528 total consensuses

Full validation:
- Missing sample directories: 0
- Missing local Racon consensuses: 0
- Missing polishing FASTQs: 0
- Missing pooled sequence IDs: 0
- Supporting-read mismatches: 0
- Sequence mismatches: 0

RECONSTRUCTED OUTPUT VALIDATION:

pooled_racon_polished_576.fasta
= EXACT BYTE-FOR-BYTE MATCH to historical file

pooled_racon_polished_576_provenance.tsv
= EXACT BYTE-FOR-BYTE MATCH to historical file

Status:
RECONSTRUCTED AND BYTE-IDENTICAL TO HISTORICAL OUTPUTS.

This reconstructed script should remain labelled RECONSTRUCTED rather
than ORIGINAL because the historical interactive source code itself
was not preserved.
