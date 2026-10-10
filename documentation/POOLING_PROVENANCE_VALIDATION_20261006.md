# Final pooling/provenance reconstruction validation
Date: 2026-10-06

## Final dataset

Final pooled FASTA:
04_global_vsearch97/min2_complete576/pooled_racon_polished_576.fasta

Final provenance:
04_global_vsearch97/min2_complete576/pooled_racon_polished_576_provenance.tsv

Validated:
- 22,528 provenance rows
- 22,528 pooled FASTA sequences
- 576 unique samples
- 22,181 ORIGINAL consensuses
- 347 RESCUE consensuses
- 568 ORIGINAL samples
- 8 RESCUE samples

## Recovered source logic

ORIGINAL records originate from:
03_ngspeciesid_racon/min2_baseline/samples/<sample_id>/

RESCUE records originate from:
03_ngspeciesid_racon/min2_rescue/samples/<sample_id>/

For each retained cluster:

1. The polished nucleotide sequence is taken from:
   racon_cl_id_<cluster_id>/consensus.fasta

2. supporting_reads equals the number of FASTQ records in:
   reads_to_consensus_<cluster_id>.fastq

3. Only clusters with a successful final Racon consensus are retained.

4. Failed Racon clusters are excluded.

## Full 22,528-record validation

Results:
- Missing sample directories: 0
- Missing local Racon consensuses: 0
- Missing polishing FASTQs: 0
- Missing pooled sequence IDs: 0
- Supporting-read mismatches: 0
- Sequence mismatches: 0

Conclusion:

ALL 22,528 FINAL PROVENANCE RECORDS VALIDATE EXACTLY AGAINST
THEIR ORIGINAL NGSpeciesID/RACON SAMPLE-LEVEL OUTPUTS.

The biological/data logic of the historical pooling operation has
therefore been recovered exactly.

The literal original interactive pooling code was not recovered.
Any replacement script must therefore be labelled RECONSTRUCTED,
not ORIGINAL.

Next step:
Create a reconstructed pooling/provenance script in a separate
reproducibility directory and test whether it regenerates the
historical final files byte-for-byte or differs only in record order /
global RACON sequence-ID assignment.
