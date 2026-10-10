# 02 - Read-length filtering and read QC (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `01_build_demux_length_qc.sh` | Counts reads per well after demultiplexing and after 350-550 bp filtering (Minibar 'unk' and 'Multiple_Matches' outputs excluded). |
| `02_summarise_read_depths.sh` | Adds per-well retention (% of reads kept by the length filter) and writes a read-depth summary. |
| `03_export_read_counts_for_supervisor.sh` | Exports the per-well read-count table as CSV (used to set per-sample abundance-ratio cut-offs). |
| `07_length_filter.pbs` | Keeps reads of 350-550 bp (expected V3-V4 amplicon) with SeqKit (seq -g -m 350 -M 550) and writes per-file read statistics. |
