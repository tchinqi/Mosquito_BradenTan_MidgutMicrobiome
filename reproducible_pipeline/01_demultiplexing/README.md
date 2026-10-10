# 01 - Demultiplexing (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `02_run_minibar.pbs` | Demultiplexes the merged Nanopore FASTQ into per-well FASTQ files with Minibar (tag edit distance 2; tags and primers trimmed). |
| `make_demux_braden.py` | Builds the Minibar demultiplexing sheet (9-bp forward/reverse tags + 338F/806R primers per well) and checks tag length, primer sequence and unique tag pairs. |
| `minibar.py` | Minibar 0.25 demultiplexer (third-party tool; unmodified copy, no header added). |
