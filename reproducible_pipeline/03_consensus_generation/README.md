# 03 - Within-sample consensus: NGSpeciesID + Racon (Atlas9)

Run order follows the script numbers. Paths inside the scripts point to the Atlas9 project folder `/hpctmp/<user>/Braden_16S_20260819/`.

| Script | What it does |
|---|---|
| `12_min2_racon_baseline_midgut.pbs` | Runs NGSpeciesID per well (--ont --consensus, abundance ratio 2.5/N so clusters need >= 2 reads, Racon x3) using the patched consensus module that skips individual Racon failures. |
| `14_rescue_8_min2_racon.pbs` | Re-runs the 8 wells whose baseline run hit Racon failures, with the same settings. |
| `RECONSTRUCTED_build_final_racon_pool.py` | RECONSTRUCTED: pools the 22,528 polished consensuses with sample/cluster/read provenance; validated byte-for-byte against the historical output. |
