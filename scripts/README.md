# Scripts

All executable scripts for the tutorial. **All comments and messages are in English**
for international collaboration.

## Script list

| Script | Purpose | Chapter | Environment |
|---|---|---|---|
| `setup_singularity_wrappers.sh` | Create wrappers for tools inside PGGB SIF | 5 | Host |
| `download_plantgarden.sh` | Download 6 assemblies from Plant GARDEN | 3 | Host |
| `qc01_stats.sh` | Basic FASTA stats (seqkit) | 4 | Anywhere |
| `qc02_visualize.py` | Plot QC metrics (4-panel overview) | 4 | Anywhere |
| `pg01_prepare_input.sh` | PanSN renaming + chromosome grouping | 5 | Anywhere |
| `pg02_run_pggb.sh` | Run PGGB per chromosome | 5 | HPC |
| `pg03_qc_graph.sh` | 4-level graph QC | 6 | HPC |
| `pg04_visualize.sh` | Visualize graph (odgi viz + draw) | 7 | HPC |

## Order of execution

```
1. setup_singularity_wrappers.sh  (one-time setup)
2. download_plantgarden.sh         (data acquisition, Chapter 3)
3. qc01_stats.sh + qc02_visualize.py  (data QC, Chapter 4)
4. pg01_prepare_input.sh           (prepare input, Chapter 5)
5. pg02_run_pggb.sh                (run PGGB, Chapter 5)
6. pg03_qc_graph.sh                (graph QC, Chapter 6)
7. pg04_visualize.sh               (visualize, Chapter 7)
```

## Design notes

- Scripts use `set -uo pipefail` (not `-e`) to allow graceful continuation
  after non-fatal errors, so users can see what worked and what did not.
- Preflight checks are performed at the start of long-running scripts
  (pg02, pg03) to catch missing dependencies early.
- PGGB output filenames embed parameter hashes; scripts use glob patterns
  like `chr09.fa.gz.*.smooth.final.og` to locate them.
- The `-B` parameter of PGGB must be explicitly specified (see Chapter 5).
