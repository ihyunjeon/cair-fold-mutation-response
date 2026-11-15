# Hoffman2 full run

Shared queue only. Never pass `-l highp` (Kappel Lab allocation; this is CAIR work) and
never pass `-q` on a GPU job.

| where | what |
|---|---|
| `$HOME/cair_fold/envs/{colabfold,esmfold}` | conda envs (~10 GB of the 57 GB home quota) |
| `$HOME/cair_fold/pipeline/` | copy of `../pipeline/*.py` |
| `$SCRATCH/cair_fold/weights/` | AF2 params + ESMFold HF cache (scratch is purged after 14 days idle; `setup_envs.sh` re-downloads) |
| `$SCRATCH/cair_fold/<run>/` | run dir: sequences.csv + colabfold/msa_wt built on the laptop |

1. Laptop: `python pipeline/build_subset.py --out RUN --per-stratum 25` and
   `python pipeline/run_colabfold.py --run-dir RUN --step msa` (MSA server needs internet;
   compute nodes may not have it).
2. `rsync -a RUN/ hoffman2:$SCRATCH/cair_fold/<run>/` and `rsync -a pipeline/ hoffman2:cair_fold/pipeline/`
3. `qsub hoffman2/setup_envs.sh`, then `qsub -hold_jid <setup id> -v RUN=... hoffman2/fold_gpu.sh`
4. Pull back: `rsync -a hoffman2:$SCRATCH/cair_fold/<run>/preds/ RUN/preds/`, then
   `python pipeline/summarize.py --run-dir RUN` locally.

Full run settings: ESMFold 5 WT replicates; AF2 all 5 models for every sequence, 5 seeds for WT.
