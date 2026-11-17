# Hoffman2 full run

**The cluster moved to Slurm on 2026-11-05.** About 94% of nodes left the Altair Grid
Engine partition; the ~30 that remain are mostly draining. A GE GPU job submitted after
that date does not dispatch — it sits in `qw` against queues marked unreachable/disabled.
Use `fold_slurm.sh` (`sbatch`), not `fold_gpu.sh` (`qsub`).

| | legacy (Grid Engine) | production (Slurm) |
|---|---|---|
| login host | `hoffman2.idre.ucla.edu` | **`hoffman2.oarc.ucla.edu`** |
| OS | CentOS 7 (glibc 2.17) | Rocky 9.8 (glibc 2.34) |
| submit | `qsub` | `sbatch` |
| script | `fold_gpu.sh` (kept for reference) | `fold_slurm.sh` |

`$HOME`, `$SCRATCH` and project dirs are shared between both partitions, so staged data
needs no re-upload. The conda envs built under CentOS 7 run fine on Rocky 9 (verified:
`jax 0.4.23 / jaxlib 0.4.23 / haiku 0.0.10` imports clean). `cuda/12.3` does **not**
exist on Rocky 9 — available are 11.7, 11.8, 12.2, 12.4, 12.8, 13.0, 13.2; use `cuda/12.2`
for jaxlib 0.4.23.

**Do not pass `-p/--partition` or `--qos`.** Slurm routes by account, and because
`pi_kappel` owns nodes, `-A pi_kappel` is sent to the shared partition `pi_shared` with
QOS `pi_shared24` automatically. Asking for `-p campus` is rejected outright:
*"account pi_kappel has its own nodes in SLURM, so its jobs run on pi_shared"*.
`--qos=highp` would move the job onto the lab's own `pi_kappel` nodes and bill the Kappel
Lab allocation — this is CAIR work, so **never `highp`**.

`pi_shared24` allows 24 h wall time and **4 concurrent GPUs per user**, which is why the
array is `0-3`. `pi_shared` has ~25 GPU nodes (1-4 GPUs each) and usually a few idle.

| where | what |
|---|---|
| `$HOME/cair_fold/envs/{colabfold,esmfold}` | conda envs (~10 GB of the 57 GB home quota) |
| `$HOME/cair_fold/pipeline/` | copy of `../pipeline/*.py` |
| `$HOME/cair_fold/logs/` | job logs — `sbatch` needs this to exist before submitting |
| `$SCRATCH/cair_fold/weights/` | AF2 params + ESMFold HF cache (scratch purges after 14 days idle; `setup_envs.sh` re-downloads) |
| `$SCRATCH/cair_fold/<run>/` | run dir: sequences.csv + colabfold/msa_wt built on the laptop |

1. Laptop: `python pipeline/build_subset.py --out RUN --per-stratum 25` and
   `python pipeline/run_colabfold.py --run-dir RUN --step msa` (MSA server needs internet;
   compute nodes may not have it).
2. `rsync -a RUN/ ijeon@hoffman2.oarc.ucla.edu:$SCRATCH/cair_fold/<run>/` and
   `rsync -a pipeline/ ijeon@hoffman2.oarc.ucla.edu:cair_fold/pipeline/`
3. `cd $HOME/cair_fold && mkdir -p logs`, then
   `sbatch --export=ALL,RUN=$SCRATCH/cair_fold/full_v1 hoffman2/fold_slurm.sh`
   (a 4-task array; raise `--array` to split further — `run_colabfold.py --shard k/n`
   splits by sequence and the a3m writes are atomic, so tasks may run concurrently).
4. Pull back: `rsync -a ijeon@hoffman2.oarc.ucla.edu:$SCRATCH/cair_fold/<run>/preds/ RUN/preds/`,
   then locally `python pipeline/summarize.py --run-dir RUN` and
   `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -I pipeline/noise_floor.py --run-dir RUN`.

Full run settings: AF2 all 5 models for every sequence, 5 seeds for WT, `--save-recycles`.
ESMFold is **not** run on the cluster — its folds are all laptop CPU/bf16 and must stay
on one device and precision, or the replicate noise floor stops being a property of the
predictor and becomes a property of the hardware. Its recycle sweep runs locally too
(`run_esmfold.py --recycle-sweep 0,1,2,3,4 --sweep-only`).
