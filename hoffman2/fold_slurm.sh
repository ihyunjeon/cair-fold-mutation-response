#!/bin/bash
# AF2 arm on the Hoffman2 SLURM partition (Rocky 9).
#
# Why this file exists: on 2026-11-05 Hoffman2 moved ~94% of its nodes off Altair Grid
# Engine onto Slurm. The GE partition keeps only ~30 nodes and the rest are draining, so
# hoffman2/fold_gpu.sh (qsub) no longer dispatches -- a GE GPU job submitted after that
# date waits forever against queues marked unreachable/disabled. $HOME and $SCRATCH are
# shared between the two partitions, so the staged run dir, MSAs and weights are already
# in place; only the scheduler changed.
#
#   ssh ijeon@hoffman2.oarc.ucla.edu          # NOT hoffman2.idre.ucla.edu (legacy GE)
#   cd $HOME/cair_fold && mkdir -p logs
#   sbatch --export=ALL,RUN=$SCRATCH/cair_fold/full_v1 hoffman2/fold_slurm.sh
#
# No -p/--partition and no --qos. pi_kappel owns nodes, so Slurm routes the job by its
# account: -A pi_kappel lands on the shared partition pi_shared (QOS pi_shared24, 24 h,
# up to 4 concurrent GPUs -- hence --array=0-3). Asking for -p campus is rejected
# outright for this account. --qos=highp would move it onto the lab's own pi_kappel
# nodes and bill the Kappel allocation; this is CAIR work, so never highp.
#
#SBATCH --job-name=cair_af2
#SBATCH --account=pi_kappel
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/cair_af2.%A_%a.log
#SBATCH --open-mode=append
set -euo pipefail

: "${RUN:?set RUN, e.g. --export=ALL,RUN=\$SCRATCH/cair_fold/full_v1}"

. /u/local/Modules/default/init/modules.sh
module load cuda/12.2          # cuda/12.3 does not exist on Rocky 9; jaxlib 0.4.23 is a CUDA 12 build

ENVS=$HOME/cair_fold/envs
W=$SCRATCH/cair_fold/weights
P=$HOME/cair_fold/pipeline

SHARD_ID=${SLURM_ARRAY_TASK_ID:-0}
SHARD_N=${SLURM_ARRAY_TASK_COUNT:-1}

echo "host=$(hostname)  job=${SLURM_JOB_ID:-na}  shard=${SHARD_ID}/${SHARD_N}  run=$RUN"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
$ENVS/colabfold/bin/python -c "import jax; d=jax.devices(); print(d); assert d[0].platform=='gpu', 'jax sees no GPU'"

# The ESMFold arm is deliberately NOT run here; see the note in fold_gpu.sh. Its folds are
# all laptop CPU/bf16 and must stay homogeneous or the replicate noise floor becomes a
# measurement of device and precision instead of the predictor.

# AF2. --wt-seeds 5 x --wt-models 5 gives, per construct, 10 within-model seed pairs per
# model (the seed floor) plus the across-model spread. --save-recycles adds the
# recycle-convergence floor, the only floor definition ESMFold also has.
$ENVS/colabfold/bin/python $P/run_colabfold.py --run-dir "$RUN" --step fold \
   --colabfold-batch $ENVS/colabfold/bin/colabfold_batch --data $W/colabfold \
   --wt-seeds 5 --mut-models 5 --save-recycles \
   --shard "${SHARD_ID}/${SHARD_N}"
echo FOLD_DONE
