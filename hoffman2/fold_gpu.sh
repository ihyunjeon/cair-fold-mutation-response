#!/bin/bash
# GPU fold job: ESMFold then AF2 on one run directory.
#   qsub -v RUN=$SCRATCH/cair_fold/full_v1 hoffman2/fold_gpu.sh
# GPU is requested through -l (gpu + cuda consumable). No -q, no highp.
#$ -cwd
#$ -o $HOME/cair_fold/logs/fold.$JOB_ID.log
#$ -j y
#$ -l gpu,cuda=1,h_rt=8:00:00,h_data=32G
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module load cuda/12.3
ENVS=$HOME/cair_fold/envs
W=$SCRATCH/cair_fold/weights
P=$HOME/cair_fold/pipeline
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv
$ENVS/esmfold/bin/python -c "import torch; assert torch.cuda.is_available(), 'torch sees no GPU'; print(torch.cuda.get_device_name())"
$ENVS/colabfold/bin/python -c "import jax; d=jax.devices(); print(d); assert d[0].platform=='gpu', 'jax sees no GPU'"

# The ESMFold arm is deliberately NOT run here. full_v1's ESMFold folds were all produced
# on laptop CPU at bf16, and run_esmfold.py skips any rep whose PDB already exists -- so
# this job would have folded only the MISSING reps, on CUDA at fp16. That puts a
# device+precision difference inside the replicate noise floor, which is the one number
# the floor has to be free of. ESMFold stays homogeneous CPU/bf16, including its
# --recycle-sweep. Uncomment only for a run whose ESMFold arm is folded here start to
# finish, and then fold every rep here, not just the gaps.
# HF_HOME=$W/hf HF_HUB_OFFLINE=1 $ENVS/esmfold/bin/python $P/run_esmfold.py \
#    --run-dir $RUN --wt-reps 5 --device cuda --lm-dtype fp16

# AF2 arm. --wt-seeds 5 at --wt-models 5 gives, per construct, 10 within-model seed pairs
# per model (the seed floor) plus the across-model spread. --save-recycles adds the
# recycle-convergence floor, the only floor definition that also exists for ESMFold.
$ENVS/colabfold/bin/python $P/run_colabfold.py --run-dir $RUN --step fold \
   --colabfold-batch $ENVS/colabfold/bin/colabfold_batch --data $W/colabfold \
   --wt-seeds 5 --mut-models 5 --save-recycles
echo FOLD_DONE
