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

HF_HOME=$W/hf HF_HUB_OFFLINE=1 $ENVS/esmfold/bin/python $P/run_esmfold.py --run-dir $RUN --wt-reps 5 --device cuda
$ENVS/colabfold/bin/python $P/run_colabfold.py --run-dir $RUN --step fold \
   --colabfold-batch $ENVS/colabfold/bin/colabfold_batch --data $W/colabfold \
   --wt-seeds 5 --mut-models 5
echo FOLD_DONE
