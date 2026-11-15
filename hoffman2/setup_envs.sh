#!/bin/bash
# One-time env build on a Hoffman2 COMPUTE node (never the login node).
#   qsub hoffman2/setup_envs.sh
# Shared queue only: no -l highp (that allocation is the Kappel Lab's; this is CAIR work).
#$ -cwd
#$ -o $HOME/cair_fold/logs/setup.$JOB_ID.log
#$ -j y
#$ -l h_rt=3:00:00,h_data=8G
#$ -pe shared 2
set -euo pipefail
. /u/local/Modules/default/init/modules.sh
module load miniforge
ENVS=$HOME/cair_fold/envs
W=$SCRATCH/cair_fold/weights
mkdir -p $ENVS $W/colabfold/params $W/hf
export PIP_CACHE_DIR=$SCRATCH/pip_cache
# Hoffman2 is CentOS 7 (glibc 2.17). Current scipy/torch wheels need glibc 2.28, so pip
# silently falls back to source builds that fail. Pin versions that ship manylinux2014
# wheels and forbid source builds so a bad pin fails fast instead.
PIN="numpy<2 scipy<1.12 pandas<2.3 matplotlib<3.9"

# ---- ColabFold 1.5.5 / AF2 (JAX CUDA 12). jax 0.4.23 + dm-haiku 0.0.10 are the versions
# 1.5.5 works with; newer jax removes jax.linear_util and breaks haiku.
[ -x $ENVS/colabfold/bin/colabfold_batch ] || {
  rm -rf $ENVS/colabfold; conda create -y -p $ENVS/colabfold python=3.11
  $ENVS/colabfold/bin/pip install --only-binary=:all: $PIN "colabfold[alphafold]==1.5.5" \
      "jax[cuda12_pip]==0.4.23" "dm-haiku==0.0.10" \
      -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
}
[ -f $W/colabfold/params/download_finished.txt ] || {
  curl -sSL https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar | tar -x -C $W/colabfold/params
  touch $W/colabfold/params/download_finished.txt
}

# ---- ESMFold via HF transformers (no openfold build needed)
[ -x $ENVS/esmfold/bin/python ] || {
  rm -rf $ENVS/esmfold; conda create -y -p $ENVS/esmfold python=3.11
  $ENVS/esmfold/bin/pip install --only-binary=:all: "torch==2.5.1" --index-url https://download.pytorch.org/whl/cu121
  $ENVS/esmfold/bin/pip install --only-binary=:all: $PIN "transformers==4.46.3" accelerate
}
HF_HOME=$W/hf $ENVS/esmfold/bin/python -c "from huggingface_hub import snapshot_download as s; print(s('facebook/esmfold_v1'))"
echo SETUP_DONE
