#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 1

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate /public4/software/conda_env/ancestryinfer
module load java

perl /public4/software/ancestryinfer/Ancestry_HMM_parallel_v7.pl hmm.cfg
