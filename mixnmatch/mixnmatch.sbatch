#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 1

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate /public4/software/conda_env/mixnmatch

perl /public4/software/mixnmatch/mixnmatch/simulate_admixed_genomes_v6.pl  hybrid.cfg
