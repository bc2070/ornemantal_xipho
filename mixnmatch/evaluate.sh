#!/usr/bin/bash
#SBATCH -p blade,gpu,himem,hugemem
#SBATCH -c 1

module load R
source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate /public4/software/conda_env/mixnmatch


scriptdir=/public4/software/mixnmatch/mixnmatch/
perl $scriptdir/post_hmm_accuracy_shell.pl ancestry-probs-par1_allchrs.tsv ancestry-probs-par2_allchrs.tsv simulated_hybrids_reads_gen50_prop_par1_0.1 0.8 $scriptdir
