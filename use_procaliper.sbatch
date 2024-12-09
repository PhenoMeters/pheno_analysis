#!/bin/sh
#SBATCH -A ptm_annotation
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p "slurm"
#SBATCH --job-name="procaliper_local"
#SBATCH -o logs/local_procaliper_analysis_out.log
#SBATCH -e logs/local_procaliper_analysis_err.txt
#SBATCH --mail-user=alexandria.im@pnnl.gov
#SBATCH --mail-type END


module purge
module load gcc/7.5.0
module load python/miniconda23.3.1
source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh
conda activate procaliper
export PATH=/people/imal967/.conda/envs/procaliper/bin:$PATH


python use_procaliper.py /rcfs/projects/proteometer/UniProt_proteomes/UP000005640_9606_Human/uniprotkb_proteome_UP000005640_AND_revi_2024_11_19.tsv /people/imal967/git_repos/pheno_analysis/procaliper_output.csv