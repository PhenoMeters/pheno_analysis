#!/bin/sh
#SBATCH -A ptm_annotation
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p "slurm"
#SBATCH --job-name="procaliper_full"
#SBATCH -o logs/procaliper_analysis_out.log
#SBATCH -e logs/procaliper_analysis_err.txt
#SBATCH --mail-user=alexandria.im@pnnl.gov
#SBATCH --mail-type END


module purge
module load gcc/7.5.0
module load python/miniconda23.3.1
source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh
conda activate procaliper
export PATH=/people/imal967/.conda/envs/procaliper/bin:$PATH


python use_procaliper.py /rcfs/projects/proteometer/human_proteome_precaculated/stability_precaculated.csv /people/imal967/git_repos/pheno_analysis/procaliper_output.csv