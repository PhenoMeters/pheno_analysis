#!/bin/sh
#SBATCH -A proteometer
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p "slurm"
#SBATCH --job-name="pockets_analysis"
#SBATCH -o pockets_analysis_out.log
#SBATCH -e pockets_analysis_err.txt
#SBATCH --mail-user=alexandria.im@pnnl.gov
#SBATCH --mail-type END


module purge
module load gcc/7.5.0
module load python/miniconda23.3.1
source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh
conda activate ml_ptm
export PATH=/people/imal967/.conda/envs/ml_ptm/bin:$PATH


python pockets.py /people/imal967/git_repos/pheno_analysis/phosphosite_for_pockets.csv /people/imal967/git_repos/pheno_analysis/merged_pockets_FINAL.csv