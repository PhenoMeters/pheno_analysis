#!/bin/sh
#SBATCH -A ptm_annotation
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p "slurm"
#SBATCH --job-name="merge_results"
#SBATCH -o logs/merge_results_out.log
#SBATCH -e logs/merge_results_err.txt
#SBATCH --mail-user=alexandria.im@pnnl.gov
#SBATCH --mail-type END


module purge
module load gcc/7.5.0
module load python/miniconda23.3.1
source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh
conda activate ml_ptm
export PATH=/people/imal967/.conda/envs/ml_ptm/bin:$PATH

python merge_results.py