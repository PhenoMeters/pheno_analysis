#!/bin/sh
#SBATCH -A proteometer
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p "slurm" 
#SBATCH --job-name="full_stability_psp"
#SBATCH -o logs/OUT_stability.log
#SBATCH -e logs/ERR_stability.txt
#SBATCH --mail-user=alexandria.im@pnnl.gov
#SBATCH --mail-type END


module purge
module load gcc/7.5.0
module load python/miniconda23.3.1
source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh
conda activate ml_ptm
export PATH=/people/imal967/.conda/envs/ml_ptm/bin:$PATH

python full_stability.py