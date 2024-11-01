#!/bin/sh
#SBATCH -A proteometer
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p "slurm"
#SBATCH --job-name="procaliper_test"
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


python use_procaliper.py /people/imal967/git_repos/pheno_analysis/phosphosite_for_pockets.csv /people/imal967/git_repos/pheno_analysis/procaliper_TEST.csv