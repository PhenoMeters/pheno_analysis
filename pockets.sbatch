#!/bin/sh
#SBATCH -A proteometer
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH --job-name="pockets_analysis"
#SBATCH -o clean_out.log
#SBATCH -e clean_err.txt
#SBATCH --mail-user=alexandria.im@pnnl.gov
#SBATCH --mail-type END
#SBATCH --cpus-per-task=64


module purge
module load gcc/7.5.0
module load python/miniconda23.3.1
source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh
conda activate ml_ptm

python pockets.py /people/imal967/git_repos/pheno_analysis/phosphosite_for_pockets.csv /people/imal967/git_repos/pheno_analysis/pockets_data.csv