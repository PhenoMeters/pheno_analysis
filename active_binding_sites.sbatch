#!/bin/sh
#SBATCH -A proteometer
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p "slurm"
#SBATCH --job-name="binding_sites_test"
#SBATCH -o logs/LOG_binding_sites_test.log
#SBATCH -e logs/ERR_binding_sites_test.txt
#SBATCH --mail-user=alexandria.im@pnnl.gov
#SBATCH --mail-type END


module purge
module load gcc/7.5.0
module load python/miniconda23.3.1
source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh
conda activate procaliper
export PATH=/people/imal967/.conda/envs/procaliper/bin:$PATH

python active_binding_sites.py 