#!/bin/sh
#SBATCH -A proteometer
#SBATCH -t 4-00:00:00
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -p "slurm"
#SBATCH --job-name="pockets_all"
#SBATCH -o pockets_all_out.log
#SBATCH -e pockets_all_err.txt
#SBATCH --mail-user=song.feng@pnnl.gov
#SBATCH --mail-type END


module purge
module load gcc/14.1.0
module load python/miniconda24.4.0
source /share/apps/python/miniconda24.4.0/etc/profile.d/conda.sh
conda activate /qfs/projects/proteometer/proteometer_conda
export PATH=/qfs/projects/proteometer/proteometer_conda/bin:$PATH

cd /qfs/projects/proteometer/pheno_analysis

python pockets_all.py /rcfs/projects/proteometer/ProtVar/predictions/pockets/2024.05.28_pockets.tsv ./results 
