import os
import glob
from pathlib import Path
import pandas as pd

### Constant SLURM Parameters ###
account = "ptm_annotation"
time = "4-00:00:00"
queue = "slurm"
node = '1'
core = "64"
# mail = "alexandria.im@pnnl.gov"  # will send a email once finished or terminated or failed
mail = "song.feng@pnnl.gov"  # will send a email once finished or terminated or failed
modules = ["gcc/7.5.0", "python/miniconda23.3.1"]
extras = [
        "source /share/apps/python/miniconda23.3.1/etc/profile.d/conda.sh",
        "conda activate /people/imal967/.conda/envs/ml_ptm",
        "export PATH=/people/imal967/.conda/envs/ml_ptm/bin:$PATH",
        ]

### INPUTS ###
interfaces_python_script = "/qfs/projects/proteometer/pheno_analysis/interfaces_all.py"
pockets_python_script = "/qfs/projects/proteometer/pheno_analysis/pockets_all.py"
working_dir = "/qfs/projects/proteometer/pheno_analysis"
stability_chunks_dir = "/qfs/projects/proteometer/pheno_analysis/split_stability_input"
interfaces_output_dir = "/qfs/projects/proteometer/pheno_analysis/interfaces_chunk_output"
pockets_output_dir = "/qfs/projects/proteometer/pheno_analysis/pockets_chunk_output"


all_stability_chunks = glob.glob("/qfs/projects/proteometer/pheno_analysis/split_stability_input/stability_chunk_*.csv")

for stability_chunk in all_stability_chunks:
    filename_only = os.path.basename(stability_chunk)[:-4]
    batch_number = filename_only.split("_")[2]

    # set changing SLURM parameters
    pockets_job_name = batch_number  +  "_multisub_pockets"
    interfaces_job_name = batch_number  +  "_multisub_interfaces"
    p_out = f"{working_dir}/logs/{batch_number}_multisub_pockets.log"
    p_err = f"{working_dir}/logs/{batch_number}_multisub_pockets.err"
    i_out = f"{working_dir}/logs/{batch_number}_multisub_interfaces.log"
    i_err = f"{working_dir}/logs/{batch_number}_multisub_interfaces.err"



    # set filenames
    pockets_output_file = pockets_output_dir + f"/pockets_output_batch_{batch_number}.csv"
    interfaces_output_file = interfaces_output_dir + f"/interfaces_output_batch_{batch_number}.csv"

    # create pockets slurm template
    pockets_slurm_template = f"#!/bin/sh\n#SBATCH -A {account}\n#SBATCH -t {time}\n#SBATCH -N {node}\n#SBATCH --cpus-per-task={core}\n#SBATCH -p {queue}"
    pockets_slurm = [pockets_slurm_template]
    pockets_slurm.append(f"\n#SBATCH -o {p_out}\n#SBATCH -e {p_err}\n#SBATCH --mail-user={mail}\n#SBATCH --mail-type END\n \nmodule purge\n")
    for module in modules:
        pockets_slurm.append(f"module load {module}\n\n")
    for extra in extras:
        pockets_slurm.append(f"{extra}\n")
    pockets_slurm.append(f"cd {working_dir} \n")
    pockets_slurm.append(f"python {pockets_python_script} {stability_chunk} {pockets_output_file} \n")


    # create interfaces slurm template
    interfaces_slurm_template = f"#!/bin/sh\n#SBATCH -A {account}\n#SBATCH -t {time}\n#SBATCH -N {node}\n#SBATCH --cpus-per-task={core}\n#SBATCH -p {queue}"
    interfaces_slurm = [interfaces_slurm_template]
    interfaces_slurm.append(f"\n#SBATCH -o {i_out}\n#SBATCH -e {i_err}\n#SBATCH --mail-user={mail}\n#SBATCH --mail-type END\n \nmodule purge\n")
    for module in modules:
        interfaces_slurm.append(f"module load {module}\n\n")
    for extra in extras:
        interfaces_slurm.append(f"{extra}\n")
    interfaces_slurm.append(f"cd {working_dir} \n")
    interfaces_slurm.append(f"python {interfaces_python_script} {stability_chunk} {interfaces_output_file} \n")
    
    # write the slurm files
    pockets_slurm_filename = f"{working_dir}/sbatch_scripts/batch_pockets_{batch_number}.sbatch"
    interfaces_slurm_filename = f"{working_dir}/sbatch_scripts/batch_interfaces_{batch_number}.sbatch"
    with open(pockets_slurm_filename, 'w') as sbatch:
        sbatch.write(''.join(pockets_slurm))
    with open(interfaces_slurm_filename, 'w') as sbatch:
        sbatch.write(''.join(interfaces_slurm))


    # submit slurm command
    pockets_submission = f"sbatch {pockets_slurm_filename}"
    interfaces_submission = f"sbatch {interfaces_slurm_filename}"
    os.system(pockets_submission)
    os.system(interfaces_submission)
    print(f"submitted job {batch_number}")

    








