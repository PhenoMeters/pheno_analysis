import pandas as pd
import numpy as np 
import glob
from biopandas.pdb import PandasPdb
import warnings
warnings.filterwarnings("ignore")
from get_distances import *
import os
from multiprocessing import Pool
import time
import math
import datetime
import pickle


'''
This function adds plddt data to the phosphosite dataset, using the stability data.
'''

def run_parallel_stability(number_of_threads: int, phosphosite_data_PATH: str) -> pd.DataFrame:
    # reading in data
    phosphosite_data = pd.read_csv(phosphosite_data_PATH) 

    # adding columns to phosphosite data
    phosphosite_data['plddt'] = np.nan
    phosphosite_data['foldx_ddg_min'] = np.nan
    phosphosite_data['foldx_ddg_max'] = np.nan
    phosphosite_data['foldx_ddg_abs_max'] = np.nan
    phosphosite_data['foldx_ddg_abs_median'] = np.nan

    unique_uniprot_psp = [phosphosite_data[phosphosite_data["uniprot_id"]==uniprot_id].copy() for uniprot_id in phosphosite_data["uniprot_id"].unique()]

    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(add_plddt_and_foldx_ddg, unique_uniprot_psp)
    finish_time = time.perf_counter()
    full_df = pd.concat(output)

    print("Program finished in {} seconds".format(str(datetime.timedelta(seconds=finish_time-start_time))))
    return(full_df)

def add_plddt_and_foldx_ddg(psp_only_uniprot: pd.DataFrame, pickle_output = "/qfs/projects/proteometer/pheno_analysis/stability_pickle_files") -> pd.DataFrame:
    uniprot = psp_only_uniprot["uniprot_id"].to_list()[0] # getting uniprot string
    stability_uniprot_path = glob.glob("/rcfs/projects/proteometer/ProtVar/predictions/stability/stability_split/*" + uniprot + "*_stability.csv")
    if stability_uniprot_path:
            
        stability_only_uniprot = pd.read_csv(stability_uniprot_path[0])
        pickle_file_path = f"{pickle_output}/{uniprot}_stability.pkl"
        stability_only_uniprot['position'] = stability_only_uniprot['position'].astype(float)
        print(uniprot)

        # check if the pickle file already exists
        if os.path.isfile(pickle_file_path):
            with open(pickle_file_path, 'rb') as handle:
                psp_data = pickle.load(handle)
            return(psp_data)
        

        for phosphosite_row_index in psp_only_uniprot.index:
            psp_res_number = psp_only_uniprot.loc[phosphosite_row_index,'res_number']
            #print(psp_res_number)
            if psp_res_number in stability_only_uniprot['position'].unique():
                stability_only_res_uniprot = stability_only_uniprot[stability_only_uniprot['position'] == psp_res_number] #isolate to only residue number
                plddt_to_add = stability_only_res_uniprot['plddt'].values[0]
                psp_only_uniprot.loc[phosphosite_row_index, 'plddt'] = plddt_to_add

                ### adding fold x data
                psp_only_uniprot.loc[phosphosite_row_index, 'foldx_ddg_min'] = stability_only_res_uniprot['foldx_ddg'].min()
                psp_only_uniprot.loc[phosphosite_row_index, 'foldx_ddg_max'] = stability_only_res_uniprot['foldx_ddg'].max()
                psp_only_uniprot.loc[phosphosite_row_index, 'foldx_ddg_abs_max'] = stability_only_res_uniprot['foldx_ddg'].abs().max()
                psp_only_uniprot.loc[phosphosite_row_index, 'foldx_ddg_abs_median'] = stability_only_res_uniprot['foldx_ddg'].abs().median()

        with open(pickle_file_path, 'wb') as handle:
            pickle.dump(psp_only_uniprot, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(psp_only_uniprot)

if __name__ == "__main__":
    psp = "/people/imal967/git_repos/pheno_analysis/phosphosite_for_pockets.csv"
    output_location = "/people/imal967/git_repos/pheno_analysis/2024.05.28_foldx_energy.csv"
    num_threads = 64
    df_to_export = run_parallel_stability(num_threads, psp)
    pd.DataFrame(df_to_export).to_csv(output_location)