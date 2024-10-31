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
    paths_to_stability = glob.glob("/rcfs/projects/proteometer/ProtVar/predictions/stability/stability_split/*_stability.csv")

    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(add_plddt_and_foldx_ddg, paths_to_stability)
    finish_time = time.perf_counter()
    full_df = pd.concat(output)

    print("Program finished in {} seconds".format(str(datetime.timedelta(seconds=finish_time-start_time))))
    return(full_df)

def add_plddt_and_foldx_ddg(stability_path:str, pickle_output = "/qfs/projects/proteometer/pheno_analysis/FULL_stability_pickle_files") -> pd.DataFrame:
    stability_df = pd.read_csv(stability_path)


    # creating new df and adding columns

    # getting uniprot string
    uniprot = stability_df["protein_acc"].to_list()[0]
    pickle_file_path = f"{pickle_output}/{uniprot}_stability.pkl" # making the name of the pickle file

    # check if the pickle file already exists
    if os.path.isfile(pickle_file_path):
        with open(pickle_file_path, 'rb') as handle:
            stability_data = pickle.load(handle)
        return(stability_data)
    
    stability_df["foldx_ddg_min"] = stability_df["foldx_ddg"]
    stability_df["foldx_ddg_max"] = stability_df["foldx_ddg"]
    stability_df["foldx_ddg_abs_max"] = stability_df["foldx_ddg"]
    stability_df["foldx_ddg_abs_median"] = stability_df["foldx_ddg"]
    agg_method = {
        "protein_acc": lambda x: x.iloc[0],
        "wild_type": lambda x: x.iloc[0],
        "position": lambda x: x.iloc[0],
        "plddt": lambda x: x.iloc[0],
        "foldx_ddg_min": lambda x: x.min(),
        "foldx_ddg_max": lambda x: x.max(),
        "foldx_ddg_abs_max": lambda x: x.abs().max(),
        "foldx_ddg_abs_median": lambda x: x.abs().median(),
    }

    new_df = stability_df.groupby(by="position", as_index=False).agg(agg_method) 

    with open(pickle_file_path, 'wb') as handle:
        pickle.dump(new_df, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(new_df)

if __name__ == "__main__":
    psp = "/people/imal967/git_repos/pheno_analysis/phosphosite_for_pockets.csv"
    output_location = "/rcfs/projects/proteometer/ProtVar/predictions/stability/stability_with_ddg/stability_with_ddg_calculations.csv"
    num_threads = 64
    df_to_export = run_parallel_stability(num_threads, psp)
    pd.DataFrame(df_to_export).to_csv(output_location)