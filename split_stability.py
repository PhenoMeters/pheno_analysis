import pandas as pd
import numpy as np 
from biopandas.pdb import PandasPdb
import warnings
warnings.filterwarnings("ignore")
from get_distances import *
from multiprocessing import Pool
import time
import datetime
import pickle

# global variables
stability_df = pd.read_csv("/rcfs/projects/proteometer/ProtVar/predictions/stability/2024.05.28_foldx_energy.csv")
output_location = "/rcfs/projects/proteometer/ProtVar/predictions/stability/stability_split"

def get_stability_only_uniprot(uniprot_id: str, output_path = output_location, stability_df = stability_df):
    name_of_output_file = "/" + uniprot_id + "_stability.csv"
    stability_df[stability_df['protein_acc'] == uniprot_id].copy().to_csv(output_path + name_of_output_file)
    return(stability_df[stability_df['protein_acc'] == uniprot_id].copy())

def split_stability_parallel(number_of_threads: int, stability_df = stability_df): 
    unique_uniprots_stability = stability_df["protein_acc"].unique()
    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(get_stability_only_uniprot, unique_uniprots_stability)
    finish_time = time.perf_counter()
    
    print("Program finished in {} seconds".format(str(datetime.timedelta(seconds=finish_time-start_time))))
    return(output)

if __name__ == "__main__":
    num_threads = 64
    export_pickle = "/rcfs/projects/proteometer/ProtVar/predictions/stability/stability_chunks.pkl"
    data_to_export = split_stability_parallel(num_threads)
    with open(export_pickle, 'wb') as handle:
        pickle.dump(data_to_export, handle, protocol=pickle.HIGHEST_PROTOCOL)
