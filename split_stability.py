import pandas as pd
import numpy as np 
import glob
from biopandas.pdb import PandasPdb
import warnings
warnings.filterwarnings("ignore")
from get_distances import *
import sys
from multiprocessing import Pool
import time
import math
import datetime

# read in the file first and keep it as a global variable

stability_df = pd.read_csv("/rcfs/projects/proteometer/ProtVar/predictions/stability/2024.05.28_foldx_energy.csv")

def get_stability_only_uniprot(uniprot_id: str, stability_df = stability_df, output_dir=):
    stability_df[stability_df['protein_acc'] == uniprot_id].copy().to_csv(output_dir)
    return(stability_df[stability_df['protein_acc'] == uniprot_id].copy())

def split_stability_parallel(number_of_threads: int, output_path: str, stability_df = stability_df): 
    unique_uniprot_stability = stability_df["protein_acc"].unique()

