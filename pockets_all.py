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
import pickle
import os


'''
this script takes in psp data and pockets data.  It will check if the phosphosites are located in any pockets, which is the nearest, and the distance to each
usage:   pockets.py <path to phosphosite data> <path to pockets data>

'''

pockets_data = pd.read_csv("/rcfs/projects/proteometer/ProtVar/predictions/pockets/2024.05.28_pockets.tsv", delimiter='\t', header=0)

# currently putting these outside the functions because i dont' know how to put this in with multithreading. but this is pretty bad practice
def run_parallel_pockets(number_of_threads, stability_data):

    # get all of the unique uniprots
    unique_uniprot_stability = [stability_data[stability_data["uniprot_id"]==uniprot_id].copy() for uniprot_id in stability_data["uniprot_id"].unique()]


    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(find_pockets_per_uniprot, unique_uniprot_stability)
    finish_time = time.perf_counter()
    
    # save the csv and output the start and end times
    print("Program finished in {} seconds - using multiprocessing with {} cores".format(str(datetime.timedelta(seconds=finish_time-start_time)), number_of_threads))
    concatenated_output = pd.concat(output)
    return(concatenated_output)

'''
this function does pockets calcuations for each uniprot tthat it is given
'''

# for each unique uniprotID...
# for uniprot in unique_uniprots:
def find_pockets_per_uniprot(uniprot_only_stability, pockets_data = pockets_data, pickle_output = "/qfs/projects/proteometer/pheno_analysis/FULL_pockets_pickle_files"):
    #print("start")
    # isolate to psp and pockets in each uniprot
    uniprot = uniprot_only_stability["protein_acc"].to_list()[0]
    pickle_file_path = f"{pickle_output}/{uniprot}.pkl"
    pocket_only_uniprot = pockets_data[pockets_data['uniprot_id'] == uniprot] # isolate pockets data to only that uniprot

    if os.path.isfile(pickle_file_path):
        with open(pickle_file_path, 'rb') as handle:
            psp_data = pickle.load(handle)
        return(psp_data)
    
    # add columns to df
    uniprot_only_stability['closest_pocket'] = ""
    uniprot_only_stability['inside_pocket'] = 0
    uniprot_only_stability['min_distance_from_pocket'] = np.nan
    uniprot_only_stability['mean_distance_from_pocket'] = np.nan

    # parse your structure here
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "-F1-*")
    print("name of pdb is:", pdb_name)
    if pdb_name:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])
        input_struct = ppdb.df['ATOM']


    # for each psp
        for phosphosite_row_index in uniprot_only_stability.index:
            if pd.notna(uniprot_only_stability.loc[phosphosite_row_index,'res_number']):
                residue_num = int(uniprot_only_stability.loc[phosphosite_row_index,'res_number']) # finding the residue number of the psp
                min_dist = np.inf # make min dist extremely high at first
                mean_dist = np.inf
                #print(residue_num)
                # use the residue # to get the coordinates in space from pdb file
                
                pocket_list = ""
                for pocket_index in pocket_only_uniprot.index : # get all the residues in all of the pockets 
                    if pd.notna(pocket_only_uniprot.loc[pocket_index,'pocket_resid']):
                        pocket_residues = pocket_only_uniprot.loc[pocket_index,'pocket_resid']

                        # format the residues
                        pocket_residues = [int(e) for e in pocket_residues[1:-1].split(",")]
                        #print(pocket_residues)
                        if residue_num in pocket_residues:
                            uniprot_only_stability.loc[phosphosite_row_index,'inside_pocket'] = 1 # if residue is in the pocket, put 1 in the inside pocket column
                            pocket_list = ','.join([pocket_list, str(pocket_only_uniprot.loc[pocket_index,'pocket_id'])])
                            uniprot_only_stability.loc[phosphosite_row_index,'closest_pocket'] = pocket_list # put unique pocketID in closest pocket
                            uniprot_only_stability.loc[phosphosite_row_index,'distance_from_pocket'] = 0.0
                            new_min_dist, new_mean_dist = find_min_and_mean_distance(input_struct, residue_num, pocket_residues)
                            if new_mean_dist < mean_dist:
                                mean_dist = new_mean_dist
                                uniprot_only_stability.loc[phosphosite_row_index,'mean_distance_from_pocket'] = mean_dist 
                            min_dist = 0.0
                        elif min_dist != 0.0:                            #print("phosphosite isn't in any pockets")
                            input_struct = ppdb.df['ATOM']
                            new_min_dist, new_mean_dist = find_min_and_mean_distance(input_struct, residue_num, pocket_residues)

                            #print("the new dist is:" , new_dist)
                            if new_mean_dist:
                                if mean_dist > new_mean_dist: # if this is the smallest distance so far, replace min_dist with new_dist
                                    pocket_to_add = (pocket_only_uniprot.loc[pocket_index,'interaction_id'].split('_'))
                                    pocket_to_add.remove(uniprot)
                                    uniprot_only_stability.loc[phosphosite_row_index,'closest_pocket'] = pocket_to_add[0]
                                    uniprot_only_stability.loc[phosphosite_row_index,'mean_distance_from_pocket'] = new_mean_dist # replace distance_from_pocket with min_dist
                                    mean_dist = new_mean_dist
                                    uniprot_only_stability.loc[phosphosite_row_index,'min_distance_from_pocket'] = new_min_dist
                                    min_dist = new_min_dist


        with open(pickle_file_path, 'wb') as handle:
            pickle.dump(uniprot_only_stability, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(uniprot_only_stability)
                    

if __name__ == "__main__":
    num_threads = 64
    stability_data = pd.read_csv("/rcfs/projects/proteometer/human_proteome_precaculated/TEST_stability_precalculated.csv")
    output_location = "/people/imal967/git_repos/pheno_analysis/TEST_pockets_precaculated.csv"
    df_to_export = run_parallel_pockets(num_threads, stability_data)
    pd.DataFrame(df_to_export).to_csv(output_location)
