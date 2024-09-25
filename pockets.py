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

# currently putting these outside the functions because i dont' know how to put this in with multithreading. but this is pretty bad practice
def run_parallel_pockets(number_of_threads, psp_data):
    #adding columns to psp data

    psp_data['closest_pocket'] = []
    psp_data['inside_pocket'] = 0
    psp_data['distance_from_pocket'] = np.nan

    # get all of the unique uniprots
    unique_uniprot_psp = [psp_data[psp_data["uniprot_id"]==uniprot_id].copy() for uniprot_id in psp_data["uniprot_id"].unique()]


    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(find_pockets_per_uniprot, unique_uniprot_psp)
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
def find_pockets_per_uniprot(psp_only_uniprot, pickle_output = "/qfs/projects/proteometer/pheno_analysis/pocket_pickle_files"):
    #print("start")
    # isolate to psp and pockets in each uniprot
    pockets_data = pd.read_csv("/people/imal967/git_repos/pheno_analysis/pockets_data.csv")
    uniprot = psp_only_uniprot["uniprot_id"].to_list()[0]
    pickle_file_path = f"{pickle_output}/{uniprot}.pkl"
    if os.path.isfile(pickle_file_path):
        with open(pickle_file_path, 'rb') as handle:
            psp_data = pickle.load(handle)
        return(psp_data)
    
    pocket_only_uniprot = pockets_data[pockets_data['uniprot_id'] == uniprot]


    # parse your structure here
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "-F1-*")
    print("name of pdb is:", pdb_name)
    if pdb_name:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])


    # for each psp
        for phosphosite_row_index in psp_only_uniprot.index:
            if pd.notna(psp_only_uniprot.loc[phosphosite_row_index,'res_number']):
                #print(psp_only_uniprot)
                #print(phosphosite_row_index)
                residue_num = int(psp_only_uniprot.loc[phosphosite_row_index,'res_number']) # finding the residue number of the psp
                min_dist = np.inf # make min dist extremely high at first
                #print(residue_num)
                # use the residue # to get the coordinates in space from pdb file
                
                pocket_list = []
                for pocket_index in pocket_only_uniprot.index : # get all the residues in all of the pockets 
                    if pd.notna(pocket_only_uniprot.loc[pocket_index,'pocket_resid']):
                        pocket_residues = pocket_only_uniprot.loc[pocket_index,'pocket_resid']

                        # format the residues
                        pocket_residues = [int(e) for e in pocket_residues[1:-1].split(",")]
                        #print(pocket_residues)
                        if residue_num in pocket_residues:
                            psp_only_uniprot.loc[phosphosite_row_index,'inside_pocket'] = 1 # if residue is in the pocket, put 1 in the inside pocket column
                            psp_only_uniprot.loc[phosphosite_row_index,'closest_pocket'] = pocket_list.append(pocket_only_uniprot.loc[pocket_index,'full_id']) # put unique pocketID in closest pocket
                            psp_only_uniprot.loc[phosphosite_row_index,'distance_from_pocket'] = 0 
                            min_dist = 0.0
                        elif min_dist != 0.0: # if the phosphosite isn't in any pockets
                            #print("phosphosite isn't in any pockets")
                            input_struct = ppdb.df['ATOM']
                            #print(input_struct)
                            print(residue_num)
                            print(pocket_residues)
                            new_dist = find_mean_distances(input_struct, residue_num, pocket_residues)
                            # print(new_dist)
                            if new_dist:
                                if min_dist > new_dist: # if this is the smallest distance so far, replace min_dist with new_dist
                                    psp_only_uniprot.loc[phosphosite_row_index,'closest_pocket'] = [pocket_only_uniprot.loc[pocket_index,'full_id']] # put unique pocketID in closest pocket
                                    psp_only_uniprot.loc[phosphosite_row_index,'distance_from_pocket'] = new_dist # replace distance_from_pocket with min_dist
                                    min_dist = new_dist 
                                    print("added smallest distance:", min_dist)


    #print(uniprot)
    with open(pickle_file_path, 'wb') as handle:
        pickle.dump(psp_only_uniprot, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(psp_only_uniprot)
                    

if __name__ == "__main__":
    num_threads = 64
    psp_data = pd.read_csv(sys.argv[1])
    output_location = sys.argv[2]
    df_to_export = run_parallel_pockets(num_threads, psp_data)
    pd.DataFrame(df_to_export).to_csv(output_location)