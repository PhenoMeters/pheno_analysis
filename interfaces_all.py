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


interfaces_data = pd.read_csv("/rcfs/projects/proteometer/ProtVar/predictions/interfaces/2024.05.28_interface_summary_5A.tsv", delimiter='\t', header=0)

def run_parallel_interfaces(number_of_threads, stability_data):

    unique_uniprot_stability = [stability_data[stability_data["protein_acc"]==uniprot_id].copy() for uniprot_id in stability_data["protein_acc"].unique()]


    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(find_interfaces_per_uniprot, unique_uniprot_stability)
    finish_time = time.perf_counter()
    
    # save the csv and output the start and end times
    print("Program finished in {} - using multiprocessing with {} cores".format(str(datetime.timedelta(seconds=finish_time-start_time)), number_of_threads))
    full_df = pd.concat(output)
    return(full_df)




'''
this function does interfaces calcuations for each uniprot tthat it is given
'''

# for each unique uniprotID...
# for uniprot in unique_uniprots:
def find_interfaces_per_uniprot(uniprot_only_stability, interfaces_data = interfaces_data, pickle_output = "/qfs/projects/proteometer/pheno_analysis/FULL_interfaces_pickle_files"):

    # read in interfaces data, uniprot human data, and get uniprot
    uniprot = uniprot_only_stability["protein_acc"].to_list()[0]
    pickle_file_path = f"{pickle_output}/{uniprot}.pkl"
    interface_only_uniprot = interfaces_data.loc[(interfaces_data['uniprot_id1'] == uniprot) | (interfaces_data['uniprot_id2'] == uniprot)] # isolate to uniprot in either 1 or 2

    if os.path.isfile(pickle_file_path):
        with open(pickle_file_path, 'rb') as handle:
            full_data = pickle.load(handle)
        return(full_data)

    # add columns to df
    uniprot_only_stability['closest_interface'] = ""
    uniprot_only_stability['inside_interface'] = 0
    uniprot_only_stability['min_distance_from_interface'] = np.nan
    uniprot_only_stability['mean_distance_from_interface'] = np.nan

    # parse your structure here
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "-F1-*") # change this to have F1 using pdb file name
    print(pdb_name)
    #print("name of pdb is:", pdb_name)
    if len(pdb_name) != 0:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])
        input_struct = ppdb.df['ATOM']


    # for each psp
        for phosphosite_row_index in uniprot_only_stability.index:
            if pd.notna(uniprot_only_stability.loc[phosphosite_row_index,'position']): # only if there is a residue # (this threw an error previously)
                residue_num = int(uniprot_only_stability.loc[phosphosite_row_index,'position']) # finding the residue number of the psp
                min_dist = np.inf # make min dist extremely high at first
                mean_dist = np.inf
                #print(residue_num)
                # use the residue # to get the coordinates in space from pdb file
                
                interface_list = ""
                for interface_index in interface_only_uniprot.index : # get all the residues in all of the interfaces 
                    if pd.notna(interface_only_uniprot.loc[interface_index,'ifresid1']) & pd.notna(interface_only_uniprot.loc[interface_index,'ifresid2']):
                        if interfaces_data.loc[interface_index,'uniprot_id1'] == uniprot:
                            interface_residues = interface_only_uniprot.loc[interface_index,'ifresid1']
                        elif interfaces_data.loc[interface_index,'uniprot_id2'] == uniprot:
                            interface_residues = interface_only_uniprot.loc[interface_index,'ifresid2']
                        
                        # check if it's inside of a interface
                        #print(interfaces_data.loc[interface_index,'interaction_id']) 
                        interface_residues = [int(e[1:]) for e in interface_residues.split(",")] # remove the first letter from each bc it includes residue type ormat the interface_residues because it's a string
                        #print(interface_residues)
                        if residue_num in interface_residues:
                            #print("found inside interface!")
                            uniprot_only_stability.loc[phosphosite_row_index,'inside_interface'] = 1 # if residue is in the interface, put 1 in the inside interface column
                            interface_to_add = (interface_only_uniprot.loc[interface_index,'interaction_id'].split('_'))
                            interface_to_add.remove(uniprot)
                            interface_list = ','.join([interface_list, interface_to_add[0]])
                            uniprot_only_stability.loc[phosphosite_row_index,'closest_interface'] = interface_list # put unique interfaceID in closest interface
                            uniprot_only_stability.loc[phosphosite_row_index,'min_distance_from_interface'] = 0.0 
                            new_mean_dist = find_mean_distance(input_struct, residue_num, interface_residues)
                            if new_mean_dist < mean_dist:
                                mean_dist = new_mean_dist
                                uniprot_only_stability.loc[phosphosite_row_index,'mean_distance_from_interface'] = mean_dist 
                            min_dist = 0.0
                        else: # if the phosphosite isn't in any interfaces
                            new_mean_dist = find_mean_distance(input_struct, residue_num, interface_residues)

                            #print("the new dist is:" , new_dist)
                            if new_mean_dist:
                                if mean_dist > new_mean_dist: # if this is the smallest distance so far, replace min_dist with new_dist
                                    interface_to_add = (interface_only_uniprot.loc[interface_index,'interaction_id'].split('_'))
                                    interface_to_add.remove(uniprot)
                                    uniprot_only_stability.loc[phosphosite_row_index,'closest_interface'] = interface_to_add[0]
                                    uniprot_only_stability.loc[phosphosite_row_index,'mean_distance_from_interface'] = new_mean_dist # replace distance_from_interface with min_dist
                                    mean_dist = new_mean_dist
                                    new_min_dist = find_min_distance(input_struct, residue_num, interface_residues)
                                    uniprot_only_stability.loc[phosphosite_row_index,'min_distance_from_interface'] = new_min_dist
                                    min_dist = new_min_dist
                                    #print("replaced old dist with", min_dist)
        with open(pickle_file_path, 'wb') as handle:
            pickle.dump(uniprot_only_stability, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(uniprot_only_stability) 



                
if __name__ == "__main__":
    num_threads = 64
    stability_data = pd.read_csv("/rcfs/projects/proteometer/human_proteome_precaculated/TEST_stability_precalculated.csv")
    output_location = "/people/imal967/git_repos/pheno_analysis/TEST_interfaces_precaculated.csv"
    df_to_export = run_parallel_interfaces(num_threads, stability_data)
    pd.DataFrame(df_to_export).to_csv(output_location)