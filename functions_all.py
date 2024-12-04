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


# read in active, binding, and dna binding sites data
functions_data = pd.read_csv("/rcfs/projects/proteometer/UniProt_proteomes/UP000005640_9606_Human/functions_UP000005640_AND_revi_2024_11_19.tsv", sep='\t')


def run_parallel_functions(number_of_threads, stability_data):

    # get all of the unique uniprots
    unique_uniprot_stability = [stability_data[stability_data["protein_acc"]==uniprot_id].copy() for uniprot_id in stability_data["protein_acc"].unique()]

    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(functional_distances, unique_uniprot_stability)
    finish_time = time.perf_counter()
    
    # save the csv and output the start and end times
    print("Program finished in {} - using multiprocessing with {} cores".format(str(datetime.timedelta(seconds=finish_time-start_time)), number_of_threads))
    concatenated_output = pd.concat(output)
    return(concatenated_output)



'''
parses local uniprotKB database. Separates ligand binding, dna binding, and active site fields into a dataframe with rows: [['uniprot_ID','active_sites','active_sites_type','binding_sites','binding_sites_ligand','dna_binding_sites']]
'''

# for each unique uniprotID...
# for uniprot in unique_uniprots:
def functional_distances(uniprot_only_stability, functions_data = functions_data, pickle_output = "/qfs/projects/proteometer/pheno_analysis/functions_pickle_files"):

   # read in interfaces data, uniprot human data, and get uniprot
    uniprot = uniprot_only_stability["protein_acc"].to_list()[0]
    pickle_file_path = f"{pickle_output}/{uniprot}.pkl"
    functions_only_uniprot = functions_data.loc[functions_data['uniprot'].str.contains(uniprot)] # isolate binding sites data to only that uniprot

    if os.path.isfile(pickle_file_path):
        with open(pickle_file_path, 'rb') as handle:
            psp_data = pickle.load(handle)
        return(psp_data)

    # add columns to df
    uniprot_only_stability['closest_binding_site'] = ""
    uniprot_only_stability['closest_binding_site_decription'] = ""
    uniprot_only_stability['closest_binding_site_min_distance'] = np.nan

    uniprot_only_stability['closest_active_site'] = ""
    uniprot_only_stability['closest_active_site_decription'] = ""
    uniprot_only_stability['closest_active_site_min_distance'] = np.nan

    uniprot_only_stability['closest_dna_binding_site'] = ""
    uniprot_only_stability['closest_dna_binding_site_min_distance'] = np.nan

    # parse your structure here
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "-F1-*") # change this to have F1 using pdb file name
    print(pdb_name)
    #print("name of pdb is:", pdb_name)
    if len(pdb_name) != 0:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])
        input_struct = ppdb.df['ATOM']
      
        # for each psp
        for stability_row_index in uniprot_only_stability.index:
            if pd.notna(uniprot_only_stability.loc[stability_row_index,'position']):
                residue_num = int(uniprot_only_stability.loc[stability_row_index,'position']) # finding the residue number of the psp
                active_min_dist = np.inf # make min dist extremely high at first
                binding_min_dist = np.inf
                dna_min_dist = np.inf
                # active site
                if functions_only_uniprot['active_site_residues'].values[0] != "[]":
                    active_site_residues = [int(e[1:]) for e in functions_only_uniprot['active_site_residues'].values[0][1:-1].split(", ")]
                    active_site_descriptions = functions_only_uniprot['active_site_description'].values[0][1:-1].split(", ")
                    if residue_num in active_site_residues:
                        uniprot_only_stability.loc[stability_row_index,'closest_active_site'] = residue_num
                        uniprot_only_stability.loc[stability_row_index,'closest_active_site_decription'] = active_site_descriptions[active_site_residues.index(residue_num)]
                        uniprot_only_stability.loc[stability_row_index,'closest_active_site_min_distance'] = 0.0
                    else:    
                        for active_site in active_site_residues:
                            active_new_min_dist = get_pairwise_distance(input_struct, residue_num, active_site)
                            if active_min_dist > active_new_min_dist:
                                uniprot_only_stability.loc[stability_row_index,'closest_active_site'] = active_site
                                uniprot_only_stability.loc[stability_row_index,'closest_active_site_decription'] = active_site_descriptions[active_site_residues.index(active_site)]
                                uniprot_only_stability.loc[stability_row_index,'closest_active_site_min_distance'] = active_new_min_dist
                                active_min_dist = active_new_min_dist
                                print(active_min_dist)
                                print(active_new_min_dist)

                # binding site
                if functions_only_uniprot['binding_site_residues'].values[0] != "[]":
                    binding_site_residues = [int(e[1:]) for e in functions_only_uniprot['binding_site_residues'].values[0][1:-1].split(", ")]
                    binding_site_descriptions = functions_only_uniprot['binding_site_description'].values[0][1:-1].split(", ")
                    if residue_num in binding_site_residues:
                        uniprot_only_stability.loc[stability_row_index,'closest_binding_site'] = residue_num
                        uniprot_only_stability.loc[stability_row_index,'closest_binding_site_decription'] = binding_site_descriptions[binding_site_residues.index(residue_num)]
                        uniprot_only_stability.loc[stability_row_index,'closest_binding_site_min_distance'] = 0.0
                    else:    
                        for binding_site in binding_site_residues:
                            binding_new_min_dist = get_pairwise_distance(input_struct, residue_num, binding_site)
                            if binding_min_dist > binding_new_min_dist:
                                uniprot_only_stability.loc[stability_row_index,'closest_binding_site'] = binding_site
                                uniprot_only_stability.loc[stability_row_index,'closest_binding_site_decription'] = binding_site_descriptions[binding_site_residues.index(binding_site)]
                                uniprot_only_stability.loc[stability_row_index,'closest_binding_site_min_distance'] = binding_new_min_dist
                                binding_min_dist = binding_new_min_dist
                
                # dna binding site
                if functions_only_uniprot['dna_binding_residues'].values[0] != "[]":
                    dna_binding_site_residues = [int(e[1:]) for e in functions_only_uniprot['dna_binding_residues'].values[0][1:-1].split(", ")]
                    if residue_num in dna_binding_site_residues:
                        uniprot_only_stability.loc[stability_row_index,'closest_dna_binding_site'] = residue_num
                        uniprot_only_stability.loc[stability_row_index,'closest_dna_binding_site_min_distance'] = 0.0
                    else:    
                        for dna_binding_site in dna_binding_site_residues:
                            dna_new_min_dist = get_pairwise_distance(input_struct, residue_num, dna_binding_site)
                            if dna_min_dist > dna_new_min_dist:
                                uniprot_only_stability.loc[stability_row_index,'closest_dna_binding_site'] = dna_binding_site
                                uniprot_only_stability.loc[stability_row_index,'closest_dna_binding_site_min_distance'] = dna_new_min_dist
                                dna_min_dist = dna_new_min_dist

    with open(pickle_file_path, 'wb') as handle:
        pickle.dump(uniprot_only_stability, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return(uniprot_only_stability) 


                    

if __name__ == "__main__":
    num_threads = 64
    stability_data = pd.read_csv(sys.argv[1])
    output_location = sys.argv[2]
    df_to_export = run_parallel_functions(num_threads, stability_data)
    pd.DataFrame(df_to_export).to_csv(output_location)
