import pandas as pd
import sys
sys.path.insert(1, '/people/imal967/git_repos/ProCaliper')
import procaliper as am
import glob
import numpy as np 
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
sites_data = pd.read_csv("/people/imal967/git_repos/pheno_analysis/binding_active_dna.csv")


def run_parallel_sites(number_of_threads, human_db):

    # get all of the unique uniprots
    unique_uniprot_psp = [human_db[human_db["Entry"]==uniprot_id].copy() for uniprot_id in human_db["Entry"].unique()]

    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(binding_active_dna_distances, unique_uniprot_psp)
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
def binding_active_dna_distances(uniprot_only_humandb, sites_data = sites_data, distances_output = "/qfs/projects/proteometer/pheno_analysis/active_binding_dna_distances"):

   # read in interfaces data, uniprot human data, and get uniprot
    uniprot = uniprot_only_humandb["protein_acc"].to_list()[0]
    intermediate_files_path = f"{distances_output}/{uniprot}_site_distances.csv"
    sites_only_uniprot = sites_data.loc[sites_data['uniprot'] == uniprot] # isolate binding sites data to only that uniprot

    if os.path.isfile(intermediate_files_path):
        distances_data = pd.read_csv(intermediate_files_path)
        return(distances_data)

    # add columns to df
    uniprot_only_humandb['closest_binding_site'] = ""
    uniprot_only_humandb['closest_binding_site_decription'] = ""
    uniprot_only_humandb['closest_binding_site_min_distance'] = ""

    uniprot_only_humandb['closest_active_site'] = ""
    uniprot_only_humandb['closest_active_site_decription'] = ""
    uniprot_only_humandb['closest_active_site_min_distance'] = ""


    uniprot_only_humandb['closest_dna_binding_site'] = ""
    uniprot_only_humandb['closest_dna_binding_site_min_distance'] = ""



    # parse your structure here
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "-F1-*") # change this to have F1 using pdb file name
    print(pdb_name)
    #print("name of pdb is:", pdb_name)
    if len(pdb_name) != 0:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])
        input_struct = ppdb.df['ATOM']

        
    # for each psp
        for humandb_row_index in uniprot_only_humandb.index:
            residue_num = int(uniprot_only_humandb.loc[humandb_row_index,'position']) # finding the residue number of the psp
            active_min_dist = np.inf # make min dist extremely high at first
            binding_min_dist = np.inf
            dna_min_dist = np.inf
            # active site
            if sites_only_uniprot['active_site_residues'].values[0] != "[]":
                active_site_residues = [int(e[1:]) for e in sites_only_uniprot['active_site_residues'].values[0][1:-1].split(", ")]
                active_site_descriptions = sites_only_uniprot['active_site_description'].values[0][1:-1].split(", ")
                for active_site in active_site_residues:
                    active_new_min_dist = find_min_distance(input_struct, residue_num, active_site)
                    if active_min_dist > active_new_min_dist:
                        uniprot_only_humandb.loc[humandb_row_index,'closest_active_site'] = active_site
                        uniprot_only_humandb.loc[humandb_row_index,'closest_active_site_decription'] = active_site_descriptions[active_site_residues.index(active_site)]
                        uniprot_only_humandb.loc[humandb_row_index,'closest_active_site_min_distance'] = active_new_min_dist
                        active_min_dist = active_new_min_dist
                        print(active_min_dist)
                        print(active_new_min_dist)

            # binding site
            if sites_only_uniprot['binding_site_residues'].values[0] != "[]":
                binding_site_residues = [int(e[1:]) for e in sites_only_uniprot['binding_site_residues'].values[0][1:-1].split(", ")]
                binding_site_descriptions = sites_only_uniprot['binding_site_description'].values[0][1:-1].split(", ")
                for binding_site in binding_site_residues:
                    binding_new_min_dist = find_min_distance(input_struct, residue_num, binding_site)
                    if binding_min_dist > binding_new_min_dist:
                        uniprot_only_humandb.loc[humandb_row_index,'closest_binding_site'] = binding_site
                        uniprot_only_humandb.loc[humandb_row_index,'closest_binding_site_decription'] = binding_site_descriptions[binding_site_residues.index(binding_site)]
                        uniprot_only_humandb.loc[humandb_row_index,'closest_binding_site_min_distance'] = binding_new_min_dist
                        binding_min_dist = binding_new_min_dist
            
            # dna binding site
            if sites_only_uniprot['dna_binding_residues'].values[0] != "[]":
                dna_site_residues = [int(e[1:]) for e in sites_only_uniprot['dna_binding_residues'].values[0][1:-1].split(", ")]
                for dna_site in dna_site_residues:
                    dna_new_min_dist = find_min_distance(input_struct, residue_num, dna_site)
                    if dna_min_dist > dna_new_min_dist:
                        uniprot_only_humandb.loc[humandb_row_index,'closest_dna_binding_site'] = dna_site
                        uniprot_only_humandb.loc[humandb_row_index,'closest_dna_binding_site_min_distance'] = dna_new_min_dist
                        dna_min_dist = dna_new_min_dist

    uniprot_only_humandb.to_csv(intermediate_files_path)
    return(uniprot_only_humandb) 


                    

if __name__ == "__main__":
    num_threads = 64
    human_db = pd.read_csv("/rcfs/projects/proteometer/human_proteome_precaculated/stability_precaculated.csv")
    output_location = "/people/imal967/git_repos/pheno_analysis/test_distances.csv"
    df_to_export = binding_active_dna_distances(num_threads, human_db)

    pd.DataFrame(df_to_export).to_csv(output_location)