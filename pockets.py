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


'''
this script takes in psp data and pockets data.  It will check if the phosphosites are located in any pockets, which is the nearest, and the distance to each
usage:   pockets.py <path to phosphosite data> <path to pockets data>

'''


def run_parallel_pockets():
    #adding columns to psp data
    output_location = pd.read_csv(sys.argv[3])

    # get all of the unique uniprots
    unique_uniprots = psp_data['uniprot_id'].unique() 

    start_time = time.perf_counter()
    with Pool(2) as pool:
        pool.map(target = find_pockets_per_uniprot, args = unique_uniprots)
    finish_time = time.perf_counter()
    
    # save the csv and output the start and end times
    psp_data.to_csv(output_location)


'''
this function does pockets calcuations for each uniprot tthat it is given
'''

# for each unique uniprotID...
# for uniprot in unique_uniprots:
def find_pockets_per_uniprot(uniprot):

    # currently putting in the arguments like this bc i dont know how to put in several inputs into multithreading
    psp_data = pd.read_csv(sys.argv[1])
    pockets_data = pd.read_csv(sys.argv[2])

    psp_data['closest_pocket'] = "NaN"
    psp_data['inside_pocket'] = 0
    psp_data['distance_from_pocket'] = "NaN"

    
    # isolate to psp and pockets in each uniprot
    psp_only_uniprot = psp_data[psp_data.uniprot_id == uniprot]
    pocket_only_uniprot = pockets_data[pockets_data.uniprot_id == uniprot]


    # parse your structure here
    pdb_path = "/rcfs/projects/proteometer/alphafold_swissprot_pdb"
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "*")
    print("name of pdb is:", pdb_name)
    if pdb_name:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])


    # for each psp
        for phosphosite_row_index in psp_only_uniprot.index:
            #print(psp_only_uniprot)
            #print(phosphosite_row_index)
            residue_num = psp_only_uniprot.loc[phosphosite_row_index,'res_number'] # finding the residue number of the psp
            min_dist = 100000000000000000000000000000000 # make min dist extremely high at first
            #print(residue_num)
            # use the residue # to get the coordinates in space from pdb file
            
            
            for pocket_index in pocket_only_uniprot.index : # get all the residues in all of the pockets 
                pocket_residues = pocket_only_uniprot.loc[pocket_index,'pocket_resid']

                # check if it's inside of a pocket
                pocket_residues = pocket_residues[1:-1].split(",") # format the pocket_residues because it's a string
                #print(pocket_residues)
                if residue_num in pocket_residues:
                    psp_data.loc[phosphosite_row_index,'inside_pocket'] = 1 # if residue is in the pocket, put 1 in the inside pocket column
                    psp_data.loc[phosphosite_row_index,'closest_pocket'] = pocket_only_uniprot.loc[pocket_index,'full_id'] # put unique pocketID in closest pocket
                    psp_data.loc[phosphosite_row_index,'distance_from_pocket'] = 0 
                    break # break because you don't want to contiue looking for pockets (and therefore overwrite the inside pocket and closest pocket)

                if psp_data.loc[phosphosite_row_index,'inside_pocket'] == 0: # if the phosphosite isn't in any pockets
                    print("phosphosite isn't in any pockets")
                    input_struct = ppdb.df['ATOM']
                    #print(input_struct)
                    new_dist = find_mean_distances(input_struct, residue_num, pocket_residues)
                    if residue_num:
                        if min_dist > new_dist: # if this is the smallest distance so far, replace min_dist with new_dist
                            psp_data.loc[phosphosite_row_index,'closest_pocket'] = pocket_only_uniprot.loc[pocket_index,'full_id'] # put unique pocketID in closest pocket
                            psp_data.loc[phosphosite_row_index,'distance_from_pocket'] = new_dist # replace distance_from_pocket with min_dist
                            min_dist = new_dist 
                            print("added smallest distance:", min_dist)
                
    else: # if we can't find the pdb file
        for phosphosite_row_index in psp_only_uniprot.index:
            psp_data.loc[phosphosite_row_index,'inside_pocket'] = 'NaN' 
            psp_data.loc[phosphosite_row_index,'closest_pocket'] = 'NaN' 
            psp_data.loc[phosphosite_row_index,'distance_from_pocket'] = 'NaN'






                    

if __name__ == "__main__":
    run_parallel_pockets()
    # first way, using multiprocessing
    print("Program finished in {} seconds - using multiprocessing".format(finish_time-start_time))
    print("---")