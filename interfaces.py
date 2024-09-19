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
this script takes in psp data and interfaces data.  It will check if the phosphosites are located in any interfaces, which is the nearest, and the distance to each
usage:   interfaces.py <path to phosphosite data> <path to interfaces data>

'''

# currently putting these outside the functions because i dont' know how to put this in with multithreading. but this is pretty bad practice
psp_data = pd.read_csv(sys.argv[1])
interfaces_data = pd.read_csv(sys.argv[2])
output_location = sys.argv[3]

def run_parallel_interfaces():
    #adding columns to psp data

    psp_data['closest_interface'] = "NaN"
    psp_data['inside_interface'] = 0
    psp_data['distance_from_interface'] = "NaN"

    # get all of the unique uniprots
    unique_uniprots = psp_data['uniprot_id'].unique() 

    start_time = time.perf_counter()
    with Pool(64) as pool:
        output = pool.map(find_interfaces_per_uniprot, unique_uniprots)
    finish_time = time.perf_counter()
    
    # save the csv and output the start and end times
    print("Program finished in {} seconds - using multiprocessing".format(finish_time-start_time))
    return(output)

'''
this function does interfaces calcuations for each uniprot tthat it is given
'''

# for each unique uniprotID...
# for uniprot in unique_uniprots:
def find_interfaces_per_uniprot(uniprot):

    # isolate to psp and interfaces in each uniprot
    psp_only_uniprot = psp_data[psp_data.uniprot_id == uniprot]
    interface_only_uniprot = interfaces_data.loc[(interfaces_data['uniprot_id1'] == uniprot) | (interfaces_data['uniprot_id2'] == uniprot)] # isolate to uniprot in either 1 or 2


    # parse your structure here
    pdb_path = "/rcfs/projects/proteometer/alphafold_swissprot_pdb"
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "*")
    #print("name of pdb is:", pdb_name)
    if pdb_name:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])


    # for each psp
        for phosphosite_row_index in psp_only_uniprot.index:
            residue_num = psp_only_uniprot.loc[phosphosite_row_index,'res_number'] # finding the residue number of the psp
            min_dist = 100000000000000000000000000000000 # make min dist extremely high at first
            #print(residue_num)
            # use the residue # to get the coordinates in space from pdb file
            
            
            for interface_index in interface_only_uniprot.index : # get all the residues in all of the interfaces 
                if pd.notna(interface_only_uniprot.loc[interface_index,'ifresid1']) & pd.notna(interface_only_uniprot.loc[interface_index,'ifresid1']):
                    if interfaces_data.loc[interface_index,'uniprot_id1'] == uniprot:
                        interface_residues = interface_only_uniprot.loc[interface_index,'ifresid1']
                    elif interfaces_data.loc[interface_index,'uniprot_id2'] == uniprot:
                        interface_residues = interface_only_uniprot.loc[interface_index,'ifresid2']
                    
                    # check if it's inside of a interface
                    interface_residues = interface_residues[1:-1].split(",") # format the interface_residues because it's a string
                    interface_residues = [e[1:] for e in interface_residues] # remove the first letter from each bc it includes residue type
                    #print(interface_residues)
                    if residue_num in interface_residues:
                        psp_data.loc[phosphosite_row_index,'inside_interface'] = 1 # if residue is in the interface, put 1 in the inside interface column
                        psp_data.loc[phosphosite_row_index,'closest_interface'] = interface_only_uniprot.loc[interface_index,'interaction_id'] # put unique interfaceID in closest interface
                        psp_data.loc[phosphosite_row_index,'distance_from_interface'] = 0 
                        break # break because you don't want to contiue looking for interfaces (and therefore overwrite the inside interface and closest interface)

                    if psp_data.loc[phosphosite_row_index,'inside_interface'] == 0: # if the phosphosite isn't in any interfaces
                        #print("phosphosite isn't in any interfaces")
                        input_struct = ppdb.df['ATOM']
                        #print(input_struct)
                        new_dist = find_mean_distances(input_struct, residue_num, interface_residues)
                        if residue_num:
                            if min_dist > new_dist: # if this is the smallest distance so far, replace min_dist with new_dist
                                psp_data.loc[phosphosite_row_index,'closest_interface'] = interface_only_uniprot.loc[interface_index,'interaction_id'] # put unique interfaceID in closest interface
                                psp_data.loc[phosphosite_row_index,'distance_from_interface'] = new_dist # replace distance_from_interface with min_dist
                                min_dist = new_dist 
                                print("added smallest distance:", min_dist)
                                print("the interface is:", interface_only_uniprot.loc[interface_index,'interaction_id'])
                
    else: # if we can't find the pdb file
        for phosphosite_row_index in psp_only_uniprot.index:
            psp_data.loc[phosphosite_row_index,'inside_interface'] = 'NaN' 
            psp_data.loc[phosphosite_row_index,'closest_interface'] = 'NaN' 
            psp_data.loc[phosphosite_row_index,'distance_from_interface'] = 'NaN'

    return(psp_data)




                    

if __name__ == "__main__":
    list_to_export = run_parallel_pockets()
    pd.DataFrame(list_to_export[1]).to_csv(output_location)