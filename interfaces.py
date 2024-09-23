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


def run_parallel_interfaces(number_of_threads, psp_data):
    #adding columns to psp data

    psp_data['closest_interface'] = ""
    psp_data['inside_interface'] = 0
    psp_data['distance_from_interface'] = np.nan

    # get all of the unique uniprots
    # unique_uniprots = psp_data['uniprot_id'].unique()
    unique_uniprot_psp = [psp_data[psp_data["uniprot_id"]==uniprot_id].copy() for uniprot_id in psp_data["uniprot_id"].unique()]

    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        output = pool.map(find_interfaces_per_uniprot, unique_uniprot_psp)
    finish_time = time.perf_counter()
    
    # save the csv and output the start and end times
    print("Program finished in {} - using multiprocessing with {} cores".format(str(datetime.timedelta(seconds=finish_time-start_time)), number_of_threads))
    concatenated_output = pd.concat(output)
    return(concatenated_output)




'''
this function does interfaces calcuations for each uniprot tthat it is given
'''

# for each unique uniprotID...
# for uniprot in unique_uniprots:
def find_interfaces_per_uniprot(psp_only_uniprot):

    # isolate to psp and interfaces in each uniprot
    interfaces_data = pd.read_csv("/rcfs/projects/proteometer/ProtVar/predictions/interfaces/2024.05.28_interface_summary_5A.tsv", delimiter='\t', header=0)
    uniprot = psp_only_uniprot["uniprot_id"].to_list()[0]
    interface_only_uniprot = interfaces_data.loc[(interfaces_data['uniprot_id1'] == uniprot) | (interfaces_data['uniprot_id2'] == uniprot)] # isolate to uniprot in either 1 or 2


    # parse your structure here
    pdb_name = glob.glob("/rcfs/projects/proteometer/alphafold_swissprot_pdb/*" + uniprot + "-F1-*") # change this to have F1 using pdb file name
    print(pdb_name)
    #print("name of pdb is:", pdb_name)
    if len(pdb_name) != 0:  
        ppdb = PandasPdb()  
        ppdb.read_pdb(pdb_name[0])


    # for each psp
        for phosphosite_row_index in psp_only_uniprot.index:
            if pd.notna(psp_only_uniprot.loc[phosphosite_row_index,'res_number']): # only if there is a residue # (this threw an error previously)
                residue_num = int(psp_only_uniprot.loc[phosphosite_row_index,'res_number']) # finding the residue number of the psp
                min_dist = np.inf # make min dist extremely high at first
                #print(residue_num)
                # use the residue # to get the coordinates in space from pdb file
                
                
                for interface_index in interface_only_uniprot.index : # get all the residues in all of the interfaces 
                    if pd.notna(interface_only_uniprot.loc[interface_index,'ifresid1']) & pd.notna(interface_only_uniprot.loc[interface_index,'ifresid2']):
                        if interfaces_data.loc[interface_index,'uniprot_id1'] == uniprot:
                            interface_residues = interface_only_uniprot.loc[interface_index,'ifresid1']
                        elif interfaces_data.loc[interface_index,'uniprot_id2'] == uniprot:
                            interface_residues = interface_only_uniprot.loc[interface_index,'ifresid2']
                        
                        # check if it's inside of a interface
                        #print(interfaces_data.loc[interface_index,'interaction_id']) 
                        interface_residues = [int(e[1:]) for e in interface_residues.split(",")] # remove the first letter from each bc it includes residue type ormat the interface_residues because it's a string
                        print(interface_residues)
                        if residue_num in interface_residues:
                            #print("found inside interface!")
                            psp_only_uniprot.loc[phosphosite_row_index,'inside_interface'] = 1 # if residue is in the interface, put 1 in the inside interface column
                            interface_to_add = (interface_only_uniprot.loc[interface_index,'interaction_id'].split('_'))
                            interface_to_add.remove(uniprot)
                            psp_only_uniprot.loc[phosphosite_row_index,'closest_interface'] = interface_to_add[0] # put unique interfaceID in closest interface
                            psp_only_uniprot.loc[phosphosite_row_index,'distance_from_interface'] = 0.0 
                            break # break because you don't want to contiue looking for interfaces (and therefore overwrite the inside interface and closest interface)
                        else: # if the phosphosite isn't in any interfaces
                            #print("phosphosite isn't in any interfaces")
                            input_struct = ppdb.df['ATOM']
                            #print(input_struct)
                            
                            #print(residue_num, interface_residues)
                            new_dist = find_mean_distances(input_struct, residue_num, interface_residues)
                            print("the new dist is:" , new_dist)
                            if new_dist:
                                if min_dist > new_dist: # if this is the smallest distance so far, replace min_dist with new_dist
                                    interface_to_add = (interface_only_uniprot.loc[interface_index,'interaction_id'].split('_'))
                                    interface_to_add.remove(uniprot)
                                    psp_only_uniprot.loc[phosphosite_row_index,'closest_interface'] = interface_to_add[0]
                                    psp_only_uniprot.loc[phosphosite_row_index,'distance_from_interface'] = new_dist # replace distance_from_interface with min_dist
                                    min_dist = new_dist
                                    print("replaced old dist with", min_dist)
    return(psp_only_uniprot)



                    

if __name__ == "__main__":
    num_threads = 64
    psp_data = pd.read_csv(sys.argv[1])
    output_location = sys.argv[2]
    df_to_export = run_parallel_interfaces(num_threads, psp_data)
    pd.DataFrame(df_to_export).to_csv(output_location)