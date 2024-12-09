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

'''
Runs Procaliper sasa and charge calculations in parallel on csv with multiple uniprotIDs
'''

def run_parallel_procaliper(number_of_threads, human_db):

    # get all of the unique uniprots
    unique_uniprot_psp = [human_db[human_db["Entry"]==uniprot_id].copy() for uniprot_id in human_db["Entry"].unique()]

    start_time = time.perf_counter()
    with Pool(number_of_threads) as pool:
        pool.map(find_sasa_charge, unique_uniprot_psp)
        # output = pool.map(find_sasa_charge, to_run)

    finish_time = time.perf_counter()
    
    # save the csv and output the start and end times
    print("Program finished in {} - using multiprocessing with {} cores".format(str(datetime.timedelta(seconds=finish_time-start_time)), number_of_threads))
    # concatenated_output = pd.concat(output)
    # return(concatenated_output)




'''
Calculates sasa values and charge values for given csv
'''

# for each unique uniprotID...
# for uniprot in unique_uniprots:
def find_sasa_charge(uniprotKB, sasa_output = "/qfs/projects/proteometer/pheno_analysis/sasa_csv", charge_output = "/qfs/projects/proteometer/pheno_analysis/charge_csv"):

    # isolate to psp 
    uniprot = uniprotKB["Entry"].to_list()[0]

    if not os.path.isfile(f"{sasa_output}/{uniprot}_sasa.csv"):

    # create protein structur
        protein_object = []
        try:
            protein_object = am.Protein.from_uniprot_row(uniprotKB.iloc[0])# feed in the 1st (and only) row of that dataframe
        except ValueError:
            print(f"Couldn't retrieve uniprotID: {uniprot}")




        if protein_object:
            print(f"processing {uniprot}")
            pdb_name = glob.glob(f"/rcfs/projects/proteometer/alphafold_swissprot_pdb/*{protein_object.data['entry']}*.pdb")
            #print("name of pdb is:", pdb_name)
            if len(pdb_name) != 0:  
                protein_object.register_local_pdb(pdb_name[0]) # add pdb path to protein object

                # calculate sasa and charge using Procaliper
                sasa_df = pd.DataFrame(protein_object.get_sasa())
                charge_df = pd.DataFrame(protein_object.get_charge())
                charge_df['UNIPROT'] = uniprot
                sasa_df['UNIPROT'] = uniprot
                charge_df['RES_NUM'] = charge_df['residue_number']
                sasa_df['RES_NUM'] = sasa_df['residue_number']

                # save full protein sasa and charge
                sasa_df.to_csv(f"{sasa_output}/{uniprot}_sasa.csv")
                charge_df.to_csv(f"{charge_output}/{uniprot}_charge.csv")

                # merge the dataframes
                # sasa_psp = pd.merge(uniprotKB, sasa_df[['all_sasa_value', 'atom_sasa_values','RES_NUM']], how = 'left', on= 'RES_NUM')
                # sasa_charge_psp = pd.merge(sasa_psp, charge_df[['charge', 'charge_method','RES_NUM']], how = 'left', on= 'RES_NUM')


        #     with open(pickle_file_path, 'wb') as handle:
        #         pickle.dump(sasa_charge_psp, handle, protocol=pickle.HIGHEST_PROTOCOL)

    
        # return(sasa_charge_psp) 



                    

if __name__ == "__main__":
    num_threads = 64
    human_db = pd.read_csv(sys.argv[1], delimiter='\t', header=0)
    output_location = sys.argv[2]
    # df_to_export = run_parallel_procaliper(num_threads, human_db)
    run_parallel_procaliper(num_threads, human_db)

    # run_parallel_procaliper(num_threads, human_db)
    # pd.DataFrame(df_to_export).to_csv(output_location)