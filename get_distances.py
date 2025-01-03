from Bio.PDB import PDBParser
import numpy as np

'''
returns distances between target protein and list of proteins
'''

# inspo from:
# https://bioinformatics.stackexchange.com/questions/783/how-can-we-find-the-distance-between-all-residues-in-a-pdb-file
# https://www.biostars.org/p/374180/

# pdb_path = "/rcfs/projects/proteometer/alphafold_swissprot_pdb"


# best metric for distance is get the distance between target and each of the list and average them out 

def find_mean_distance(structure, target_residue, residue_list):
    residue1 = target_residue
    distances = []
    for residue2 in residue_list:
        distance_to_add = get_pairwise_distance(structure, residue1, residue2)
        # print(distance_to_add)
        distances.append(distance_to_add)
    return np.nanmean(distances)


def find_min_distance(structure, target_residue, residue_list):
    residue1 = target_residue
    distances = []
    for residue2 in residue_list:
        distance_to_add = get_pairwise_distance(structure, residue1, residue2)
        # print(distance_to_add)
        distances.append(distance_to_add)
    return np.nanmin(distances)

def find_min_and_mean_distance(structure, target_residue, residue_list):
    residue1 = target_residue
    distances = []
    for residue2 in residue_list:
        distance_to_add = get_pairwise_distance(structure, residue1, residue2)
        # print(distance_to_add)
        distances.append(distance_to_add)
    return np.nanmin(distances), np.nanmean(distances)




def get_pairwise_distance(structure, residue_num_1, residue_num_2):
    '''
        returns distance between 2 residues, using the CA atom from each residues
    '''
    #print("Begin get_pairwise_distance")
    atom_1 = structure.query('residue_number == @residue_num_1 and atom_name == "CA"') # restrict to only res #1 and alpha carbon 
    atom_2 = structure.query('residue_number == @residue_num_2 and atom_name == "CA"') # restrict to only res #2 and alpha carbon 
    #print("Pre if statement")
    if atom_1['x_coord'].values and atom_1['y_coord'].values and atom_1['z_coord'].values and atom_2['x_coord'].values and atom_2['y_coord'].values and atom_2['z_coord'].values: # if all of these values aren't empty
        #print(atom_1['x_coord'].values, atom_1['y_coord'].values, atom_1['z_coord'].values)
        #print(atom_2['x_coord'].values, atom_2['y_coord'].values, atom_2['z_coord'].values)
        x_p, y_p, z_p = atom_1['x_coord'].values, atom_1['y_coord'].values, atom_1['z_coord'].values
        x_q, y_q, z_q  = atom_2['x_coord'].values, atom_2['y_coord'].values, atom_2['z_coord'].values
        distance = np.sqrt((x_p - x_q)**2 + (y_p - y_q)**2 + (z_p - z_q)**2)
        #print(distance)
    else: 
        distance = np.array([np.nan])
    return distance
