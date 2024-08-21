from Bio.PDB import PDBParser
improt numpy as np

'''
returns distances between target protein and list of proteins
'''

# inspo from:
# https://bioinformatics.stackexchange.com/questions/783/how-can-we-find-the-distance-between-all-residues-in-a-pdb-file
# https://www.biostars.org/p/374180/

# pdb_path = "/rcfs/projects/proteometer/alphafold_swissprot_pdb"


# best metric for distance is get the distance between target and each of the list and average them out 

def find_mean_distances(structure, target_residue, residue_list, pdb_path):
    # create parser SF: this  parser need to be outside of loops
    parser = PDBParser()
    # read structure from file
    structure = parser.get_structure('PHA-L', '1fat.pdb')

    model = structure[0]
    chain = model['A']

    # this example uses only the first residue of a single chain.
    # it is easy to extend this to multiple chains and residues.
    residue1 = target_residue
    distances = []
    for residue2 in residue_list:
        distances.append(get_pairwise_distance(structure, residue1, residue2))
        # ceremove  NAs from distances
    return np.mean(distances)

def get_pairwise_distance(structure, residu_e1, residue_2):
                    try:
                    distance = residue1['CA'] - residue2['CA']
                except KeyError:
                    ## no CA atom, e.g. for H_NAG
                    continue
                if distance < 6:
                    print(residue1, residue2, distance)
    dist = XXX
    return dist