from Bio.PDB import PDBParser

'''
returns distances between target protein and list of proteins
'''

# inspo from:
# https://bioinformatics.stackexchange.com/questions/783/how-can-we-find-the-distance-between-all-residues-in-a-pdb-file
# https://www.biostars.org/p/374180/

# pdb_path = "/rcfs/projects/proteometer/alphafold_swissprot_pdb"


def find_distances(uniprot_id, target_residue, residue_list, pdb_path):
    # create parser
    parser = PDBParser()
    # read structure from file
    structure = parser.get_structure('PHA-L', '1fat.pdb')

    model = structure[0]
    chain = model['A']

    # this example uses only the first residue of a single chain.
    # it is easy to extend this to multiple chains and residues.
    for residue1 in chain:
        for residue2 in chain:
            if residue1 != residue2:
                # compute distance between CA atoms
                try:
                    distance = residue1['CA'] - residue2['CA']
                except KeyError:
                    ## no CA atom, e.g. for H_NAG
                    continue
                if distance < 6:
                    print(residue1, residue2, distance)
            # stop after first residue
            break
