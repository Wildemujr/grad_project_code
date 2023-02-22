import sys, os
import numpy as np
from Bio.PDB import PDBParser, Selection
from pprint import pprint as pp
import open3d as o3d
import copy


def get_shortest_atom_length(structure_A, structure_B):
    return min(len(structure_A), len(structure_B))



def extract_atom_coordinates(structure):
    
    coordinates = [] 
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())
    
    return coordinates


def alpha_carbon_coords(structure):
    
    CA_coordinates = np.array([]) 
    
    for model in structure:
        for chain in model:
            for residue in chain:
                CA_coordinates = np.append(CA_coordinates, residue['CA'].get_vector())
    
    return CA_coordinates



def main():
    
    CWD = os.getcwd()
    file_A = f"{CWD}/{sys.argv[1]}" 
    file_B = f"{CWD}/{sys.argv[2]}"


    # * Use PDB Parser to read in the PDB files
    parser = PDBParser()
    ref_structure = parser.get_structure("reference", file_A)
    sample_structure = parser.get_structure("sample", file_B)


    # * Get the alpha carbon coordinates for (all residues
    ref_structure_ca = alpha_carbon_coords(ref_structure)
    sample_structure_ca = alpha_carbon_coords(sample_structure)

    pp(ref_structure_ca)


    # * Extractinig the coordinates from the structures objects and 
    # * converting the list of numpy arrays in a single numpy array.
    # ref_struct_coords = np.stack( extract_atom_coordinates(ref_structure) )
    # sample_struct_coords = np.stack( extract_atom_coordinates(sample_structure) )


    # # * Convert generator object to list
    # ref_structure_atoms = list(ref_structure.get_atoms())
    # sample_structure_atoms = list(sample_structure.get_atoms())

    



    # # * Get the shortest length of atoms and use that to index 
    # # * the structures so they are the same length
    # shortest_length = get_shortest_atom_length(ref_structure_atoms, sample_structure_atoms)
    # ref_structure_atoms = ref_structure_atoms[:shortest_length]
    # sample_structure_atoms = sample_structure_atoms[:shortest_length]



 
    return 0 



if __name__ == "__main__":
    main()
    