import sys, os
import numpy as np
from Bio.PDB import PDBParser, Selection, PDBIO
from pprint import pprint as pp
import copy
import mypy
from typing import List
from math import sqrt




def get_shortest_atom_length(structure_A, structure_B):
    return min(len(structure_A), len(structure_B))


def calculate_rmsd(A, B):
    # * Ensure the matrices are the same size. Correct if not.
    try:
        assert A.shape == B.shape
        N = A.shape[0]
    except AssertionError:
        min_dim = min(A.shape, B.shape)[0]
        A_trunc, B_trunc = A[:min_dim, :], B[:min_dim, :]
        N = min_dim
        
    return sqrt(np.sum((A_trunc - B_trunc)**2) / N)



def extract_atom_coordinates(structure):
    
    coordinates = [] 
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())
    
    return coordinates


def all_atom_coords(structure):
    
    all_coordinates = np.empty((0, 3), dtype=np.float32)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord = atom.get_coord()
                    all_coordinates = np.append(all_coordinates, [coord], axis=0)
    
    return all_coordinates


def apply_atom_coords(structure, new_coordinates):
    
    i = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.set_coord(new_coordinates[i,:])
                    i+=1
    return 0


def alpha_carbon_coords(structure):
    
    CA_coordinates = np.array([]) 
    
    for model in structure:
        for chain in model:
            for residue in chain:
                CA_coordinates = np.append(CA_coordinates, residue['CA'].get_vector())
    
    return CA_coordinates


def kabsch_aln(A, B, scale = False):
    
    # * Ensure the matrices are the same size. Correct if not.
    try:
        assert A.shape == B.shape
        N = A.shape[0]
        print("Matrices are the same size. Moving on...\n")
        
    except AssertionError:
        print("Matrices are not the same size. Fixing now...\n")
        min_dim = min(A.shape, B.shape)[0]
        A_trunc, B_trunc = A[:min_dim, :], B[:min_dim, :]
        N = min_dim
    

    # * Calculate the centroid (using column mean) of each matrix
    # * You'll have the mean for each x, y, and z.
    centroid_A, centroid_A_trunc = np.mean(A, axis=0), np.mean(A_trunc, axis=0)
    centroid_B, centroid_B_trunc = np.mean(B, axis=0), np.mean(B_trunc, axis=0)
    
    # * Center the matrices conducting an elemenet-wise subtraction 
    # * of the centroid from each row of the matrix.
    A_trunc_centered = A_trunc - centroid_A_trunc
    B_trunc_centered = B_trunc - centroid_B_trunc
    
    # * Calculate the covariance matrix. Further, check whether
    # * or not to scale the matrices. @ sign here is shorthand for
    # * referencing matrix multiplication. We're essentially calculating
    # * the dot product between the columns of the two matrices.
    if scale: 
        H = ( np.transpose(A_trunc_centered) @ B_trunc_centered ) / N
    else:
        H = ( np.transpose(A_trunc_centered) @ B_trunc_centered )
    

    # * Calculate the SVD of the covariance matrix to obtain the
    # * rotation matrix.
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ np.identity(3) @ U.T
    d = np.linalg.det(Vt.T @ U.T)
    
       
    # * See if the determinate is negative.
    # * If so, then correct our rotation matrix to ensure a right-handed coordinate system.
    # * If the determinant is positive, then we are good to go and do nothing.
    if d < 0:
        S = np.identity(3)
        S[2,:] *= -1
        R = Vt.T @ S @ U.T
    
    # * If scale is true, then calculate the scale factor and translation vector.
    # * if scale is false, then set the scale factor to 1 and the translation vector
    # * is the difference between centroid A and the rotated centroid B.
    if scale:
        varA = np.var(A, axis=0).sum()
        c = (1 / varA) * (np.sum(S))
        t = centroid_A_trunc - (c * R @ centroid_B_trunc)
        
    else:
        c = 1
        t = centroid_A_trunc - R @ centroid_B_trunc
        
    return c, R, t




def main():
    
    CWD = os.getcwd()
    file_A = f"{CWD}/{sys.argv[1]}" 
    file_B = f"{CWD}/{sys.argv[2]}"


    # * Use PDB Parser to read in the PDB files
    parser = PDBParser()
    ref_structure = parser.get_structure("reference", file_A)
    sample_structure = parser.get_structure("sample", file_B)

    # * Get the alpha carbon coordinates for all residues
    ref_structure_ca = alpha_carbon_coords(ref_structure)
    sample_structure_ca = alpha_carbon_coords(sample_structure)

    # * Convert the numpy list of vectors to a 2D numpy array.
    # * Then, reshape the array to be a matrix with 3 rows and N columns.
    refStruct_ca_vectors: np.ndarray([List[float]]) = np.array([ list(vec) for vec in ref_structure_ca] )
    refStruct_ca_matrix: np.ndarray([List[float]]) = refStruct_ca_vectors.reshape(-1, 3)

    sampleStruct_ca_vectors: np.ndarray([List[float]]) = np.array([ list(vec) for vec in sample_structure_ca] )
    sampleStruct_ca_matrix: np.ndarray([List[float]]) = sampleStruct_ca_vectors.reshape(-1, 3)

    iterations = 100
    i = 0
    best_RMSD = 10000
    previous_RMSD = 10000

    scale_factor, R, t = kabsch_aln(refStruct_ca_matrix, sampleStruct_ca_matrix, True)
    ref_structure_coords = all_atom_coords(ref_structure)
    sample_structure_coords = all_atom_coords(sample_structure)     


    while i < iterations:

        # * Perform the Kabsch algorithm and apply the transformation to the sample structure.
        scale_factor, R, t = kabsch_aln(refStruct_ca_matrix, sampleStruct_ca_matrix, True)
        # kabsch_aln(refStruct_ca_matrix, sampleStruct_ca_matrix, False)
        sample_struct_transformed = scale_factor * (sampleStruct_ca_matrix @ R) + t


        # * Calculate the RMSD between the two structures.
        rmsd = calculate_rmsd(refStruct_ca_matrix, sample_struct_transformed)
        print(f"RMSD: {rmsd}")
    
        if rmsd < best_RMSD and rmsd < previous_RMSD:
            best_RMSD = rmsd
            previous_RMSD = rmsd

            # * Apply the transformation matrix to the sample structure.
            sample_structure.transform(R,t)

            # * Output the transformed structure to a PDB file.
            io = PDBIO()
            io.set_structure(sample_structure)
            io.save(f"transformed_structure_scaled_{i}.pdb")
            
            sample_structure_coords = all_atom_coords(sample_structure) @ R
            apply_atom_coords(sample_structure, sample_structure_coords)
            i+=1
            
        else:
            sample_structure_coords = all_atom_coords(sample_structure) @ R
            apply_atom_coords(sample_structure, sample_structure_coords) 
            i+=1
    

    return 0 



if __name__ == "__main__":
    main()
    