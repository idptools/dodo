# functionality for trying to add in all atoms to structures. 
# probably won't be perfect, but could help with visualization. 
import numpy as np
from dodo.parameters import AMINO_ACID_ATOMS_AF2
from scipy.spatial.transform import Rotation
from dodo.pdb import array



def approximate_amino_acid_coordinates(amino_acid, alpha_carbon_coordinate, coord_dict=AMINO_ACID_ATOMS_AF2):
    """
    Approximate the 3D coordinates of other atoms in a specific amino acid based on its identity
    and the alpha carbon's position.

    Parameters
    ------------
    amino_acid : str 
        The identity of the amino acid (e.g., 'A', 'L').
    alpha_carbon_coordinate : tuple 
        3D coordinates of the alpha carbon (Ca) (x, y, z).

    Returns:
        dict: A dictionary containing approximate 3D coordinates of common atoms (Ca, C, O, N, H).
    """
    if amino_acid not in coord_dict:
        raise ValueError(f"Amino acid '{amino_acid}' not found in the database.")

    if len(alpha_carbon_coordinate) != 3:
        raise ValueError("Alpha carbon coordinate must be a 3D coordinate (x, y, z).")

    # Get pre-defined coordinates for the specified amino acid
    amino_acid_atoms = coord_dict[amino_acid]

    # Calculate the positions of atoms relative to the alpha carbon
    relative_positions = {atom: alpha_carbon_coordinate + np.array(coord) for atom, coord in amino_acid_atoms.items()}

    return relative_positions


def rotate_translate_with_distance_preservation(coords_list, original_alpha_coords, new_alpha_coords):
    '''
    Function to try to use original atomic coordinates, original coordinates for alpha carbons, and 
    new coordinates for alpha carbons to place the atoms from a protein back such that the angle of
    each atom relative to CA is preserved and the distance is also preserved. Kind of a hack to get 
    the atoms back to where they should be after a rotation. Should work OK for the IDRs...


    '''
    # Initialize an empty list to store the transformed coordinates
    transformed_coords_list = []

    for i in range(len(coords_list)):
        amino_acid_coords = coords_list[i]
        alpha_coords = original_alpha_coords[i]
        new_alpha = new_alpha_coords[i]

        # Calculate the translation vector from the original alpha carbon to the new alpha carbon
        translation_vector = new_alpha - alpha_coords

        # Initialize a dictionary to store the translated and scaled atom coordinates
        transformed_amino_acid = {}

        # Calculate the translation factor
        translation_factor = np.linalg.norm(translation_vector)

        # Translate the amino acid's atoms based on the translation vector
        for atom, coord in amino_acid_coords.items():
            transformed_amino_acid[atom] = coord + translation_vector

        # Append the translated coordinates to the list
        transformed_coords_list.append(transformed_amino_acid)

    return transformed_coords_list


# because I'm just translating FDs, going to try to input their positions and use the relative
# position of the FDs to populate the atomic positions of the IDRs. 

def calculate_rotation_matrix(p1, p2, q1, q2):
    """
    Calculate a rotation matrix to rotate vector p1p2 to q1q2.
    :param p1: Initial position vector 1.
    :param p2: Initial position vector 2.
    :param q1: Target position vector 1.
    :param q2: Target position vector 2.
    :return: Rotation matrix.
    """
    v1 = p2 - p1
    v2 = q2 - q1
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)

    axis = np.cross(v1, v2)
    angle = np.arccos(np.dot(v1, v2))

    rotation_matrix = Rotation.from_rotvec(angle * axis)
    return rotation_matrix.as_matrix()

def align_amino_acids(amino_acid_1_coords, amino_acid_1_initial_coords, amino_acid_2_initial_coords, amino_acid_2_alpha_final):
    # Calculate the translation vector to move amino_acid_2's alpha carbon to the final position
    translation_vector = amino_acid_2_alpha_final - amino_acid_2_initial_coords['CA']

    # Translate amino_acid_2_initial_coords
    translated_aa2_coords = {
        atom_name: atom_coord + translation_vector for atom_name, atom_coord in amino_acid_2_initial_coords.items()
    }

    # Calculate the rotation matrix
    rotation_matrix = calculate_rotation_matrix(
        amino_acid_1_initial_coords['CA'], translated_aa2_coords['CA'], amino_acid_1_coords['CA'], translated_aa2_coords['N']
    )

    # Apply the rotation to all coordinates of amino_acid_2
    aligned_aa2_coords = {
        atom_name: np.dot(rotation_matrix, atom_coord - translated_aa2_coords['CA']) + amino_acid_2_alpha_final
        for atom_name, atom_coord in translated_aa2_coords.items()
    }

    return aligned_aa2_coords


def translate_single_amino_acid_fd(initial_coords, final_alpha_coords):
    # Calculate the translation vector
    final_alpha_coords=np.array(final_alpha_coords)
    initial_CA = initial_coords['CA']
    initial_CA_arr = np.array(initial_CA)

    translation_vector = final_alpha_coords - initial_CA_arr

    # Translate all atoms based on the translation vector
    translated_coords = {
        atom_name: atom_coord + translation_vector for atom_name, atom_coord in initial_coords.items()
    }
    return translated_coords




def add_C_terminal_IDR_atoms(pdb_dict, idr_name):
    #need to take out IDR name later and make dimilar function for N terminal IDR...
    reg_interest=p53_dict['regions_dict'][idr_name]
    for aaind in range(reg_interest[0], reg_interest[1]):
        established_aa_ind = aaind-1
        final_established_all_atom = pdb_dict['final_coords_atoms_added'][established_aa_ind]
        starting_established_all_atom = pdb_dict['coord_by_aa_ind_with_name'][established_aa_ind]
        amino_acid_2_initial=pdb_dict['coord_by_aa_ind_with_name'][aaind]
        amino_acid_2_ca = np.array(pdb_dict['ca_coords'][aaind])
        new_all_atom_coords=align_amino_acids(final_established_all_atom, starting_established_all_atom,
                amino_acid_2_initial, amino_acid_2_ca)
        cur_all_atom=pdb_dict['final_coords_atoms_added']
        cur_all_atom[aaind]=new_all_atom_coords
        pdb_dict['final_coords_atoms_added']=cur_all_atom  
    return pdb_dict  


def add_all_atom_info(pdb_dict, regions_dict):
    '''
    function to add in all_atom_info after a structure has been translated.
    Basically it will translate the coordinates of folded domains across space
    to match the new coordinates of the CA atoms. For IDRs, because their relative
    orientation is almost certanily not what it originally was, just guesses using
    idealized bond lengths and angles. 
    This will be the final function for adding in atoms... not going great but I guess we will see. 

    # the all atom stuff for the FDs as they are moved is saved in 
    pdb_dict['final_coords_atoms_added']. Need to add the first FD to 
    the dict too because everything is moved relative to that, so the 
    all atom coords should be good there. Once we have that, we just need
    to figure out how to get the IDR atoms added back. Pretty close!

    *** Really should just use more basic approach.
        1. determine change in initial angle to final angle for a pair of alpha carbons. 
        2. Take the atoms, translate them the same distance as the alpha carbon, then
        change their angle to match Ca1->Ca2 changes.
        3. repeat for all the atoms.
        4. For IDRs, just use idealized bond lengths and angles.*If there is no FD to anchor*


    '''
    # get all ca coords
    ca_coords = pdb_dict['ca_coords']
    pass