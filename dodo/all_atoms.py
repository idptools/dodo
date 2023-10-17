# functionality for trying to add in all atoms to structures. 
# probably won't be perfect, but could help with visualization. 

'''
Current status is it's hot garbage. Going to disable access for 
user facing functionality for the time being. 
'''

import numpy as np
from dodo.parameters import AMINO_ACID_ATOMS_AF2
from dodo.pdb_tools import array
import math

def get_res_dist(xyz1, xyz2):
    '''
    get distance between two coordinates in 3D space

    Parameters
    ----------
    xyz1 : list
        list of x,y,z coordinates for first point
    xyz2 : list
        list of x,y,z coordinates for second point

    Returns
    -------
    float
        distance between two points
    '''
    xdiff = xyz1[0]-xyz2[0]
    ydiff = xyz1[1]-xyz2[1]
    zdiff = xyz1[2]-xyz2[2]
    return math.sqrt((xdiff*xdiff)+(ydiff*ydiff)+(zdiff*zdiff))


def approximate_amino_acid_coordinates(amino_acid, alpha_carbon_coordinate, coord_dict=AMINO_ACID_ATOMS_AF2):
    """
    SORT OF USELESS 

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


def translate_atoms(original_alpha_carbon, final_alpha_carbon, amino_acid_atoms):
    # Calculate the translation vector
    translation_vector = final_alpha_carbon - original_alpha_carbon
    
    # Translate all atoms by the same vector
    translated_atoms = amino_acid_atoms + translation_vector
    
    return translated_atoms

def rotate_atoms(original_alpha_carbon_1, original_alpha_carbon_2, final_alpha_carbon_1, final_alpha_carbon_2, amino_acid_atoms):
    # Calculate the original and final vectors between alpha carbons
    v1_original = original_alpha_carbon_2 - original_alpha_carbon_1
    v1_final = final_alpha_carbon_2 - final_alpha_carbon_1

    # Calculate the rotation axis and angle
    rotation_axis = np.cross(v1_original, v1_final)
    norm_rotation_axis = np.linalg.norm(rotation_axis)
    if norm_rotation_axis < 1e-8:
        # No rotation is needed if the vectors are parallel or anti-parallel
        return amino_acid_atoms

    rotation_axis /= norm_rotation_axis
    cos_theta = np.dot(v1_original, v1_final) / (np.linalg.norm(v1_original) * np.linalg.norm(v1_final))
    sin_theta = np.sqrt(1 - cos_theta**2)
    theta = np.arccos(cos_theta)

    # Create the rotation matrix
    rotation_matrix = np.zeros((3, 3))
    rotation_matrix[0, 0] = cos_theta + rotation_axis[0]**2 * (1 - cos_theta)
    rotation_matrix[0, 1] = rotation_axis[0] * rotation_axis[1] * (1 - cos_theta) - rotation_axis[2] * sin_theta
    rotation_matrix[0, 2] = rotation_axis[0] * rotation_axis[2] * (1 - cos_theta) + rotation_axis[1] * sin_theta
    rotation_matrix[1, 0] = rotation_axis[0] * rotation_axis[1] * (1 - cos_theta) + rotation_axis[2] * sin_theta
    rotation_matrix[1, 1] = cos_theta + rotation_axis[1]**2 * (1 - cos_theta)
    rotation_matrix[1, 2] = rotation_axis[1] * rotation_axis[2] * (1 - cos_theta) - rotation_axis[0] * sin_theta
    rotation_matrix[2, 0] = rotation_axis[0] * rotation_axis[2] * (1 - cos_theta) - rotation_axis[1] * sin_theta
    rotation_matrix[2, 1] = rotation_axis[1] * rotation_axis[2] * (1 - cos_theta) + rotation_axis[0] * sin_theta
    rotation_matrix[2, 2] = cos_theta + rotation_axis[2]**2 * (1 - cos_theta)

    # Apply the rotation to amino_acid_atoms
    rotated_atoms = np.dot(amino_acid_atoms - original_alpha_carbon_1, rotation_matrix.T) + final_alpha_carbon_1

    return rotated_atoms

def calculate_polar_angle(point1, point2):
    """
    Calculate the polar angle in radians between two 3D coordinates and the positive x-axis.

    Args:
        point1 (tuple or list): The coordinates of the first point (x1, y1, z1).
        point2 (tuple or list): The coordinates of the second point (x2, y2, z2).

    Returns:
        float: The polar angle in radians between the two points and the positive x-axis.
    """
    if len(point1) != 3 or len(point2) != 3:
        raise ValueError("Both points must have 3 coordinates (x, y, z)")

    x1, y1, z1 = point1
    x2, y2, z2 = point2

    # Calculate the differences in x and y coordinates
    delta_x = x2 - x1
    delta_y = y2 - y1

    # Calculate the polar angle using arctan2
    polar_angle = math.atan2(delta_y, delta_x)

    return polar_angle

def calculate_azimuthal_angle(point1, point2):
    """
    Calculate the azimuthal angle in radians between two 3D coordinates and the positive x-axis in the xy-plane.

    Args:
        point1 (tuple or list): The coordinates of the first point (x1, y1, z1).
        point2 (tuple or list): The coordinates of the second point (x2, y2, z2).

    Returns:
        float: The azimuthal angle in radians between the two points and the positive x-axis in the xy-plane.
    """
    if len(point1) != 3 or len(point2) != 3:
        raise ValueError("Both points must have 3 coordinates (x, y, z)")

    x1, y1, z1 = point1
    x2, y2, z2 = point2

    # Calculate the differences in x and y coordinates
    delta_x = x2 - x1
    delta_y = y2 - y1

    # Calculate the azimuthal angle using arctan2
    azimuthal_angle = math.atan2(delta_y, delta_x)

    return azimuthal_angle

def calculate_radial_distance(point1, point2):
    """
    Calculate the radial distance (r) between two 3D coordinates.

    Args:
        point1 (tuple or list): The coordinates of the first point (x1, y1, z1).
        point2 (tuple or list): The coordinates of the second point (x2, y2, z2).

    Returns:
        float: The radial distance (r) between the two points.
    """
    if len(point1) != 3 or len(point2) != 3:
        raise ValueError("Both points must have 3 coordinates (x, y, z)")

    x1, y1, z1 = point1
    x2, y2, z2 = point2

    # Calculate the differences in x, y, and z coordinates
    delta_x = x2 - x1
    delta_y = y2 - y1
    delta_z = z2 - z1

    # Calculate the radial distance using the Euclidean distance formula
    radial_distance = math.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

    return radial_distance

def spherical_to_cartesian(starting_point, radial_distance, polar_angle, azimuthal_angle):
    """
    Convert spherical coordinates (r, θ, φ) to Cartesian coordinates (x, y, z).

    Args:
        starting_point (tuple or list): The coordinates of the starting point (x0, y0, z0).
        radial_distance (float): The radial distance (r) from the starting point.
        polar_angle (float): The polar angle (θ) in radians.
        azimuthal_angle (float): The azimuthal angle (φ) in radians.

    Returns:
        tuple: The 3D Cartesian coordinates (x, y, z).
    """
    if len(starting_point) != 3:
        raise ValueError("The starting point must have 3 coordinates (x0, y0, z0)")

    x0, y0, z0 = starting_point

    # Calculate the Cartesian coordinates
    x = x0 + radial_distance * math.sin(polar_angle) * math.cos(azimuthal_angle)
    y = y0 + radial_distance * math.sin(polar_angle) * math.sin(azimuthal_angle)
    z = z0 + radial_distance * math.cos(polar_angle)

    return (x, y, z)

def rotate_atom_position(atom_position, alpha_carbon_position, delta_theta, delta_phi):
    """
    Adjust the final position of an atom based on changes in angles between alpha carbons.

    Args:
        atom_position (tuple or list): The initial coordinates of the atom (x, y, z).
        alpha_carbon_position (tuple or list): The coordinates of the alpha carbon (x0, y0, z0).
        delta_theta (float): Change in polar angle (θ) in radians.
        delta_phi (float): Change in azimuthal angle (φ) in radians.

    Returns:
        tuple: The adjusted 3D coordinates of the atom.
    """
    if len(atom_position) != 3 or len(alpha_carbon_position) != 3:
        raise ValueError("Both atom_position and alpha_carbon_position must have 3 coordinates (x, y, z)")

    x, y, z = atom_position
    x0, y0, z0 = alpha_carbon_position

    # Translate the atom_position to a local coordinate system centered at the alpha carbon
    x_local = x - x0
    y_local = y - y0
    z_local = z - z0

    # Create the rotation matrix based on delta_theta and delta_phi
    R_theta = np.array([[1, 0, 0],
                        [0, math.cos(delta_theta), -math.sin(delta_theta)],
                        [0, math.sin(delta_theta), math.cos(delta_theta)]])

    R_phi = np.array([[math.cos(delta_phi), -math.sin(delta_phi), 0],
                      [math.sin(delta_phi), math.cos(delta_phi), 0],
                      [0, 0, 1]])

    # Rotate the local coordinates
    rotated_local_coords = np.dot(R_phi, np.dot(R_theta, np.array([x_local, y_local, z_local])))

    # Translate the rotated local coordinates back to the global coordinate system
    x_final = x0 + rotated_local_coords[0]
    y_final = y0 + rotated_local_coords[1]
    z_final = z0 + rotated_local_coords[2]

    return (x_final, y_final, z_final)

'''
New strategy to try is just generate functions to place each 
atom based on the alpha carbon location

start with FDS and work backwards based on their backbone angles.
It should work to at leat get N-CA-C--....

For now, just making it so going form CA - all atom or all atom - CA
has the correct approximate bond. 

'''
def calculate_C_current_aa(current_alpha_carbon, next_alpha_carbon, next_N, angle_degrees=111.0):
    """
    Calculate the coordinates of the C atom based on the locations of the alpha carbon (Cα) of the current amino acid,
    the alpha carbon (Cα) of the next amino acid, and the N atom of the next amino acid, while accounting for the specified
    dihedral angle (in degrees).
    """
    
    # Convert the angle from degrees to radians
    angle_radians = np.deg2rad(angle_degrees)
    
    # Calculate the vector between the current alpha carbon (Cα) and the next alpha carbon (Cα)
    vector_CA_CA_next = next_alpha_carbon - current_alpha_carbon
    
    # Normalize the vector
    vector_CA_CA_next /= np.linalg.norm(vector_CA_CA_next)
    
    # Calculate the vector between the current alpha carbon (Cα) and the N atom of the next amino acid
    vector_CA_N_next = next_N - current_alpha_carbon
    
    # Normalize the vector
    vector_CA_N_next /= np.linalg.norm(vector_CA_N_next)
    
    # Calculate the cross product between the two normalized vectors
    cross_product = np.cross(vector_CA_CA_next, vector_CA_N_next)
    
    # Normalize the cross product vector
    cross_product /= np.linalg.norm(cross_product)
    
    # Calculate the rotation matrix for the specified angle around the cross product vector
    rotation_matrix = np.array([
        [np.cos(angle_radians) + cross_product[0]**2 * (1 - np.cos(angle_radians)),
         cross_product[0] * cross_product[1] * (1 - np.cos(angle_radians)) - cross_product[2] * np.sin(angle_radians),
         cross_product[0] * cross_product[2] * (1 - np.cos(angle_radians)) + cross_product[1] * np.sin(angle_radians)],
        
        [cross_product[1] * cross_product[0] * (1 - np.cos(angle_radians)) + cross_product[2] * np.sin(angle_radians),
         np.cos(angle_radians) + cross_product[1]**2 * (1 - np.cos(angle_radians)),
         cross_product[1] * cross_product[2] * (1 - np.cos(angle_radians)) - cross_product[0] * np.sin(angle_radians)],
        
        [cross_product[2] * cross_product[0] * (1 - np.cos(angle_radians)) - cross_product[1] * np.sin(angle_radians),
         cross_product[2] * cross_product[1] * (1 - np.cos(angle_radians)) + cross_product[0] * np.sin(angle_radians),
         np.cos(angle_radians) + cross_product[2]**2 * (1 - np.cos(angle_radians))]
    ])
    
    # Calculate the coordinates of the C atom of the current amino acid using the rotation matrix
    C_coordinates = current_alpha_carbon + np.dot(rotation_matrix, vector_CA_CA_next)
    
    return C_coordinates

def calculate_next_N_coordinates(current_N, current_CA, current_C, next_alpha_carbon, angle_degrees=120.0):
    """
    Calculate the coordinates of the N atom of the next amino acid based on the locations of the N, CA, and C
    atoms of the current amino acid, and the alpha carbon (Cα) of the next amino acid, while accounting for the specified
    angle (in degrees).
    """
    
    # Convert the angle from degrees to radians
    angle_radians = np.deg2rad(angle_degrees)
    
    # Calculate the vector between the current N and CA atoms
    vector_N_CA = current_CA - current_N
    
    # Normalize the vector
    vector_N_CA /= np.linalg.norm(vector_N_CA)
    
    # Calculate the vector between the current C and CA atoms
    vector_C_CA = current_CA - current_C
    
    # Normalize the vector
    vector_C_CA /= np.linalg.norm(vector_C_CA)
    
    # Calculate the cross product between the two normalized vectors
    cross_product = np.cross(vector_C_CA, vector_N_CA)
    
    # Normalize the cross product vector
    cross_product /= np.linalg.norm(cross_product)
    
    # Calculate the rotation matrix for the specified angle around the cross product vector
    rotation_matrix = np.array([
        [np.cos(angle_radians) + cross_product[0]**2 * (1 - np.cos(angle_radians)),
         cross_product[0] * cross_product[1] * (1 - np.cos(angle_radians)) - cross_product[2] * np.sin(angle_radians),
         cross_product[0] * cross_product[2] * (1 - np.cos(angle_radians)) + cross_product[1] * np.sin(angle_radians)],
        
        [cross_product[1] * cross_product[0] * (1 - np.cos(angle_radians)) + cross_product[2] * np.sin(angle_radians),
         np.cos(angle_radians) + cross_product[1]**2 * (1 - np.cos(angle_radians)),
         cross_product[1] * cross_product[2] * (1 - np.cos(angle_radians)) - cross_product[0] * np.sin(angle_radians)],
        
        [cross_product[2] * cross_product[0] * (1 - np.cos(angle_radians)) - cross_product[1] * np.sin(angle_radians),
         cross_product[2] * cross_product[1] * (1 - np.cos(angle_radians)) + cross_product[0] * np.sin(angle_radians),
         np.cos(angle_radians) + cross_product[2]**2 * (1 - np.cos(angle_radians))]
    ])
    
    # Calculate the coordinates of the N atom of the next amino acid using the rotation matrix
    next_N_coordinates = next_alpha_carbon + np.dot(rotation_matrix, vector_N_CA)
    
    return next_N_coordinates


def add_necessary_C_N(PDBParserObj):
    '''
    function to add in N or C between CA and all atom
    to make bonds not ... off.
    '''
    # get the all atom coords
    all_coords = PDBParserObj.all_atom_coords_by_index
    # get the region dict
    region_info = PDBParserObj.regions_dict
    #get the idrs
    IDR_regions = [region for region in list(region_info.keys()) if 'idr' in region]
    # get the IDR indices
    IDR_residue_indices=[]
    for region in IDR_regions:
        cur_region=region_info[region]
        IDR_residue_indices.extend([aa for aa in range(cur_region[0], cur_region[1])])
    # add loops
    loop_regions = [region for region in list(region_info.keys()) if 'loop' in region]
    if loop_regions != []:
        for region in loop_regions:
            cur_region = region_info[region]
            cur_loop = cur_region[1]
            for loops in cur_loop:
                # add loops and connect back.
                IDR_residue_indices.extend([aa for aa in range(loops[0]-1, loops[1]+1)])

    # add in the Cs
    for aa in range(0, len(all_coords)-2):
        if aa != 0:
            if aa in IDR_residue_indices and aa+1 not in IDR_residue_indices:
                #calculate_C_current_aa(current_alpha_carbon, next_alpha_carbon, next_N
                next_values = all_coords[aa+1]
                next_n=np.array(next_values['N'])
                next_ca=np.array(next_values['CA'])
                current_ca=np.array(all_coords[aa]['CA'])
                add_C = calculate_C_current_aa(current_ca, next_ca, next_n)
                PDBParserObj.all_atom_coords_by_index[aa]['CA']=current_ca
                PDBParserObj.all_atom_coords_by_index[aa]['C']=add_C
            if aa-1 not in IDR_residue_indices and aa in IDR_residue_indices:
                current_values = all_coords[aa-1]
                current_n=np.array(current_values['N'])
                current_ca=np.array(current_values['CA'])
                current_c=np.array(current_values['C'])
                next_ca=np.array(all_coords[aa]['CA'])
                add_N = calculate_next_N_coordinates(current_n, current_ca, current_c, next_ca)
                PDBParserObj.all_atom_coords_by_index[aa]={}
                PDBParserObj.all_atom_coords_by_index[aa]['N']=np.array(add_N)
                PDBParserObj.all_atom_coords_by_index[aa]['CA']=np.array(next_ca)
    return PDBParserObj


