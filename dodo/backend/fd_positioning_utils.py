'''
code for determining the translation vector and rotation matrix
necessary to translate and rotate a folded domain such that a specified 
coordinate on that folded domain is facing towards a specified coordinate
on another folded domain and the two coordinates are a specified distance apart. 

Code could probably be refactored but going to leave it alone for now.
'''

import numpy as np

from dodo.backend.dodo_math import get_random_coordinates_on_sphere_surface

def get_rot_matrix_fd2_coord_towards_fd1(FD2_coord, FD2_COM, FD1_coord):
    '''
    Function to rotate FD2 to face FD1_coord while keeping
    FD2_COM in the same location. This function will return
    the rotation matrix that can be used to rotate FD2.
    '''
     # now we need to rotate FD2 such that it's COM stays in the same location
    # but FD2_coord is facing towards FD1_coord
    # Get current and target vectors
    current_vector = FD2_coord - FD2_COM
    target_vector = FD1_coord - FD2_coord

    # Normalize vectors
    current_vector = current_vector / np.linalg.norm(current_vector)
    target_vector = target_vector / np.linalg.norm(target_vector)

    # Calculate rotation matrix using Rodrigues' formula
    rotation_axis = np.cross(current_vector, target_vector)
    if np.all(np.abs(rotation_axis) < 1e-10):
        if np.dot(current_vector, target_vector) < 0:
            # Vectors are antiparallel, rotate 180° around any perpendicular axis
            rotation_axis = np.array([1, 0, 0]) if not np.allclose(current_vector, [1, 0, 0]) else np.array([0, 1, 0])
            angle = np.pi
        else:
            # Vectors are parallel, no rotation needed
            return None
    else:
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        cos_angle = np.clip(np.dot(current_vector, target_vector), -1.0, 1.0)
        angle = np.arccos(cos_angle)

    # Create rotation matrix using Rodrigues' formula
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                  [rotation_axis[2], 0, -rotation_axis[0]],
                  [-rotation_axis[1], rotation_axis[0], 0]])
    
    rotation_matrix = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

    return rotation_matrix


def position_folded_domain(FD1_index, FD2_index, chain_obj,
                         objective_final_distance):
    '''
    This function translates and rotates FD2 such that
    a line can be formed from FD1_center_of_mass through
    FD1_coord and FD2_center_of_mass and FD2_coord. This
    line will be the line that FD2_coord will be facing
    towards. FD2 will be translated and rotated such that
    FD2_coord is a specified distance from FD1_coord.

    Parameters
    ----------
    FD1_index : index for FD1 in the chain object
    FD2_index : index for FD2 in the chain object
    chain_obj : Chain object that is associated with the FDs
    objective_final_distance : float
        The distance between FD1_coord and FD2_coord

    Returns
    -------
    FD2 : Domain object
        The translated and rotated FD2
    '''
    # get FD1 and FD2
    FD1 = chain_obj.domains[FD1_index]
    FD2 = chain_obj.domains[FD2_index]

    # now we need to get the COM for FD1 and FD2
    FD2_COM = FD2.identify_center_of_mass()

    # get the coord that needs to be facing towards the 
    # other FD and be a specified distance away
    FD1_coord = FD1.get_coordinates_array()[-1]
    FD2_coord = FD2.get_coordinates_array()[0]

    # Try random positions on sphere surface until we find one without clashes
    potential_coordinates = get_random_coordinates_on_sphere_surface(FD1_coord, objective_final_distance, 500)

    success = False
    for target_position in potential_coordinates:
        # Calculate translation vector to move FD2_coord to target position
        translation_vector = target_position - FD2_coord
        # Try the translation
        FD2.translate(translation_vector)
        
        # this will return False if no clashes, otherwise will return the index of the clashing domain
        domain_clashes = chain_obj.check_if_domain_clashes(FD2_index)
        if domain_clashes == False:
            success = True
            break
        
        # Undo translation if it caused clashes
        FD2.translate(-translation_vector)
    
    # Raise error if no clash-free position was found
    if not success:
        raise ValueError("Could not find clash-free position for FD2")

    # get rot matrix for rotation
    rotation_matrix = get_rot_matrix_fd2_coord_towards_fd1(FD2_coord, FD2_COM, FD1_coord)

    # if rot matrix not none, apply rotation
    if rotation_matrix is not None:
        # Apply rotation
     FD2.rotate(rotation_matrix)

    # finally do the last translation of FD2 to make sure that FD2_coord
    # is the specified distance from FD1_coord
    FD2_coord = FD2.get_coordinates_array()[0]
    distance = np.linalg.norm(FD1_coord - FD2_coord)
    translation_vector = FD1_coord + (objective_final_distance/distance)*(FD2_coord-FD1_coord) - FD2_coord
    FD2.translate(translation_vector)

    return FD2


def position_folded_domain_linear(FD1_index, FD2_index, chain_obj,
                         objective_final_distance):
    '''
    This function translates and rotates FD2 such that
    a line can be formed from FD1_center_of_mass through
    FD1_coord and FD2_center_of_mass and FD2_coord. This
    line will be the line that FD2_coord will be facing
    towards. FD2 will be translated and rotated such that
    FD2_coord is a specified distance from FD1_coord.

    Like position_folded_domain, but positions everything linearly.

    Parameters
    ----------
    FD1_index : index for FD1 in the chain object
    FD2_index : index for FD2 in the chain object
    chain_obj : Chain object that is associated with the FDs
    objective_final_distance : float
        The distance between FD1_coord and FD2_coord

    Returns
    -------
    FD2 : Domain object
        The translated and rotated FD2
    '''
    # get FD1 and FD2
    FD1 = chain_obj.domains[FD1_index]
    FD2 = chain_obj.domains[FD2_index]

    # translate FD2 such that it's COM is at 0,0,0.
    FD2.translate(-FD2.identify_center_of_mass())

    # Translate FD2 along x-axis by the objective distance
    # with respect to FD1_coord
    FD1_coord = FD1.get_coordinates_array()[-1]
    FD2_coord = FD2.get_coordinates_array()[0]
    translation_vector = FD1_coord + np.array([objective_final_distance, 0, 0]) - FD2_coord
    FD2.translate(translation_vector)
    
    # get new FD2 coord, FD2_COM, and FD1_coord
    FD2_coord = FD2.get_coordinates_array()[0]
    FD2_COM = FD2.identify_center_of_mass()
    
    # get rot matrix for rotation
    rotation_matrix = get_rot_matrix_fd2_coord_towards_fd1(FD2_coord, FD2_COM, FD1_coord)

    # if rot matrix not none, apply rotation
    if rotation_matrix is not None:
        # Apply rotation
        FD2.rotate(rotation_matrix)

    # now translate FD2 such that FD2_coord is the specified distance from FD1_coord
    FD2_coord = FD2.get_coordinates_array()[0]
    distance = np.linalg.norm(FD1_coord - FD2_coord)
    translation_vector = FD1_coord + (objective_final_distance/distance)*(FD2_coord-FD1_coord) - FD2_coord
    FD2.translate(translation_vector)

    return FD2
