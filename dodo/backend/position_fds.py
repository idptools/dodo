'''
code for determining the translation vector and rotation matrix
necessary to translate and rotate a folded domain such that a specified 
coordinate on that folded domain is facing towards a specified coordinate
on another folded domain and the two coordinates are a specified distance apart. 
'''

import numpy as np

def position_folded_domain(FD1_coord, FD1_center_of_mass, 
                         FD2_coord, FD2_center_of_mass, 
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
    FD1_coord : np.array, shape=(3,)
        coordinate that we want to face towards FD2_coord
    FD1_center_of_mass : np.array, shape=(3,)
        center of mass of folded domain 1
    FD2_coord : np.array, shape=(3,)
        coordinate that we want to face towards FD1_coord
    FD2_center_of_mass : np.array, shape=(3,)
        center of mass of folded domain 2
    objective_final_distance : float
        objective final distance between FD1_coord and FD2_coord

    Returns
    -------
    rotation : np.array, shape=(3, 3)
        rotation matrix necessary to apply to fd2 to get it in the correct position
    translation : np.array, shape=(3,)
        translation vector necessary to apply to fd2 to get it in the correct position
    '''
    # Calculate vectors
    v1 = FD1_coord - FD1_center_of_mass
    v2 = FD2_coord - FD2_center_of_mass

    # Normalize vectors
    v1_norm = v1 / np.linalg.norm(v1)
    v2_norm = v2 / np.linalg.norm(v2)

    # Calculate rotation axis and angle
    rotation_axis = np.cross(v2_norm, -v1_norm)
    if np.all(rotation_axis == 0):
        # Vectors are parallel, rotation axis is arbitrary perpendicular vector
        rotation_axis = np.array([1, 0, 0]) if np.allclose(v1_norm, [0, 1, 0]) else np.cross(v1_norm, [0, 1, 0])
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    cos_angle = np.dot(v2_norm, -v1_norm)
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))

    # Create rotation matrix using Rodrigues' rotation formula
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                  [rotation_axis[2], 0, -rotation_axis[0]],
                  [-rotation_axis[1], rotation_axis[0], 0]])
    rotation = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

    # Calculate the position where FD2_coord should end up
    target_position = FD1_coord + v1_norm * objective_final_distance

    # Calculate where FD2_coord will be after rotation (relative to center of mass)
    rotated_relative_pos = np.dot(rotation, FD2_coord - FD2_center_of_mass)
    
    # Calculate translation vector that moves the rotated coordinate to target
    translation = target_position - (rotated_relative_pos + FD2_center_of_mass)

    return rotation, translation


def position_fds_along_x_axis(FD1_coord, FD1_center_of_mass, 
                              FD2_coord, FD2_center_of_mass, 
                              objective_final_distance):
    '''
    This function translates and rotates FD2 such that
    a line can be formed from FD1_center_of_mass through
    FD1_coord and FD2_center_of_mass and FD2_coord. This
    line will be the line that FD2_coord will be facing
    towards. FD2 will be translated and rotated such that
    FD2_coord is a specified distance from FD1_coord. FD2
    will be positioned such that FD1_coord and FD2_coord
    are on the x-axis.

    Parameters
    ----------
    FD1_coord : np.array, shape=(3,)
        coordinate that we want to face towards FD2_coord
    FD1_center_of_mass : np.array,
        center of mass of folded domain 1
    FD2_coord : np.array, shape=(3,)
        coordinate that we want to face towards FD1_coord
    FD2_center_of_mass : np.array, shape=(3,)
        center of mass of folded domain 2

    Returns
    -------
    rotation : np.array, shape=(3, 3)
        rotation matrix necessary to apply to fd2 to get it in the correct position
    translation : np.array, shape=(3,)
        translation vector necessary to apply to fd2 to get it in the correct position
    '''
    # Calculate vectors
    # Calculate vectors
    v1 = FD1_coord - FD1_center_of_mass
    v2 = FD2_coord - FD2_center_of_mass

    # Get initial rotation and translation
    rotation1, translation1 = position_folded_domain(FD1_coord, FD1_center_of_mass, 
                                                   FD2_coord, FD2_center_of_mass, 
                                                   objective_final_distance)

    # Apply initial rotation and translation to FD2_coord
    rotated_FD2_coord = np.dot(rotation1, FD2_coord - FD2_center_of_mass) + translation1 + FD2_center_of_mass

    # Calculate vector to rotate onto x-axis
    vector = rotated_FD2_coord - FD1_coord
    vector_norm = vector / np.linalg.norm(vector)
    x_axis = np.array([1, 0, 0])

    # Calculate rotation to align with x-axis
    rotation_axis = np.cross(vector_norm, x_axis)
    if np.all(rotation_axis == 0):
        rotation2 = np.eye(3)
    else:
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        cos_angle = np.dot(vector_norm, x_axis)
        angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
        K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                      [rotation_axis[2], 0, -rotation_axis[0]],
                      [-rotation_axis[1], rotation_axis[0], 0]])
        rotation2 = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

    # Combine rotations and translation
    rotation = np.dot(rotation2, rotation1)
    translation = translation1

    return rotation, translation
