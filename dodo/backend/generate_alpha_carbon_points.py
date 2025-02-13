'''
Code to generate positions for alpha carbons.
The approach is to basically use points on a circle formed by 
a cone with the apex at the coordinates for the second alpha carbon
and the base defined by an angle that I've determined based on AF2
structures.

This generates many possible points that reasonably could be alpha carbon positions
based on two previous alpha carbons. 

The main reason for this over a random walk is it lets us avoid truly
weird alpha carbon positions that would be impossible in reality. Downstream,
I'm hoping this makes rebuilding the backbone easier.
'''

import numpy as np
from dodo.backend.dodo_math import find_transform, reverse_transform, generate_circle_points


def generate_ca_points(ca1, ca2, points_per_angle=10, num_angles=None, 
                       min_angle=91, max_angle=161, angles=None, ideal_angle=125):
    '''
    Note: In IDRs, this can vary form around 90 to around ~160. Based on 
    AF2 structures, I've seen angles with mean at around 125 and range being 75 to 179 with
    a standard deviation of ~26. Going to do mean +/- std dev for now.'''
    ca1 = np.asarray(ca1)
    ca2 = np.asarray(ca2)
    translation, rotation = find_transform(ca1, ca2)

    if num_angles is None:
        num_angles = max_angle - min_angle + 1
    
    if angles is None:
        angles = np.linspace(min_angle, max_angle, num_angles)
    
    # sort angles by how close they are to the ideal angle
    angles = angles[np.argsort(np.abs(angles - ideal_angle))]
    
    # Pre-compute constants
    angles_rad = np.deg2rad(180 - angles)
    cos_angles = np.cos(angles_rad)
    sin_angles = np.sin(angles_rad)
    
    # Constants
    CA_DISTANCE = 3.856
    
    # Vectorized calculations
    h_values = CA_DISTANCE * cos_angles
    r_values = CA_DISTANCE * sin_angles
    
    ca_points = []
    for h, r in zip(h_values, r_values):
        points = generate_circle_points(
            [CA_DISTANCE, 0, 0],
            [CA_DISTANCE + h, 0, 0],
            r,
            points_per_angle
        )
        # Vectorized transform for all points
        transformed_points = np.array([reverse_transform(p, translation, rotation) for p in points])
        ca_points.extend(transformed_points)
    
    return np.array(ca_points)


def generate_ca_points_vectorized(ca_pairs, points_per_angle=10, num_angles=None,
                                min_angle=91, max_angle=161, angles=None, ideal_angle=125):
    '''
    Vectorized version that generates CA points for multiple pairs simultaneously.
    
    Parameters
    ----------
    ca_pairs : np.array, shape=(N,2,3)
        Array of N pairs of CA coordinates, where each pair contains two 3D points
    points_per_angle : int
        Number of points to generate per angle
    num_angles : int, optional
        Number of angles to use between min_angle and max_angle
    min_angle : float
        Minimum angle in degrees
    max_angle : float
        Maximum angle in degrees
    angles : array-like, optional
        Specific angles to use (in degrees)
    ideal_angle : float
        The ideal angle to sort generated angles by

    Returns
    -------
    np.array, shape=(N,K,3)
        Array of K possible next CA positions for each of the N input pairs
    '''
    ca_pairs = np.asarray(ca_pairs)
    if ca_pairs.ndim != 3 or ca_pairs.shape[1:] != (2, 3):
        raise ValueError("ca_pairs must have shape (N,2,3)")
    
    # get angles
    if angles is None:
        if num_angles is not None:
            if (max_angle - min_angle) + 1 < num_angles:
                all_angle_num = num_angles + 1
            else:
                all_angle_num = max_angle - min_angle + 1
        else:
            all_angle_num = max_angle - min_angle + 1
        all_angles = np.linspace(min_angle, max_angle, all_angle_num)
        if num_angles is not None:
            if num_angles < all_angle_num:
                angles = np.random.choice(all_angles, num_angles, replace=False)
        else:
            angles = all_angles

    # account for if angles is a single angle.
    if isinstance(angles, (int, float)):
        angles = np.array([angles])

    
    
    # Constants
    CA_DISTANCE = 3.856
    angles_rad = np.deg2rad(180 - angles)
    cos_angles = np.cos(angles_rad)
    sin_angles = np.sin(angles_rad)
    
    # Calculate h and r values for all angles
    h_values = CA_DISTANCE * cos_angles
    r_values = CA_DISTANCE * sin_angles
    
    # Get transformations for all pairs
    translations = -ca_pairs[:, 0]  # Shape: (N,3)
    # Fix: Remove np.newaxis to prevent extra dimension
    ca2_translated = ca_pairs[:, 1] + translations  # Shape: (N,3)
    
    # Calculate rotations for all pairs
    
    # Normalize translated vectors
    norms = np.linalg.norm(ca2_translated, axis=1, keepdims=True)
    v = ca2_translated / norms
    x_axis = np.array([1., 0., 0.])
    
    # Calculate cos_theta values
    cos_theta = np.dot(v, x_axis)
    
    # Initialize rotations as identity matrices
    rotations = np.tile(np.eye(3), (len(ca_pairs), 1, 1))
    
    # Calculate rotation axes and angles for pairs that need rotation
    mask = ~np.isclose(cos_theta, 1.0)
    
    rotation_axis = np.cross(v[mask], x_axis)
    rotation_axis_norm = np.linalg.norm(rotation_axis, axis=1, keepdims=True)
    
    # Normalize rotation axes
    rotation_axis_normalized = rotation_axis / rotation_axis_norm
    
    sin_theta = np.sqrt(1 - cos_theta[mask]**2)
    
    # Construct K matrices
    K = np.zeros((len(rotation_axis), 3, 3))
    K[:, 0, 1] = -rotation_axis_normalized[:, 2]
    K[:, 0, 2] = rotation_axis_normalized[:, 1]
    K[:, 1, 0] = rotation_axis_normalized[:, 2]
    K[:, 1, 2] = -rotation_axis_normalized[:, 0]
    K[:, 2, 0] = -rotation_axis_normalized[:, 1]
    K[:, 2, 1] = rotation_axis_normalized[:, 0]
    
    # Calculate Rodrigues rotation matrices
    identity = np.tile(np.eye(3), (len(rotation_axis), 1, 1))
    rotations[mask] = identity + sin_theta[:, None, None] * K + (1 - cos_theta[mask])[:, None, None] * np.matmul(K, K)
    
    # Transform points back
    all_transformed_points = []
    for h, r in zip(h_values, r_values):
        circle_points = generate_circle_points(
            [CA_DISTANCE, 0, 0],
            [CA_DISTANCE + h, 0, 0],
            r,
            points_per_angle
        )
        transformed_points = np.matmul(circle_points, np.transpose(rotations, (0, 2, 1))) - translations[:, None, :]
        all_transformed_points.append(transformed_points)
    transformed_points = np.concatenate(all_transformed_points, axis=1)
    
    return transformed_points


