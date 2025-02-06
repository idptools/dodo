'''
Code to generate positions for alpha carbons.
The approach is to basically use points on a circle formed by 
a cone with the apex at the coordinates for the second alpha carbon
and the base defined by an angle that I've determined based on AF2
structures. Translation and rotation are used to make the math a little easier
and then that is applied to the final point to reverse the original transformation.

This generates many possible points that reasonably could be alpha carbon positions
based on two previous alpha carbons. The idea here is that I can then check for clashing
and distance constraints and choose the next alpha carbon based on that. This will let me
have less weird angles in my generated alpha carbons which should make it easier to build the
backbone of the protein down the line. It's still pretty fast.
'''

import numpy as np

def find_transform(coord1, coord2):
    p1 = np.asarray(coord1)
    p2 = np.asarray(coord2)
    
    translation = -p1
    p2_translated = p2 + translation
    
    # Calculate rotation matrix to align p2 with x-axis
    v = p2_translated / np.linalg.norm(p2_translated)
    x_axis = np.array([1., 0., 0.])
    
    # Optimize rotation calculation
    cos_theta = np.dot(v, x_axis)
    if np.isclose(cos_theta, 1.0):
        return translation, np.eye(3)
    if np.isclose(cos_theta, -1.0):
        return translation, -np.eye(3)
    
    rotation_axis = np.cross(v, x_axis)
    rotation_axis /= np.linalg.norm(rotation_axis)
    
    # Pre-compute values for Rodrigues formula
    sin_theta = np.sqrt(1 - cos_theta**2)
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                  [rotation_axis[2], 0, -rotation_axis[0]],
                  [-rotation_axis[1], rotation_axis[0], 0]])
    
    rotation = np.eye(3) + sin_theta * K + (1 - cos_theta) * (K @ K)
    return translation, rotation

def apply_transform(coord, translation, rotation):
    return rotation @ (np.asarray(coord) + translation)

def reverse_transform(coord, translation, rotation):
    return (rotation.T @ np.asarray(coord)) - translation

def generate_circle_points(origin, center, radius, num_points):
    origin = np.asarray(origin)
    center = np.asarray(center)
    
    normal = center - origin
    normal /= np.linalg.norm(normal)
    
    # Pre-compute angles for all points at once
    angles = np.linspace(0, 2*np.pi, num_points, endpoint=False)
    
    # Optimize orthogonal vector calculation
    v1 = np.array([1., 0., 0.]) if not np.allclose(normal, [1., 0., 0.]) else np.array([0., 1., 0.])
    u = np.cross(normal, v1)
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    
    # Vectorized point generation
    cos_theta = np.cos(angles)
    sin_theta = np.sin(angles)
    
    # Generate all points at once using broadcasting
    points = center + radius * (u[:, np.newaxis] * cos_theta + v[:, np.newaxis] * sin_theta).T
    return points

def generate_ca_points(ca1, ca2, points_per_angle=10, num_angles=10, 
                       min_angle=99, max_angle=151, angles=None, ideal_angle=125):
    '''
    Note: In IDRs, this can vary form around 90 to around ~160. Based on 
    AF2 structures, I've seen angles with mean at around 125 and range being 75 to 179 with
    a standard deviation of ~26. Going to do mean +/- std dev for now.'''
    ca1 = np.asarray(ca1)
    ca2 = np.asarray(ca2)
    translation, rotation = find_transform(ca1, ca2)
    
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

