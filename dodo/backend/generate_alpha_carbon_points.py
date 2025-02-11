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


