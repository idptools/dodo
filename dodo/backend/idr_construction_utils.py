'''
various functions we want to use to
construct the IDR. Specifically, this functionality
will focus on taking the candidate coordinates from
the generate_alpha_carbon_points.py function and determine
which one is the best to use for the next alpha carbon
position
'''

import numpy as np

def find_points_within_sphere(points_of_interest, 
                              sphere_center, 
                              sphere_radius):
    '''
    This function checks if a point is within a sphere.
    
    Parameters
    ----------
    points_of_interest : np.array, shape=(N, 3)
        points we want to check if they are within the sphere
    sphere_center : np.array, shape=(3,)
        center of the sphere
    sphere_radius : float
        radius of the sphere
    
    Returns
    -------
    bool
        True if the point is within the sphere, False otherwise
    '''
    # calculate distances
    distances = np.linalg.norm(points_of_interest - sphere_center, axis=1)

    # filter points to only those within the sphere
    within_sphere = distances <= sphere_radius
    return within_sphere

def find_point_closest_to_sphere_surface(points_of_interest, 
                                         sphere_center, 
                                         sphere_radius):
    '''
    This function finds the point closest to the surface of a sphere.
    
    Parameters
    ----------
    points_of_interest : np.array, shape=(N, 3)
        points we want to check if they are within the sphere
    sphere_center : np.array, shape=(3,)
        center of the sphere
    sphere_radius : float
        radius of the sphere
    
    Returns
    -------
    np.array, shape=(3,)
        point closest to the sphere surface
    '''
    # calculate distances
    distances = np.linalg.norm(points_of_interest - sphere_center, axis=1)

    # find the point closest to the sphere surface
    closest_point = points_of_interest[np.argmin(np.abs(distances - sphere_radius))]
    return closest_point


    
