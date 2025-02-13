'''
various functions we want to use to
construct the IDR. 

New strategy is going to be to build a lot of IDRs and then then position necessary 
FDs to get the right distances and then see which IDRs do not clash with the rest of the chain.
This will be a lot faster than the old strategy of building the IDRs and then checking the distances.
'''

import random

import numpy as np
import time

from dodo.backend.dodo_structures import Atom, Monomer, Domain, Chain, Complex
from dodo.backend.generate_alpha_carbon_points import generate_ca_points_vectorized
import dodo_math
from dodo.backend import utils


def build_candidate_idrs(num_residues, num_idrs, idr_expansion='standard',
                         start_coord=None, end_coord=None, 
                         points_to_check_clashing=None,
                         end_point_dist_range=(3.7, 4.0),
                         batch_size=50,
                         max_attempts=10):
    """
    build a bunch of candidate IDRs using the generate_alpha_carbon_points
    functionality. Uses vectorized numpy magic to keep things fast.
    
    Parameters
    ----------
    num_residues : int
        number of residues in the IDR
    num_idrs : int
        number of IDRs to generate
    idr_expansion : str
        How expanded we want the IDR to be. Options are 'compact',
        'standard', 'expanded'. Default is 'standard'.
    start_coord : np.ndarray
        3D coordinate of the start of the IDR. If None, will use (0,0,0)
    end_coord : np.ndarray
        3D coordinate of the end of the IDR. 
        If None, will use a random point based on IDR expansion and num_residues
    points_to_check_clashing : np.ndarray
        3D coordinates of points to check for clashes. Shape should be (M, 3)
    end_point_dist_range : tuple
        range of distances allowed for within_sphere to cur_sphere_center at the end of the IDR
        default is (3.7, 4.0)
    batch_size : int
        number of IDRs to try generating in each batch
    max_attempts : int
        maximum number of attempts to generate the requested number of IDRs

    Returns
    -------
    np.ndarray
        array of candidate IDRs. Shape is (num_idrs, num_residues, 3)    
    """
    if idr_expansion != 'standard':
        raise NotImplementedError('Only standard IDR expansion is supported for now.')

    # set total distance for end-to-end
    total_distance = 1.4 * num_residues
    
    # Initialize array to store successful IDRs
    successful_idrs = []
    attempts = 0

    while len(successful_idrs) < num_idrs and attempts < max_attempts:
        # Initialize the starting coordinates
        if start_coord is None:
            start_coord = np.zeros(3)
        start_coords = np.tile(start_coord, (batch_size, 1))

        # Initialize the ending coordinates or generate them
        if end_coord is not None:
            end_coords = np.tile(end_coord, (batch_size, 1))
        else:
            # Generate random end points at total_distance
            end_coords = dodo_math.get_ending_points(
                start_coords,
                np.full(batch_size, total_distance),
                clash_points=points_to_check_clashing,
                num_points_to_check=1000
            )

        # Initialize array to store all coordinates for this batch
        all_coords = np.zeros((batch_size, num_residues, 3))
        all_coords[:, 0] = start_coords
        all_coords[:, -1] = end_coords

        # Initialize separate clash points array for each IDR
        if points_to_check_clashing is None:
            clash_points_per_idr = [np.array([start_coords[i], end_coords[i]]) for i in range(batch_size)]
        else:
            clash_points_per_idr = [np.vstack([points_to_check_clashing, start_coords[i], end_coords[i]]) 
                                   for i in range(batch_size)]
            
        # now get radii for sphere
        start_radius = total_distance
        position_mod = [x for x in [15.8868, 13.6407, 11.0914, 8.6688, 6.441, 3.89] if x < start_radius]
        starting_sphere_radii = np.linspace(start_radius, position_mod[0], (num_residues-1)-len(position_mod))
        sphere_radii = np.concatenate((starting_sphere_radii, position_mod))

        # Generate intermediate points
        for i in range(1, num_residues-1):
            # Get the previous points
            prev_points = all_coords[:, i-1]
            
            # Generate multiple candidate points
            potential_points = dodo_math.generate_random_points_by_distance_mult(
                prev_points, 
                np.full(batch_size, 3.8),  # CA-CA distance
                1000  # number of candidates per point
            )

            # Process each IDR separately
            valid_idrs = []
            for j in range(len(all_coords)):
                # Check for clashes
                non_clashing = dodo_math.find_points_not_clashing(
                    potential_points[j:j+1],
                    clash_points_per_idr[j],
                    3.0
                )

                if len(non_clashing) == 0:
                    continue

                # Get current sphere center and radius
                cur_sphere_center = end_coords[j]
                cur_sphere_radius = sphere_radii[i]

                # Find points within acceptable distance to previous point
                distances_to_prev = np.linalg.norm(non_clashing - prev_points[j], axis=1)
                valid_distance_mask = np.abs(distances_to_prev - 3.8) < 0.1
                valid_points = non_clashing[valid_distance_mask]

                if len(valid_points) == 0:
                    continue

                best_point = dodo_math.find_point_closest_to_sphere_surface(
                    valid_points, 
                    cur_sphere_center, 
                    cur_sphere_radius
                )

                if i > 1:
                    dist_to_prev = np.linalg.norm(best_point - all_coords[j, i-1])
                    if abs(dist_to_prev - 3.8) > 0.1:
                        continue

                if i == num_residues-2:
                    end_point_dist = np.linalg.norm(best_point - cur_sphere_center)
                    if end_point_dist < end_point_dist_range[0] or end_point_dist > end_point_dist_range[1]:
                        continue

                all_coords[j, i] = best_point
                clash_points_per_idr[j] = np.vstack([clash_points_per_idr[j], best_point])
                valid_idrs.append(j)

            # Keep only the IDRs that were successfully built up to this point
            if len(valid_idrs) == 0:
                break
                
            all_coords = all_coords[valid_idrs]
            clash_points_per_idr = [clash_points_per_idr[j] for j in valid_idrs]
            end_coords = end_coords[valid_idrs]

        # Add successful IDRs from this batch to our collection
        successful_idrs.extend([coords for coords in all_coords])
        attempts += 1

        # If we have enough IDRs, we're done
        if len(successful_idrs) >= num_idrs:
            break

    if len(successful_idrs) < num_idrs:
        raise Exception(f'Only generated {len(successful_idrs)} valid IDRs after {max_attempts} attempts')

    # Convert to numpy array and return only the requested number of IDRs
    return np.array(successful_idrs[:num_idrs])
