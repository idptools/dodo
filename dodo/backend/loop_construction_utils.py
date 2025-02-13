'''
Code for building the loops. 

The primary function here takes in the chain 
as well as the index for the loop domain and then
from there builds all loops for that domain. 
'''

import numpy as np
import os


# import dodo math stuff
from dodo.backend.dodo_math import find_points_within_sphere, find_point_closest_to_sphere_surface, calculate_distance, find_furthest_coordinate, find_points_not_clashing
# import CA generation code.
from dodo.backend.generate_alpha_carbon_points import generate_ca_points
# import utility functions
from dodo.backend.utils import write_dumby_pdb


# DELETE THESE IMPORTS LATER
# import reader, delete this later.
from dodo_readers import Reader
from dodo_structures import Atom, Monomer, Domain, Chain, Complex
from dodo.backend.domain_identification import get_fds_loops_idrs, assign_domains_from_dict




def build_loops(chain, loop_domain_index, num_conformations=1, num_attempts=10):
    '''
    This function builds the loops for a given domain. 
    It takes in the chain, the index of the loop domain,
    and the number of conformations to generate.
    
    Parameters
    ----------
    chain : Chain
        Chain object that we are building loops for
    loop_domain_index : int
        index of the loop domain
    num_conformations : int
        number of conformations to generate for the loops
        
    Returns
    -------
    None
    '''
    # get the loop domain
    loop_domain = chain.domains[loop_domain_index]
    # get the coords for the loop domain by index of the monomers
    coords_by_index = loop_domain.get_coords_by_monomer_index_with_atom_dict()
    # get the loops
    all_loop_indices = loop_domain.get_loop_indices()
    # this should never happen. 
    if all_loop_indices == []:
        raise ValueError('No loops found for this domain! Check why this is happening. Please raise an issue on Github!')

    # get the coordinates for everything currently rebuilt. 
    # this will be used to check for clashes.
    all_coords_to_check = chain.get_coords_of_rebuilt_domains()

    # now we will be building things to be within a sphere of a known distance.
    # this distance is equal to the distance between the CA atoms of the two
    # monomers that are on either side of the loop.
    # thus we can eliminate coords not within this distance + 2 angsrtoms.

    # now iterate over the loops, which are in a nested list format.
    for loop_start_end in all_loop_indices:
        # get loop indices
        loop_indices = [i for i in range(loop_start_end[0], loop_start_end[1]+1)]

        # add empty vlues to coords_by_index to populate later. 
        for index in range(loop_indices[0], loop_indices[-1]+1):
            coords_by_index[index] = {'CA': None}
            

        # make a copy of all_coords_to_check so we don't mess up the original
        coords_to_check = all_coords_to_check.copy()

        start_loop_index = loop_indices[0]
        end_loop_index = loop_indices[-1]
        
        # now we build coords by being within a sphere that shrinks iteratively.
        # The distance will start the distance from the two coords we just generated.
        # it will then shrink from there.

        start_radius = calculate_distance(np.array(coords_by_index[start_loop_index-1]['CA']), 
                                                 np.array(coords_by_index[end_loop_index+1]['CA']))
        
        # remove points outside of sphere from coords to check
        coords_to_check = find_points_within_sphere(coords_to_check, np.array(coords_by_index[end_loop_index+1]['CA']), start_radius+2)
        
        # now use np.linspace to get the radii of the spheres from start + 3.8 to 3.8
        position_mod = [x for x in [15.8868, 13.6407, 11.0914, 8.6688, 6.441, 3.89] if x < start_radius]
        starting_sphere_radii = np.linspace(start_radius, position_mod[0], len(loop_indices)-len(position_mod))
        sphere_radii = np.concatenate((starting_sphere_radii, position_mod))

        # set the sphere center as the end point we are building towards.
        sphere_center = coords_by_index[end_loop_index+1]['CA']

        # now we iterate over the loop indices and build the coords.
        for num, index in enumerate(loop_indices):
            # get cur sphere radius
            cur_sphere_radius=sphere_radii[num]
            # get the CA for the previous 2 monomers
            previous_2_CA = np.array(coords_by_index[index-2]['CA'])
            previous_CA = np.array(coords_by_index[index-1]['CA'])


            # get possible alpha carbon positions
            possible_next_CA_coords = generate_ca_points(previous_2_CA, previous_CA)
            # get points within the sphere    
            potential_coords = find_points_within_sphere(possible_next_CA_coords, sphere_center, cur_sphere_radius)

            # get coords not clashing with coords_to_check
            not_clashing = find_points_not_clashing(potential_coords, coords_to_check, 
                                                        clash_distance=3)

            if len(not_clashing)==0:
                not_clashing = find_points_not_clashing(potential_coords, coords_to_check, 
                                                        clash_distance=2.5)
            if len(not_clashing)==0:
                not_clashing = find_points_not_clashing(potential_coords, coords_to_check, 
                                                        clash_distance=2.0)
            
            if len(not_clashing)==0:
                raise ValueError('No non-clashing points found for loop monomer!')

            # if this isn't the final coord we are building, choose a random coord.
            chosen_coord = find_point_closest_to_sphere_surface(potential_coords, sphere_center, cur_sphere_radius)

            # add the coord to coords_by_index
            coords_by_index[index]['CA'] = chosen_coord
            coords_to_check = np.append(coords_to_check, [chosen_coord], axis=0)

        # now we repopulate the monomers in the chain.
        # for index in loop_indices:
        # first get the last atom in the previous monomer
        last_monomer = loop_domain.monomers[loop_indices[0]-1]
        last_atoms = last_monomer.atoms
        next_atom_id = last_atoms[list(last_atoms.keys())[-1]].atom_id

        loop_ind_num_to_monomer_id = loop_domain.monomer_ind_to_aa
        for index in loop_indices:
            next_atom_id+=1
            # get the CA coord
            coord = coords_by_index[index]['CA']
            # create a new CA atom
            atom=Atom(atom_id=next_atom_id, atom_name='CA', 
                      x_coord=coord[0], y_coord=coord[1], z_coord=coord[2],
                      element='C')
            
            # create a new monomer
            monomer_name=loop_ind_num_to_monomer_id[index]
            
            # get charge
            if monomer_name in ['ASP', 'GLU']:
                charge = -1
            elif monomer_name in ['LYS', 'ARG']:
                charge = 1
            else:
                charge = 0
            # make new monomer
            new_monomer = Monomer(monomer_number=index, 
                                  monomer_name=monomer_name,
                                  occupancy=1.0,
                                  b_factor=0,
                                  charge=charge)
            # add atom to monomer
            new_monomer.add_atom(atom)

            # add the monomer to the chain
            loop_domain.add_monomer(new_monomer)

    # return the loop domain
    return loop_domain


