'''
Final backend module for the DODO project.
This module contains everything necessary for actually reconstructing
the structure of a protein. This module will ultimately 
be imported and used for all of the frontend functionality.

There will be two general functionalities.
1) rebuild structure
    This allows repositioning of FDs within the structure relative to one another.
2) rebuild assembly
    This will rebuild IDRs in an assembly but will not move anything relative to one another.
'''


import numpy as np
import os

from dodo.backend.generate_alpha_carbon_points import generate_ca_points
import idr_construction_utils as idr_utils
import position_fds
import utils
from dodo_readers import Reader
from dodo_structures import Atom, Monomer, Domain, Chain, Complex
from dodo.backend.domain_identification import get_fds_loops_idrs, assign_domains_from_dict

def rebuild_structure(path_to_pdb, num_conformations=1,
                      idr_expansion = 'standard',
                      linear_positioning=False,
                      manual_domain_assignment=None):
    '''
    This function will rebuild the structure of a protein.
    It will take in a path to a PDB file and return a Complex object
    that represents the protein structure. The complex object is then
    used to rebuild everything. Importantly, we can choose to
    only reposition the FDs once and then build multiple IDRs, which
    are inherently random to some extent. This lets us build multiple 
    configurations of the IDRs for a single structure. 

    Parameters
    ----------
    path_to_pdb : str
        path to the PDB file

    num_conformations : int
        number of conformations to generate for the IDRs

    idr_expansion : str
        How expanded we want the IDR to be. Options are 'compact',
        'standard', 'expanded'. Default is 'standard'.
        Standard will (eventually) use the predicted end-to-end
        distance for the IDR, but we can do that later. 

    manual_domain_assignment : dict
        Dictionary mapping each domain to the indices of the
        residues that belong to that domain. This will
        let people override the automatic domain assignment.
    
    linear_positioniong : bool
        If True, the FDs will be positioned linearly. If False,
        the FDs will be not be constrained to be linearly positioned. 
        Default is False.

    Returns
    -------
    Complex
        Complex object representing the protein structure
    '''
    # read in the PDB file and make a DodoComplex
    DodoComplex = Reader(path_to_pdb).return_complex()
    
    if manual_domain_assignment==None:
        # now identify the domains.
        DodoComplex = get_fds_loops_idrs(DodoComplex)
    else:
        DodoComplex = assign_domains_from_dict(DodoComplex, manual_domain_assignment)

    # DELETE NEW_COORDS LATER!
    new_coords=[]    

    # now we want to iterate over the chains of the domain. 
    for chain in DodoComplex.get_chains():
        # get current domains
        all_domains = chain.get_domains()
        # get domains that are fds, regardless of if they have loops. 
        fd_indices = [d.domain_id for d in all_domains if d.domain_type=='FD' or d.domain_type=='FD_with_loop']
        # get the loop indices
        loop_indices = [d.domain_id for d in all_domains if d.domain_type=='loop']
        # get the IDR indices
        idr_indices = [d.domain_id for d in all_domains if d.domain_type=='IDR']
        # if we only have 1 FD, we don't need to do anything.
        if len(fd_indices)<=1:
            continue
        # if we have more than one FD, we need to reposition them.
        else:
            # get FD1, translate it such that the final coord in FD1 is at the origin. 
            FD1 = chain.get_domain(fd_indices[0])
            FD1.translate_domain_to_origin()


            # DELETE LATER
            # DELETE LATER
            # DELETE LATER
            # DELETE LATER
            # DELETE NEW_COORDS LATER!
            new_coords.extend(list(FD1.get_coordinates_array()))
            
            # iterate over FD indices. We always move the next FD relative to the previous FD,
            # so we need to make sure that we don't go out of bounds.
            for i in range(0, len(fd_indices)-1):
                # get FD1
                FD1 = chain.get_domain(fd_indices[i])
                
                # if this is not the final domain in the chain, make sure
                # that the next domain is an IDR.
                if i!=len(fd_indices)-1:
                    if all_domains[FD1.domain_id+1].domain_type!='IDR':
                        continue
                
                '''
                NOTE!
                We need to modify this later to be predicted distance for standard,
                but for now we will just use a distance of 1.5 times the length of the
                next domain. Also later need to add in compact and expanded.
                '''

                if idr_expansion!='standard':
                    raise NotImplementedError('This functionality is not yet implemented.')

                distance=1.4*len(all_domains[FD1.domain_id+1].sequence)

                # reposition the FDs
                if linear_positioning==False:
                    position_fds.position_folded_domain(fd_indices[i], fd_indices[i+1], chain, distance)
                else:
                    position_fds.position_folded_domain_linear(fd_indices[i], fd_indices[i+1], chain, distance)

                # calculate distance between the two FDs
                # DELETE LATER
                # DELETE LATER
                # DELETE LATER
                # DELETE LATER
                FD1 = chain.get_domain(fd_indices[i])
                FD2 = chain.get_domain(fd_indices[i+1])
                new_coords.extend(list(FD2.get_coordinates_array()))

        # now we need to build loops
        if len(loop_indices) > 0:
            for loop_index in loop_indices:
                loop = chain.get_domain(loop_index)
                # keep going from here!
                loop.build_loop()

        # write the new coordinates to a PDB file
        utils.write_dumby_pdb(new_coords, f'{os.getcwd()}/dodo/data/testing_translation.pdb')
            

rebuild_structure(f'{os.getcwd()}/dodo/data/p300.cif', linear_positioning=False)
