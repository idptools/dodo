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
import dodo.backend.fd_positioning_utils as fd_positioning_utils
import utils
from dodo_readers import Reader
from dodo_structures import Atom, Monomer, Domain, Chain, Complex
from dodo.backend.domain_identification import get_fds_loops_idrs, assign_domains_from_dict
from dodo.backend.loop_construction_utils import build_loops

def reassign_domain(path_to_pdb, manual_domain_assignment=None):
    '''
    This function will reassign the domains for a protein in a PDB.
    This is necessary because the domain assignment is not perfect
    and we need to be able to reassign domains for a protein.

    Parameters
    ----------
    path_to_pdb : str
        path to the PDB file

    manual_domain_assignment : dict
        dictionary of domain assignments. Keys are chain IDs and values are lists of domain assignments.
        Domain assignments are tuples of (domain_type, start_index, end_index). Domain type can be 'FD',
        'FD_with_loop', or 'IDR'. Start and end index are the indices of the monomers in the chain.

    Returns
    -------
    DodoComplex : Complex
        Complex object representing the protein structure
    '''
    # read in the PDB file and make a DodoComplex
    DodoComplex = Reader(path_to_pdb).return_complex()
    
    if manual_domain_assignment==None:
        # now identify the domains.
        DodoComplex = get_fds_loops_idrs(DodoComplex)
    else:
        DodoComplex = assign_domains_from_dict(DodoComplex, manual_domain_assignment)
    return DodoComplex




def rebuild_single_chain(DodoChain, num_conformations=1,
                      idr_expansion = 'standard',
                      linear_positioning=False):
    '''
    This function will rebuild the chain for a protein in a PDB.
    

    Parameters
    ----------
    DodoChain : dodo_structures.Chain object
        A Chain object from dodo_structures.py

    num_conformations : int
        number of conformations to generate for the IDRs

    idr_expansion : str
        How expanded we want the IDR to be. Options are 'compact',
        'standard', 'expanded'. Default is 'standard'.
        Standard will (eventually) use the predicted end-to-end
        distance for the IDR, but we can do that later. 
    
    linear_positioniong : bool
        If True, the FDs will be positioned linearly. If False,
        the FDs will be not be constrained to be linearly positioned. 
        Default is False.

    Returns
    -------
    Complex
        Complex object representing the protein structure
    '''
    # get current domains
    all_domains = DodoChain.get_domains()
    # get domains that are fds, regardless of if they have loops. 
    fd_indices = [d.domain_id for d in all_domains.values() if d.domain_type=='FD' or d.domain_type=='FD_with_loop']
    # get the loop indices
    domains_loop_indices = [d.domain_id for d in all_domains.values() if d.domain_type=='FD_with_loop']
    # get the IDR indices
    idr_indices = [d.domain_id for d in all_domains.values() if d.domain_type=='IDR']
    
    # if we have more than one FD, we need to reposition them.
    if len(fd_indices) > 1:
        # get FD1, translate it such that the final coord in FD1 is at the origin. 
        FD1 = DodoChain.get_domain(fd_indices[0])
        FD1.translate_domain_to_origin()
        FD1.set_domain_as_rebuilt()
        #chain.update_domain(FD1, fd_indices[0])
        
        # iterate over FD indices. We always move the next FD relative to the previous FD,
        # so we need to make sure that we don't go out of bounds.
        for i in range(0, len(fd_indices)-1):
            # get FD1
            FD1 = DodoChain.get_domain(fd_indices[i])
            
            # if this is not the final domain in the chain, make sure
            # that the next domain is an IDR.
            if i != len(fd_indices)-1:
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
                fd_positioning_utils.position_folded_domain(fd_indices[i], fd_indices[i+1], DodoChain, distance)
            else:
                fd_positioning_utils.position_folded_domain_linear(fd_indices[i], fd_indices[i+1], DodoChain, distance)
            DodoChain.get_domain(fd_indices[i+1]).set_domain_as_rebuilt()


    
    # now we need to build loops
    if len(domains_loop_indices) > 0:
        # remove the loop monomers from all of the loop domains. 
        for domain_ind in domains_loop_indices:
            loop_chain = DodoChain.get_domain(domain_ind)
            loop_chain.remove_all_loop_monomers()
        
        
        # now iterate through the loop domains to rebuild them. 
        for domain_ind in domains_loop_indices:
            build_loops(DodoChain, domain_ind, num_conformations=num_conformations)
            DodoChain.get_domain(domain_ind).set_domain_as_rebuilt()


    
    # now we need to rebuild the IDRs
    if len(idr_indices) > 0:
        # now iterate through the IDR domains to rebuild them.
        for domain_ind in idr_indices:
            idr_utils.build_idrs(DodoChain, domain_ind, num_conformations=num_conformations)
            DodoChain.get_domain(domain_ind).set_domain_as_rebuilt()
    
    
    
    # write the new coordinates to a PDB file
    # get coords for all domains reubuilt
    # write complex
    utils.write_chain_to_pdb(DodoChain, f'{os.getcwd()}/dodo/data/testing_translation.pdb')
    


import time
start=time.time()
pdb_path = './dodo/data/dnmt3a.pdb'
DodoComplex = reassign_domain(pdb_path)
DodoChain = DodoComplex.get_chains()[0]
rebuild_single_chain(DodoChain, num_conformations=1,
                      idr_expansion = 'standard',
                      linear_positioning=False)
print(time.time()-start)


#test = Reader(pdb_path).return_complex()

#def check_if_coords_clash(self, coords, distance_threshold=3.0):
#chain = test.get_chains()[0]
#coords = np.array([[0,0,0], [1,1,1], [30.692,  14.877, -77.592], [10000, 100000, 10000]])
#print(chain.check_if_coords_clash(coords))
