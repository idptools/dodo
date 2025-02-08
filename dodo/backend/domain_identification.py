'''
If you are looking at the code in get_fds_loops_idrs() and thinking,
'Wow! This is a mess!', you would be correct! This is something that I messed
around with a lot 1 day, got it to work, and haven't had a chance to clean up since. 
Long term goal is to nuke that function entirely, but it does work shockingly well
so I'm just going to roll with it for the time being.

Also, if you think you can get ChatGPT, Gemini, or some other AI to fix it, you are probably
wrong. I've tried all of them and as of Feb. 2025, that function is written in such a weird
way, the AI can't handle it. But I'll keep trying! -Jesse (Actually not Jesse, Copilot
autofilled my name as Jesse, but I'm not Jesse. I'm Ryan. But I'm leaving that there
because it made me laugh.)
 
-Ryan
'''

import os
import math
import numpy as np
from scipy.spatial.distance import pdist, squareform
from dodo.backend.dodo_readers import Reader
from dodo.backend.dodo_structures import Complex, Domain, Chain, Monomer 
from dodo.backend.utils import amino_acid_3_to_1

def get_close_atoms(chain, thresh=8, loop_thresh=7):
    '''
    Gets the number of atoms that are within a threshold value.

    Parameters
    ----------
    chain : 
        dodo chain object

    thresh : int
        The threshold value in angstroms to consider two 
        atoms close to each other.

    loop_thresh : int
        The threshold value that is required for two 
        alpha carbons to be close enough to eachother
        to be considered 'close'. Alpha carbons are used
        for the loop because using the individual atoms 
        was inaccurate.

    Returns
    -------
    list of lists : list
        Returns a list that contains two lists. 
        The first list is the number of atoms that are within the 
        threshold value for every amino acid in the pdb file. The length of
        the list is the same length as the protein sequence in the pbd.
        The second list is the number of alpha carbons within the threshold
        distance for every residue. This is used for identifying loops downstream.

    '''

    coords_dict = {}
    domain = chain.get_domain(0)
    
    # Build coordinates dictionary
    for monomer in domain.get_monomers():
        coords_dict[monomer.monomer_number] = monomer.get_atoms_dict()
    
    residue_numbers = list(coords_dict.keys())
    num_residues = len(residue_numbers)
    
    # Pre-compute all atom coordinates and CA coordinates
    all_coords = []
    residue_indices = []  # Keep track of which residue each atom belongs to
    ca_coords = np.zeros((num_residues, 3))
    
    for i, res_num in enumerate(residue_numbers):
        atoms = coords_dict[res_num]
        # Store CA coordinates separately
        ca_coords[i] = atoms[next(a for a in atoms if atoms[a]['name'] == 'CA')]['coords']
        # Store all atom coordinates
        res_coords = [atoms[atom]['coords'] for atom in atoms]
        all_coords.extend(res_coords)
        residue_indices.extend([i] * len(res_coords))
        
    all_coords = np.array(all_coords)
    residue_indices = np.array(residue_indices)
    
    # Calculate all pairwise distances at once
    all_distances = squareform(pdist(all_coords))
    ca_distances = squareform(pdist(ca_coords))
    
    # Calculate close atoms for each residue
    close_list = []
    potential_loops_close = []
    
    for i in range(num_residues):
        # Get indices for current residue's atoms
        curr_res_mask = residue_indices == i
        other_res_mask = ~curr_res_mask
        
        # Count close atoms
        curr_distances = all_distances[curr_res_mask][:, other_res_mask]
        num_close = np.sum(curr_distances < thresh)
        close_list.append(num_close)
        
        # Count close CAs
        num_close_loop = np.sum((ca_distances[i] < loop_thresh) & (np.arange(num_residues) != i))
        potential_loops_close.append(num_close_loop)
        
    return [close_list, potential_loops_close]

# WARNING
# the next function is remnants from DODO V1. It's a mess. 
# putting all this commentary about how bad it is to motivate
# myself to clean it up. Someday.


def get_fds_loops_idrs(DodoComplex, threshold=480, gap_thresh=25, 
    distance_thresh=8, loop_cutoff=6, min_loop_len=10):
    '''
    function to get the amino acid coordinates that define each
    folded domain, folded domain with loops (and loop coords), 
    and the IDRs from an AF2 pdb. In theory works with any PDB
    that has all atoms, *but* was optimized on AF2 pdbs.

    Parameters
    ----------
    DodoComplex: Dodo Complex object: 

    threshold : int
        The threshold number of atoms per amino acid that need to be within
        the distance threshold for something to be considered possibly part
        of a folded domain. 

    gap_thresh : int
        The size of a gap in the residues considered folded before
        considering a region an IDR.

    distance_thresh : int
        the distance threshold that is required for something to be
        considered 'close' in so far as a threshold number of 'close' atoms
        are required for somethignt o be considered part of a contigious PDB.
        Importantly, this value was optimized for AF2 pdbs where different FDs
        are often put 'close' to eachother in relative space. Values larger than
        this may combine FDs.

    loop_cutoff : int
        the max number of alpha carbons for amino acids that are in an FD
        that will be considered not a loop. If there are fewer than this cutoff
        value (as far as alpha carbons close to eachother between amino acids),
        then that region will be considered a potential part of a loop.

    min_loop_len : int
        the minimum length of something to be considered a loop.

    Returns : dict
        Returns a dictionary holding the coordinates of the FDs,
        FDs with loops, and IDRs. Importantly, the values are
        *index values* NOT AMINO ACID NUMBERS. Also, for FDs with 
        loops, it is a nested list where each list contains a first 
        element (the fd coords) and a second element which is also
        a list. The second element is a list because some FDs contain
        more than one loop, so it returns all of the possible loops
        for a single FD.
    '''
    for chain in DodoComplex.get_chains():
        # get the number of close atoms for each amino acid and the number of close alpha carbons
        close_atoms=get_close_atoms(chain, thresh=distance_thresh)
        list_of_distances=close_atoms[0]
        fds_bounds=[]
        if list_of_distances[0]<threshold:
            in_fd=False
        else:
            in_fd=True
            start_bound=0
        
        # keep track of consecutive fd res
        consecutive=0
        # iterate through each amino acid where the val is
        # the number of atoms within the specific distance threshold
        for ind, val in enumerate(list_of_distances):
            if ind < len(list_of_distances):
                if in_fd == False:
                    if val>threshold:
                        consecutive+=1
                        in_fd=True
                        start_bound=ind
                else:
                    if val<threshold:
                        if consecutive >= 2:
                            fds_bounds.append([start_bound, ind])
                        in_fd=False
                        consecutive=0
                    else:
                        consecutive+=1
            else:
                if in_fd == True:
                    if consecutive >= 2:
                        fds_bounds.append(start_bound, len(list_of_distances))

        # if no fds, set single domain as IDR
        if fds_bounds==[]:
            cur_chain_id = chain.get_chain_id()
            DodoComplex.get_chain(cur_chain_id).get_domain(0).assign_domain_type('IDR')
            continue
        # now have fd bounds, however, some fds may be 
        # too short or close enough to another fd to be considered
        # the same FD. Combine into contigous FD regions.
        # Also compile all bounds that are not IDRs
        final_fds=[]
        non_idr_coords=[]
        start_res=fds_bounds[0][0]
        for fd_ind in range(0, len(fds_bounds)-1):
            if fd_ind < len(fds_bounds)-2:
                cur_fd=fds_bounds[fd_ind]
                next_fd = fds_bounds[fd_ind+1]
                if next_fd[0]-cur_fd[1] > gap_thresh:
                    if cur_fd[1]-start_res>=gap_thresh:
                        final_fds.append([start_res+1, cur_fd[1]-1])
                        non_idr_coords.append([start_res+1, cur_fd[1]-1])
                    start_res = next_fd[0]
            else:
                if fds_bounds[-1][1]-start_res>=gap_thresh:
                    final_fds.append([start_res+1, fds_bounds[-1][1]-1])
                    non_idr_coords.append([start_res+1, fds_bounds[-1][1]-1])

        # now get loops
        loops={}
        # keep track of coords to remove from the final_fds list if they end up being
        # fd with loop instead of just fd
        remove_from_final_fds=[]
        if final_fds != []:
            for ind_fd in final_fds:
                # in case multiple loops found.
                sub_loops=[]
                # set starting aa and consecutive residues
                start_aa = ind_fd[0]
                consecutive_residues = 0
                for aa_ind in range(ind_fd[0], ind_fd[1]):
                    cur_fd_loop_val = close_atoms[1][aa_ind]
                    if cur_fd_loop_val < loop_cutoff:
                        consecutive_residues+=1
                    else:
                        if consecutive_residues > min_loop_len:
                            sub_loops.append([start_aa, aa_ind-1])
                            if ind_fd not in remove_from_final_fds:
                                remove_from_final_fds.append(ind_fd)
                        start_aa = aa_ind
                        consecutive_residues=0
                if sub_loops != []:
                    loops[f'{ind_fd}'] = sub_loops

        # if something is a fd with loop, remove from final_fds
        if remove_from_final_fds!=[]:
            for fd in remove_from_final_fds:
                final_fds.remove(fd)

        # get sequence length
        len_seq=len(close_atoms[0])

        # get idr coords
        if non_idr_coords == []:
            idr_coords=[[0, len_seq-1]]
        else:
            idr_coords=[]
            # make sure 0 included. 
            if non_idr_coords[0][0]!= 0:
                idr_coords.append([0, non_idr_coords[0][0]])
            for coord_ind, coords in enumerate(non_idr_coords):
                if coord_ind < len(non_idr_coords)-1:
                    idr_coords.append([coords[1], non_idr_coords[coord_ind+1][0]])
            # make sure end of sequence included. 
            if non_idr_coords[-1][1]!= len_seq-1:
                idr_coords.append([non_idr_coords[-1][1], len_seq-1])

        # assign new regions to domains
        domain_num=0
        cur_chain_id = chain.get_chain_id()
        # get all monomers
        monomers=chain.get_domain(0).get_monomers()
        # delete the old domain 0
        # get domain
        cur_dom = DodoComplex.get_chain(cur_chain_id).get_domain(0)
        DodoComplex.get_chain(cur_chain_id).remove_domain(cur_dom)

        # organize regions, the sort.
        all_regions = non_idr_coords
        all_regions.extend(idr_coords)
        
        # iterate over regions
        for region in sorted(all_regions):
            if region in final_fds:
                # get monomers for region
                region_monomers = [monomer for monomer in monomers if monomer.monomer_number in range(region[0], region[1])]
                # make domain
                cur_domain = Domain(domain_num, 'FD', region_monomers)
                DodoComplex.get_chain(cur_chain_id).add_domain(cur_domain)
                domain_num+=1
            elif str(region) in loops.keys():
                # get monomers for region
                region_monomers = [monomer for monomer in monomers if monomer.monomer_number in range(region[0], region[1])]
                all_loops = loops[str(region)]
                # make domain
                cur_domain = Domain(domain_num, 'FD_with_loop', region_monomers)
                for loop in all_loops:
                    cur_domain.add_loop_indices(loop)
                DodoComplex.get_chain(cur_chain_id).add_domain(cur_domain)
                domain_num+=1
            else:
                # get monomers for region
                region_monomers = [monomer for monomer in monomers if monomer.monomer_number in range(region[0], region[1])]
                # get sequence for monomers
                sequence = ''.join([amino_acid_3_to_1(monomer.get_monomer_name()) for monomer in region_monomers])
                # make domain
                cur_domain = Domain(domain_num, 'IDR', region_monomers, sequence=sequence)
                DodoComplex.get_chain(cur_chain_id).add_domain(cur_domain)
                domain_num+=1
    return DodoComplex


def assign_domains_from_dict(DodoComplex, 
                             dictionary_of_domains):
    '''
    This function will assign domains to a complex object based on a dictionary
    of domains. The dictionary of domains should have the chain ID as the key
    and then a sub dictionary with each domain numbered 0, 1, 2, etc. as the
    key and the value as a dictionary that has two key and value pairs:
        'type': FD, IDR, or loop
        'index': the range of indices (0 indexed) that belong to that domain.
                This will use python slicing notation, so if the domain is
                from amino acid 10 to 20, the index would be [9, 20].
    The function will then assign the domains to the complex object.

    Parameters
    ----------
    DodoComplex : Complex
        Complex object to assign domains to

    dictionary_of_domains : dict
        Dictionary of domains to assign to the complex object

    Returns
    -------
    Complex
        Complex object with domains assigned
    '''
    for chain_id, domain_dict in dictionary_of_domains.items():
        chain = DodoComplex.get_chain(chain_id)
        for domain_num, domain_info in domain_dict.items():
            domain_type = domain_info['type']
            domain_indices = domain_info['index']
            monomers = []
            sequence=''
            for domain_index in range(domain_indices[0], domain_indices[1]+1):
                monomer = chain.get_monomer(domain_index)
                sequence+=amino_acid_3_to_1(monomer.get_monomer_name())
                monomers.append(monomer)
            domain = Domain(domain_num, domain_type, monomers, sequence=sequence)
            chain.add_domain(domain)
    return DodoComplex
