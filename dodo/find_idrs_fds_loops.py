'''
Functionality to automate the coarsegraining process using inter-atomic distances. Requires
an AF2 structure (or all atom structure). Was developed using AF2 structures.

Can also use metapredict, but if you need this to be accurate I don't recomend using metapredict
because being off by only a few residues can be an issue for simulations or even visualization
in some cases.

Can identify IDRs, FDs, and FDs with loops (if using non-metapredict approach).
'''

import math
import metapredict as meta
import numpy as np
import os



def get_close_atoms(pdb_file_dict, thresh=8, loop_thresh=7):
    '''
    Gets the number of atoms that are within a threshold value.

    Parameters
    ----------
    pdb_file_dict : dict
        a dictionary of the structure info. Keys are:
        all_atom_coords : list
            a list of the coordinates of all atoms in the structure.
        ca_coords : list
            a list of the coordinates of all CA atoms in the structure.
        all_atoms_types : list
            a list of the atom types of all atoms in the structure.
        sequence : string
            the sequence of the structure in one letter code.
        ca_three_letter_seq : list
            a list of the three letter amino acid codes of all CA atoms in the structure.
        all_atom_three_letter_seq : list
            a list of the three letter amino acid codes of all atoms in the structure.
        coords_by_aa : dict
            For each amino acid, coords for every atom. Needed for the functionality to
            identify loops vs. IDRs

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
    # load chain, get sequence
    sequence=pdb_file_dict['sequence']
    all_coords=pdb_file_dict['coords_by_aa']
    CA_atoms=pdb_file_dict['ca_coords']
    
    # list to hold number of close atoms for each of the amino acids
    close_list=[]
    # list to hold number of close alpha carbons for each amino acid
    potential_loops_close=[]
    for res_ind_1 in range(1, len(sequence)):
        num_close=0
        num_close_loop=0
        for res_ind_2 in range(1, len(sequence)):
            if res_ind_1 != res_ind_2:
                atoms_1 = all_coords[res_ind_1]
                atoms_2 = all_coords[res_ind_2]
                for xyz1 in atoms_1:
                    for xyz2 in atoms_2:
                        dist = math.sqrt(((xyz1[0]-xyz2[0])**2)+((xyz1[1]-xyz2[1])**2)+((xyz1[2]-xyz2[2])**2))
                        if dist < thresh:
                            num_close+=1
                caxyz1=CA_atoms[res_ind_1]
                caxyz2=CA_atoms[res_ind_2]
                ca_dist= math.sqrt(((caxyz1[0]-caxyz2[0])**2)+((caxyz1[1]-caxyz2[1])**2)+((caxyz1[2]-caxyz2[2])**2))
                if ca_dist < loop_thresh:
                    num_close_loop+=1
        close_list.append(num_close)
        potential_loops_close.append(num_close_loop)

    # return the two lists.
    return [close_list, potential_loops_close]


def get_fds_loops_idrs(pdb_file_dict, threshold=480, gap_thresh=25, 
    distance_thresh=8, loop_cutoff=6, min_loop_len=10):
    '''
    function to get the amino acid coordinates that define each
    folded domain, folded domain with loops (and loop coords), 
    and the IDRs from an AF2 pdb. In theory works with any PDB
    that has all atoms, *but* was optimized on AF2 pdbs.

    Parameters
    ----------
    pdb_file_dict : dict
        a dictionary of the structure info. Keys are:
        all_atom_coords : list
            a list of the coordinates of all atoms in the structure.
        ca_coords : list
            a list of the coordinates of all CA atoms in the structure.
        all_atoms_types : list
            a list of the atom types of all atoms in the structure.
        sequence : string
            the sequence of the structure in one letter code.
        ca_three_letter_seq : list
            a list of the three letter amino acid codes of all CA atoms in the structure.
        all_atom_three_letter_seq : list
            a list of the three letter amino acid codes of all atoms in the structure.
        coords_by_aa : dict
            For each amino acid, coords for every atom. Needed for the functionality to
            identify loops vs. IDRs

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

    # get the number of close atoms for each amino acid and the number of close alpha carbons
    close_atoms=get_close_atoms(pdb_file_dict, thresh=distance_thresh)
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
                    fds_bounds.append(start_fd, len(list_of_distances))

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
    len_seq=len(close_atoms[0])+1

    # get idr coords
    if non_idr_coords == []:
        idr_coords=[0, len_seq-1]
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

    # now make dict of regions
    regions_dict={}
    region_num=1
    all_regions = non_idr_coords
    all_regions.extend(idr_coords)
    for region in sorted(all_regions):
        if region in final_fds:
            regions_dict[f'folded_{region_num}']=region
        if str(region) in loops.keys():
            regions_dict[f'fd_with_loop_{region_num}']=[region, loops[f'{region}']]
        if region in idr_coords:
            regions_dict[f'idr_{region_num}']=region
        region_num+=1

    # return regions in order with names and coordinates
    return regions_dict


def get_fds_idrs_from_metapredict(pdb_file_dict):
    '''
    Function to identify IDRs using metapredict. 
    Cannot identify loops consistently and is not
    as accurate as the function that uses atom positions. However
    this approach is MUCH faster.

    Parameters
    -----------
    a dictionary of the structure info. Keys are:
        all_atom_coords : list
            a list of the coordinates of all atoms in the structure.
        ca_coords : list
            a list of the coordinates of all CA atoms in the structure.
        all_atoms_types : list
            a list of the atom types of all atoms in the structure.
        sequence : string
            the sequence of the structure in one letter code.
        ca_three_letter_seq : list
            a list of the three letter amino acid codes of all CA atoms in the structure.
        all_atom_three_letter_seq : list
            a list of the three letter amino acid codes of all atoms in the structure.
        coords_by_aa : dict
            For each amino acid, coords for every atom. Needed for the functionality to
            identify loops vs. IDRs

    Returns
    -------
    dict
        dict of specific regions. 
    '''
    seq=pdb_file_dict['sequence']
    seq_disorder=meta.predict_disorder_domains(seq)
    idrs=seq_disorder.disordered_domain_boundaries
    fds = seq_disorder.folded_domain_boundaries
    start_idr={}
    start_fd={}
    if idrs != []:
        for idr in idrs:
            start_idr[idr[0]]=idr
    if fds != []:
        for fd in fds:
            start_fd[fd[0]]=fd


    fin_dict={}
    dom_num=0

    for ind, a in enumerate(seq):
        if ind in list(start_idr.keys()):
            fin_dict[f'idr_{dom_num}']=start_idr[ind]
            dom_num+=1
        if ind in list(start_fd.keys()):
            fin_dict[f'folded_{dom_num}']=start_fd[ind]
            dom_num+=1

    # adjust final value because... somethign weird with metapredict indexing...        
    fin_dict_adjust_val = list(fin_dict.keys())[-1]
    fin_dict_list = fin_dict[fin_dict_adjust_val]
    fin_dict_list[1]=len(seq)-1
    fin_dict[fin_dict_adjust_val]=fin_dict_list

    return fin_dict


