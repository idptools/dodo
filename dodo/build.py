# user facing stuff for making AF2 structures look less silly.

import os

from dodo.dodo_exceptions import dodoException
from dodo.get_af2_structure import get_af2_pdb
from dodo.structure_constructor import build_structure
from dodo.find_idrs_fds_loops import get_fds_loops_idrs, get_fds_idrs_from_metapredict
from dodo.dodo_tools import af2_lines_to_structure_dict, plot_structure, write_pdb


def pdb_from_name(protein_name, out_path='', mode='predicted', 
    use_metapredict=False, graph=False, beta_by_region=True, silent=False, verbose=False,
    attempts_per_region=50, attempts_per_coord=5000, bond_length=3.8,
    clash_dist=3.4, min_bond_dist=2.8, max_bond_dist=4.4, debugging=False):
    """
    Function to take in an pdb structure and
    return something where the IDRs are.... more realistic.

    """
    # make a dict with info we need for everything.
    cur_pdb=get_af2_pdb(protein_name, save=False, silent=silent)
    structure_dict = af2_lines_to_structure_dict(cur_pdb)
    # use metapredict or atom positions...
    if use_metapredict:
        if verbose==True:
            print('Predicting regions using Metapredict V2.')
        region_info = get_fds_idrs_from_metapredict(structure_dict)
    else:
        if verbose==True:
            print('Predicting folded reigons, loops, and IDRs using AF2 structure information.')        
        region_info = get_fds_loops_idrs(structure_dict)
    # build the new structure
    if verbose==True:
        print('Building new structure.')    
    new_structure = build_structure(structure_dict, region_info, mode=mode,
                                    attempts_per_region=attempts_per_region,
                                    attempts_per_coord=attempts_per_coord,
                                    bond_length=bond_length, clash_dist=clash_dist,
                                    min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist,
                                    silent=silent, verbose=verbose, debugging=debugging)
    if beta_by_region:
        betas=[]
        for region in region_info:
            coords = region_info[region]
            if region==list(region_info.keys())[-1]:
                add_val=1
            else:
                add_val=0
            if 'idr' in region:
                for val in range(coords[0], coords[1]+add_val):
                    betas.append(100)
            elif 'folded' in region:
                for val in range(coords[0], coords[1]+add_val):
                    betas.append(66)
            elif 'loop' in region:
                all_vals=[aa for aa in range(coords[0][0], coords[0][1]+add_val)]
                fd_region = [aa for aa in range(coords[0][0], coords[0][1]+add_val)]
                loop_region = []
                for loop in coords[1]:
                    for val in range(loop[0], loop[1]):
                        loop_region.append(val)
                        fd_region.remove(val)
                for val in all_vals:
                    if val in fd_region:
                        betas.append(66)
                    elif val in loop_region:
                        betas.append(22)
                    else:
                        raise dodoException('Something went wrong with the beta regions')
            else:
                raise dodoException('something went wrong with beta regions..')
    else:
        betas=None

    if graph==True:
        plot_structure(new_structure['ca_coords'], region_info, save_path=out_path)
    else:
        if out_path=='':
            return new_structure
        else:
            xyz_list = new_structure['ca_coords']
            sequence = structure_dict['ca_three_letter_seq']
            write_pdb(xyz_list, out_path, residue_names=sequence, beta=betas)


def pdb_from_pdb(path_to_pdb, mode='predicted', out_path='', 
    use_metapredict=False, graph=False, beta_by_region=True, silent=False, verbose=False,
    attempts_per_region=50, attempts_per_coord=5000, bond_length=3.8,
    clash_dist=3.4, min_bond_dist=2.8, max_bond_dist=4.4, debugging=False):
    """
    Function to take in an AF2 pdb structure and
    return something where the IDRs are.... more realistic.

    """
    if os.path.isfile(path_to_pdb)==False:
        raise dodoException('The path you provided to the pdb file does not exist.')
    # make a dict with info we need for everything.
    cur_pdb=open(path_to_pdb, 'r').read().split('\n')
    structure_dict = af2_lines_to_structure_dict(cur_pdb)
    # if coarse-grained already, just use the ca atoms. 
    if len(structure_dict['ca_coords'])==len(structure_dict['all_atom_coords']):
        if use_metapredict==False:
            if silent != True:
                print('Warning! You selected to not use metapredict for identifying IDRs. To use the other method, all atom structural info is needed! Falling back to metapredict.')
            use_metapredict=True

    # use metapredict or atom positions...
    if use_metapredict:
        if verbose==True:
            print('Predicting regions using Metapredict V2.')
        region_info = get_fds_idrs_from_metapredict(structure_dict)
    else:
        if verbose==True:
            print('Predicting folded reigons, loops, and IDRs using AF2 structure information.')        
        region_info = get_fds_loops_idrs(structure_dict)
    # build the new structure
    if verbose==True:
        print('Building new structure.')    
    new_structure = build_structure(structure_dict, region_info, mode=mode,
                                    attempts_per_region=attempts_per_region,
                                    attempts_per_coord=attempts_per_coord,
                                    bond_length=bond_length, clash_dist=clash_dist,
                                    min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist,
                                    silent=silent, verbose=verbose, debugging=debugging)
    if beta_by_region:
        betas=[]
        for region in region_info:
            coords = region_info[region]
            if region==list(region_info.keys())[-1]:
                add_val=1
            else:
                add_val=0
            if 'idr' in region:
                for val in range(coords[0], coords[1]+add_val):
                    betas.append(100)
            elif 'folded' in region:
                for val in range(coords[0], coords[1]+add_val):
                    betas.append(66)
            elif 'loop' in region:
                all_vals=[aa for aa in range(coords[0][0], coords[0][1]+add_val)]
                fd_region = [aa for aa in range(coords[0][0], coords[0][1]+add_val)]
                loop_region = []
                for loop in coords[1]:
                    for val in range(loop[0], loop[1]):
                        loop_region.append(val)
                        fd_region.remove(val)
                for val in all_vals:
                    if val in fd_region:
                        betas.append(66)
                    elif val in loop_region:
                        betas.append(22)
                    else:
                        raise dodoException('Something went wrong with the beta regions')
            else:
                raise dodoException('something went wrong with beta regions..')
    else:
        betas=None
    if graph==True:
        plot_structure(new_structure['ca_coords'], region_info, save_path=out_path)
    else:
        if out_path=='':
            return new_structure
        else:
            xyz_list = new_structure['ca_coords']
            sequence = structure_dict['ca_three_letter_seq']
            write_pdb(xyz_list, out_path, residue_names=sequence, beta=betas)



