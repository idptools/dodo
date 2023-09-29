# user facing stuff for making AF2 structures look less silly.

import os

from dodo.dodo_exceptions import dodoException
from dodo.get_af2_structure import get_af2_pdb_lines
from dodo.structure_constructor import build_structure, build_idr_from_sequence, predict_e2e, get_res_dist
from dodo.find_idrs_fds_loops import get_fds_loops_idrs, get_fds_idrs_from_metapredict
from dodo.dodo_tools import plot_structure
from dodo.pdb_tools import PDBParser, write_pdb, write_pdb_coords_seq, array
from dodo import parameters

def pdb_from_name(protein_name, out_path='', mode='predicted', CONECT_lines=True,
    include_FD_atoms=True, use_metapredict=False, graph=False, silent=False, verbose=True,
    attempts_per_region=100, attempts_per_coord=5000, 
    bond_length=parameters.CA_bond_length,
    clash_dist=parameters.CA_clash_dist, 
    min_bond_dist=parameters.min_CA_bond_dist,
    max_bond_dist=parameters.max_CA_bond_dist, 
    debugging=False):
    """
    Function to take in an pdb structure and
    return something where the IDRs are.... more realistic.

    """
    # make sure mode is reasonable before getting rolling.
    if mode not in list(parameters.modes.keys()):
        raise dodoException('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')
    # make a dict with info we need for everything.
    cur_pdb=get_af2_pdb_lines(protein_name, silent=silent)
    PDBParserObj = PDBParser(cur_pdb)
    # use metapredict or atom positions...
    if use_metapredict:
        if verbose==True:
            print('Predicting regions using Metapredict V2.')
        PDBParserObj = get_fds_idrs_from_metapredict(PDBParserObj)
    else:
        if verbose==True:
            print('Predicting folded reigons, loops, and IDRs using AF2 structure information.')        
        PDBParserObj = get_fds_loops_idrs(PDBParserObj)
    # build the new structure
    if verbose==True:
        print('Building new structure.') 
        print(PDBParserObj.regions_dict)   
    PDBParserObj = build_structure(PDBParserObj, mode=mode,
                                    attempts_per_region=attempts_per_region,
                                    attempts_per_coord=attempts_per_coord,
                                    bond_length=bond_length, clash_dist=clash_dist,
                                    min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist,
                                    silent=silent, verbose=verbose, debugging=debugging)
    region_info = PDBParserObj.regions_dict

    # make lists based on what to include...
    xyz_list=[]
    residue_names=[]
    residue_indices=[]
    atom_indices=[]
    beta_vals=[]
    atom_names=[]
    CONECT_COORDS=[]
    
    all_coords = PDBParserObj.all_atom_coords_by_index
    if include_FD_atoms==False:
        for aa in all_coords:
            cur_atoms = all_coords[aa]
            for atom in cur_atoms:
                if atom=='CA':
                    xyz_list.append(all_coords[aa][atom])
                    residue_names.append(PDBParserObj.sequence_3aa_by_index[aa]['CA'])
                    residue_indices.append(aa)
                    atom_indices.append(aa)
                    atom_names.append('CA')
                    beta_vals.append(PDBParserObj.beta_vals_by_index[aa])
                    if aa < len(PDBParserObj.sequence)-1:
                        CONECT_COORDS.append([aa, aa+1])

    else:
        atom_count=0
        all_connect_atoms={}
        final_aa = len(PDBParserObj.sequence)-1
        for aa in all_coords:
            cur_atoms = all_coords[aa]
            for atom in cur_atoms:
                xyz_list.append(all_coords[aa][atom])
                residue_names.append(PDBParserObj.sequence_3aa_by_index[aa][atom])
                residue_indices.append(aa)
                atom_indices.append(atom_count)
                atom_names.append(atom)
                beta_vals.append(PDBParserObj.beta_vals_by_index[aa])
                '''
                need to fix this. Doesn't always go CA - N - C, etc. Numbers chang.e
                '''
                if atom_count < PDBParserObj.number_atoms-1:
                    if atom == 'CA':
                        all_connect_atoms[atom_count]=atom
                    if atom == 'C':
                        all_connect_atoms[atom_count]=atom
                    if atom=='N':
                        all_connect_atoms[atom_count]=atom
                atom_count+=1
        # make conect lines
        atoms_to_connect = list(all_connect_atoms.keys())
        for atom_ind in range(0, len(atoms_to_connect)-1):
            CONECT_COORDS.append([atoms_to_connect[atom_ind], atoms_to_connect[atom_ind+1]])

    if graph==True:
        plot_structure(structure_dict['ca_coords'], region_info, save_path=out_path)
    else:
        if out_path=='':
            raise dodoException('Please specify an output path.')
        else:
            write_pdb(xyz_list, out_path, atom_indices=atom_indices,
            atom_names=atom_names, residue_indices=residue_indices,
            residue_names=residue_names, beta=beta_vals,CONECT_LINES=CONECT_COORDS)



def pdb_from_pdb(path_to_pdb, out_path='', mode='predicted', CONECT_lines=True,
    include_FD_atoms=True, use_metapredict=False, graph=False, silent=False, verbose=True,
    attempts_per_region=100, attempts_per_coord=5000, 
    bond_length=parameters.CA_bond_length,
    clash_dist=parameters.CA_clash_dist, 
    min_bond_dist=parameters.min_CA_bond_dist,
    max_bond_dist=parameters.max_CA_bond_dist, 
    debugging=False, regions_dict=None):
    """
    Function to take in an AF2 pdb structure and
    return something where the IDRs are.... more realistic.

    """

    # make sure mode is reasonable before getting rolling.
    if mode not in list(parameters.modes.keys()):
        raise dodoException('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')    
    
    if os.path.isfile(path_to_pdb)==False:
        raise dodoException('The path you provided to the pdb file does not exist.')
    # make a dict with info we need for everything.
    cur_pdb=open(path_to_pdb, 'r').read().split('\n')
    PDBParserObj = PDBParser(cur_pdb)
    # use metapredict or atom positions...
    if regions_dict==None:
        if use_metapredict:
            if verbose==True:
                print('Predicting regions using Metapredict V2.')
            PDBParserObj = get_fds_idrs_from_metapredict(PDBParserObj)
        else:
            if verbose==True:
                print('Predicting folded reigons, loops, and IDRs using AF2 structure information.')        
            PDBParserObj = get_fds_loops_idrs(PDBParserObj)
        # build the new structure
        if verbose==True:
            print('Building new structure.')   
            print(PDBParserObj.regions_dict) 
    else:
        PDBParserObj.regions_dict=regions_dict
    PDBParserObj = build_structure(PDBParserObj, mode=mode,
                                    attempts_per_region=attempts_per_region,
                                    attempts_per_coord=attempts_per_coord,
                                    bond_length=bond_length, clash_dist=clash_dist,
                                    min_bond_dist=min_bond_dist, max_bond_dist=max_bond_dist,
                                    silent=silent, verbose=verbose, debugging=debugging)
    region_info = PDBParserObj.regions_dict

    # make lists based on what to include...
    xyz_list=[]
    residue_names=[]
    residue_indices=[]
    atom_indices=[]
    beta_vals=[]
    atom_names=[]
    CONECT_COORDS=[]
    
    all_coords = PDBParserObj.all_atom_coords_by_index
    if include_FD_atoms==False:
        for aa in all_coords:
            cur_atoms = all_coords[aa]
            for atom in cur_atoms:
                if atom=='CA':
                    xyz_list.append(all_coords[aa][atom])
                    residue_names.append(PDBParserObj.sequence_3aa_by_index[aa]['CA'])
                    residue_indices.append(aa)
                    atom_indices.append(aa)
                    atom_names.append('CA')
                    beta_vals.append(PDBParserObj.beta_vals_by_index[aa])
                    if aa < len(PDBParserObj.sequence)-1:
                        CONECT_COORDS.append([aa, aa+1])

    else:
        atom_count=0
        all_connect_atoms={}
        final_aa = len(PDBParserObj.sequence)-1
        for aa in all_coords:
            cur_atoms = all_coords[aa]
            for atom in cur_atoms:
                xyz_list.append(all_coords[aa][atom])
                residue_names.append(PDBParserObj.sequence_3aa_by_index[aa][atom])
                residue_indices.append(aa)
                atom_indices.append(atom_count)
                atom_names.append(atom)
                beta_vals.append(PDBParserObj.beta_vals_by_index[aa])
                '''
                need to fix this. Doesn't always go CA - N - C, etc. Numbers chang.e
                '''
                if atom_count < PDBParserObj.number_atoms-1:
                    if atom == 'CA':
                        all_connect_atoms[atom_count]=atom
                    if atom == 'C':
                        all_connect_atoms[atom_count]=atom
                    if atom=='N':
                        all_connect_atoms[atom_count]=atom
                atom_count+=1
        # make conect lines
        atoms_to_connect = list(all_connect_atoms.keys())
        for atom_ind in range(0, len(atoms_to_connect)-1):
            CONECT_COORDS.append([atoms_to_connect[atom_ind], atoms_to_connect[atom_ind+1]])

    if graph==True:
        plot_structure(structure_dict['ca_coords'], region_info, save_path=out_path)
    else:
        if out_path=='':
            raise dodoException('Please specify an output path.')
        else:
            write_pdb(xyz_list, out_path, atom_indices=atom_indices,
            atom_names=atom_names, residue_indices=residue_indices,
            residue_names=residue_names, beta=beta_vals,CONECT_LINES=CONECT_COORDS)



def idr(sequence, out_path='', mode='predicted', attempts_all_coords=250000,
    bond_length=parameters.CA_bond_length, clash_dist=parameters.CA_clash_dist, 
    min_bond_dist=parameters.min_CA_bond_dist, max_bond_dist=parameters.max_CA_bond_dist, 
    start_coord=(0,0,0), CONECT_lines=True, graph=False):
    '''
    Function to make IDR coordinates for an IDR with a specific e2e
    based on mode and sequence
    '''

    # make sure mode is reasonable before getting rolling.
    if mode not in list(parameters.modes.keys()):
        raise dodoException('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')        

    # set mode
    run_prediction=False
    if mode =='super_compact':
        length_multiplier = 0.3
    elif mode == 'compact':
        length_multiplier=0.55
    elif mode == 'normal':
        length_multiplier=0.8
    elif mode == 'expanded':
        length_multiplier = 1.05
    elif mode == 'super_expanded':
        length_multiplier = 1.3
    elif mode == 'max_expansion':
        length_multiplier=1.65
    elif mode == 'predicted':
        run_prediction=True
    else:
        raise Exception('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')
    # get len IDR
    idr_len = len(sequence)
    if CONECT_lines:
        CONECT_coords=[]
        for val in range(1, len(sequence)):
            CONECT_coords.append([val, val+1])

    else:
        CONECT_coords=None
    if run_prediction:
        # get predicted e2e
        e2e = predict_e2e(sequence)
    else:
        e2e=length_multiplier*idr_len
    idr_coords = build_idr_from_sequence(e2e, sequence,
    bond_length=bond_length, clash_dist=clash_dist, 
    attempts_all_coords=attempts_all_coords, min_bond_dist=min_bond_dist,
    max_bond_dist=max_bond_dist, start_coord=start_coord, all_atoms=False)
    if graph==False:
        if out_path=='':
            raise dodoException('Please specify an output path.')
        write_pdb_coords_seq(idr_coords, sequence, out_path, CONECT_LINES=CONECT_coords)
    else:
        plot_structure(idr_coords, {'idr': [0, len(sequence)]})





