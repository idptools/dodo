# user facing stuff for making AF2 structures look less silly.

import os

from dodo.dodo_exceptions import dodoException
from dodo.get_af2_structure import get_af2_pdb_lines
from dodo.idr_constructor import build_structure, build_idr_from_sequence
from dodo.find_idrs_fds_loops import get_fds_loops_idrs, get_fds_idrs_from_metapredict
from dodo.dodo_tools import plot_structure
from dodo.pdb_tools import PDBParser, write_pdb, array, save_pdb_from_PDBParserObj
from dodo import parameters



def pdb_from_name(protein_name, out_path='', mode='predicted', 
    linear_placement=False, CONECT_lines=True, include_FD_atoms=True, 
    use_metapredict=False, graph=False, verbose=True, attempts_per_region=40, 
    attempts_per_coord=2000, num_models=1, beta_for_FD_IDR=False, just_fds=False):
    """
    Function to take in the name of a protein and then return an AF2 PDB
    with modified disordered regions. 

    Parameters
    ----------
    protein_name : str
        Name of protein to get PDB for. Specifying organism name helps.
        Ex. 'Arabidopsis ARF19'
    out_path : str
        Path to save the PDB to.
    mode : str
        Mode to use for building the structure. Options are 'super_compact',
        'compact', 'normal', 'expanded', 'super_expanded', 'max_expansion', 
        and 'predicted'. 'predicted' uses ALBATROSS to predict IDR end-to-end distance.
    linear_placement : bool
        whether to place the folded domains across a linear axis. 
        Default : False
    CONECT_lines : bool
        Whether to include CONECT lines in the pdb. Default is True.
    include_FD_atoms : bool
        Whether to include all atoms from the AF2 structure. Default is True.
    use_metapredict : bool
        Whether to use Metapredict to predict IDRs. Default is False.
        This is less accurate than using the actual structure but is MUCH faster.
    graph : bool
        Whether to graph the structure. Default is False.
        If you graph the output, it will not save.
    verbose : bool
        Whether to print out info about what is happening. Default is True.
    attempts_per_region : int
        Number of attempts to make per region. Default is 40.
    attempts_per_coord : int
        Number of attempts to make per coordinate. Default is 2000.
    num_models : int
        Number of times to rebuild the IDRs while holding the FDs contant.
        Default : 1
    beta_for_FD_IDR : bool
        Whether to overwrite beta vals where IDR=100 and FD=0 for 
        visualization purposes.
        Default is False
    just_fds : bool
        Whether to just dump the folded domains as individual PDBs.
        Default : False. 

    Returns
    -------
    if graph==False:
        None. Just saves the PDB.
    if graph==True:
        None. Just graphs the structure.

    """
    # don't let people waste their time trying to graph multiple models. 
    if graph==True and num_models!=1:
        num_models=1
        print('Cannot graph more than one model. Falling back to one model.')

    # make sure mode is reasonable before getting rolling.
    if mode not in list(parameters.modes.keys()):
        raise dodoException('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')
    
    # make a dict with info we need for everything.
    PDBParserObj=get_af2_pdb_lines(protein_name, verbose=verbose)
    
    # use metapredict or atom positions, update the regions info in PDBParserObj.
    if use_metapredict:
        if verbose==True:
            print('Predicting regions using Metapredict V2.')
        PDBParserObj = get_fds_idrs_from_metapredict(PDBParserObj)
    else:
        if verbose==True:
            print('Predicting folded reigons, loops, and IDRs using AF2 structure information.')        
        PDBParserObj = get_fds_loops_idrs(PDBParserObj)
    
    # Print some info
    if verbose==True:
        print('Building new structure.') 
        temp=''
        for reg in PDBParserObj.regions_dict:
            cur_reg = PDBParserObj.regions_dict[reg]
            if 'idr' in reg:
                reg_name = 'IDR'
            elif 'folded' in reg:
                reg_name = 'FD'
            else:
                reg_name = 'FD w/ loop(s)'
            if 'loop' in reg:
                fd = cur_reg[0]
                loops = cur_reg[1]
                loop_locs=''
                for loop in loops:
                    loop_locs+=f'{loop[0]+1} to {loop[1]+1}, '
                temp+=f"{reg_name}: FD {fd[0]+1} to {fd[1]+1} with loop(s) at {loop_locs}"
            else:
                temp+=f"{reg_name}: {cur_reg[0]+1} to {cur_reg[1]+1}, "
        print(temp[:len(temp)-2])

    if len(PDBParserObj.regions_dict)==1 and 'idr_1' in PDBParserObj.regions_dict:
        return pdb_from_sequence(PDBParserObj.sequence, out_path=out_path, mode=mode, 
            attempts_per_res=1000, attempts_per_idr=50, end_coord=(0,0,0), 
                CONECT_lines=CONECT_lines, graph=graph, num_models=num_models)


    # build new structure  
    if just_fds==False:
        PDBParserObj = build_structure(PDBParserObj, mode=mode, 
                                        linear_placement=linear_placement,
                                        attempts_per_region=attempts_per_region,
                                        attempts_per_coord=attempts_per_coord,
                                        verbose=verbose, num_models=num_models,
                                        beta_for_FD_IDR=beta_for_FD_IDR)
    else:
        PDBParserObj={'1': PDBParserObj}
        if num_models != 1:
            print('WARNING! When specifying just_fds==True, you cannot have multiple models. Defaulting to 1.')

    # if graphing, graph it up
    if graph==True:
        # get region info for graphing. Get here because even if user input their own
        # regions, the PDBParserObj will have it at this point. 
        region_info = PDBParserObj[1].regions_dict
        all_coords = PDBParserObj[1].all_atom_coords_by_index
        # get the CA coords for graphing.
        ca_coords = []
        for aa in all_coords:
            cur_atoms = all_coords[aa]
            for atom in cur_atoms:
                if atom == 'CA':
                    ca_coords.append(all_coords[aa][atom])
        plot_structure(ca_coords, region_info, save_path=out_path)
    else:
        # otherwise save that stuffs.
        if out_path=='':
            raise dodoException('Please specify an output path.')
        else:
            for model_number in PDBParserObj:
                if model_number==1:
                    add_mode='w'
                else:
                    add_mode='a'
                if model_number==num_models:
                    last_model=True
                else:
                    last_model=False
                save_pdb_from_PDBParserObj(PDBParserObj[model_number], out_path=out_path,
                    include_FD_atoms=include_FD_atoms, CONECT_lines=CONECT_lines,
                    add_mode=add_mode, last_model=last_model, model_num=model_number,
                    just_fds=just_fds)


def pdb_from_pdb(path_to_pdb, out_path='', mode='predicted', 
    linear_placement=False, CONECT_lines=True, include_FD_atoms=True, 
    use_metapredict=False, graph=False, verbose=True, attempts_per_region=40, 
    attempts_per_coord=2000, regions_dict=None, num_models=1, 
    beta_for_FD_IDR=False, just_fds=False):
    """
    Function to take in the path to an AF2 pdb structure and return the structure
    with modified disordered regions. 

    Parameters
    ----------
    path_to_pdb : str
        absolute path to the PDB. Include filename and extension.
    out_path : str
        Path to save the PDB to. Include filename and extension.
    mode : str
        Mode to use for building the structure. Options are 'super_compact',
        'compact', 'normal', 'expanded', 'super_expanded', 'max_expansion', 
        and 'predicted'. 'predicted' uses ALBATROSS to predict IDR end-to-end distance.
    linear_placement : bool
        whether to place the folded domains across a linear axis. 
        Default : False
    CONECT_lines : bool
        Whether to include CONECT lines in the pdb. Default is True.
    include_FD_atoms : bool
        Whether to include all atoms from the AF2 structure. Default is True.
    use_metapredict : bool
        Whether to use Metapredict to predict IDRs. Default is False.
        This is less accurate than using the actual structure but is MUCH faster.
    graph : bool
        Whether to graph the structure. Default is False.
        If you graph the output, it will not save.
    verbose : bool
        Whether to print out info about what is happening. Default is True.
    attempts_per_region : int
        Number of attempts to make per region. Default is 20.
    attempts_per_coord : int
        Number of attempts to make per coordinate. Default is 2000.
    num_models : int
        Number of times to rebuild the IDRs while holding the FDs contant.
        Default : 1
    beta_for_FD_IDR : bool
        Whether to overwrite beta vals where IDR=100 and FD=0 for 
        visualization purposes.
        Default is False
    just_fds : bool
        Whether to just dump the folded domains as individual PDBs.
        Default : False. 

    Returns
    -------
    if graph==False:
        None. Just saves the PDB.
    if graph==True:
        None. Just graphs the structure.
    """

    # don't let people waste their time trying to graph multiple models. 
    if graph==True and num_models!=1:
        num_models=1
        print('Cannot graph more than one model. Falling back to one model.')

    # make sure mode is reasonable before getting rolling.
    if mode not in list(parameters.modes.keys()):
        raise dodoException('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')    
    
    # make sure we have a PDB to open.
    if os.path.isfile(path_to_pdb)==False:
        raise dodoException('The path you provided to the pdb file does not exist.')

    # make a dict with info we need for everything.
    cur_pdb=open(path_to_pdb, 'r').read().split('\n')
    PDBParserObj = PDBParser(cur_pdb)

    # use metapredict or atom positions if user didn't specify the regions...
    if regions_dict==None:
        if use_metapredict:
            if verbose==True:
                print('Predicting regions using Metapredict V2.')
            PDBParserObj = get_fds_idrs_from_metapredict(PDBParserObj)
        else:
            if verbose==True:
                print('Predicting folded reigons, loops, and IDRs using AF2 structure information.')        
            PDBParserObj = get_fds_loops_idrs(PDBParserObj)
    else:
        PDBParserObj.regions_dict=regions_dict
        for reg in regions_dict:
            if 'folded' in reg:
                PDBParserObj.FD_coords[reg]=regions_dict[reg]
            if 'idr' in reg:
                PDBParserObj.IDR_coords[reg]=regions_dict[reg]
            if 'loop' in reg:
                PDBParserObj.FD_loop_coords[reg]=regions_dict[reg]
    # build the new structure
    if verbose==True:
        print('Building new structure.') 
        temp=''
        for reg in PDBParserObj.regions_dict:
            cur_reg = PDBParserObj.regions_dict[reg]
            if 'idr' in reg:
                reg_name = 'IDR'
            elif 'folded' in reg:
                reg_name = 'FD'
            else:
                reg_name = 'FD w/ loop(s)'
            if 'loop' in reg:
                fd = cur_reg[0]
                loops = cur_reg[1]
                loop_locs=''
                for loop in loops:
                    loop_locs+=f'{loop[0]+1} to {loop[1]+1}, '
                temp+=f"{reg_name}: FD {fd[0]+1} to {fd[1]+1} with loop(s) at {loop_locs}"
            else:
                temp+=f"{reg_name}: {cur_reg[0]+1} to {cur_reg[1]+1}, "
        print(temp[:len(temp)-2])


    if len(PDBParserObj.regions_dict)==1 and 'idr_1' in PDBParserObj.regions_dict:
        return pdb_from_sequence(PDBParserObj.sequence, out_path=out_path, mode=mode, 
            attempts_per_res=1000, attempts_per_idr=50, end_coord=(0,0,0), 
                CONECT_lines=CONECT_lines, graph=graph, num_models=num_models)


    # build the structure.
    if just_fds==False:
        PDBParserObj = build_structure(PDBParserObj, mode=mode,
                                        linear_placement=linear_placement,
                                        attempts_per_region=attempts_per_region,
                                        attempts_per_coord=attempts_per_coord,
                                        verbose=verbose, num_models=num_models,
                                        beta_for_FD_IDR=beta_for_FD_IDR)
    else:
        PDBParserObj={'1': PDBParserObj}
        if num_models != 1:
            print('WARNING! When specifying just_fds==True, you cannot have multiple models. Defaulting to 1.')

    if graph==True:
        # get region info for graphing. Get here because even if user input their own
        # regions, the PDBParserObj will have it at this point. 
        region_info = PDBParserObj[1].regions_dict
        all_coords = PDBParserObj[1].all_atom_coords_by_index
        # get the CA coords for graphing.
        ca_coords = []
        for aa in all_coords:
            cur_atoms = all_coords[aa]
            for atom in cur_atoms:
                if atom == 'CA':
                    ca_coords.append(all_coords[aa][atom])
        plot_structure(ca_coords, region_info, save_path=out_path)
    else:
        if out_path=='':
            raise dodoException('Please specify an output path.')
        else:
            for model_number in PDBParserObj:
                if model_number==1:
                    add_mode='w'
                else:
                    add_mode='a'
                if model_number==num_models:
                    last_model=True
                else:
                    last_model=False
                save_pdb_from_PDBParserObj(PDBParserObj[model_number], out_path=out_path,
                    include_FD_atoms=include_FD_atoms, CONECT_lines=CONECT_lines,
                    add_mode=add_mode, last_model=last_model, model_num=model_number,
                    just_fds=just_fds)



def pdb_from_sequence(sequence, out_path='', mode='predicted', 
    attempts_per_res=1000, attempts_per_idr=50, end_coord=(0,0,0), 
    CONECT_lines=True, graph=False, num_models=1):
    '''
    Function to generate a PDB (or graph the IDR) using a sequence as the input.

    Parameters
    ----------
    sequence : str
        Sequence to build the IDR for.
    out_path : str
        Path to save the PDB to. Include filename and extension.
    mode : str
        Mode to use for building the structure. Options are 'super_compact',
        'compact', 'normal', 'expanded', 'super_expanded', 'max_expansion', 
        and 'predicted'. 'predicted' uses ALBATROSS to predict IDR end-to-end distance.
    attempts_per_res : int
        number of attempts to place an individual residue
        Default is 1000
    attempts_per_idr : int
        Number of attempts to make entire IDR. If attempts per res hits max value,
        it resets and retries the build this number of times. Default is 50.
    end_coord : tuple
        End coordinate to use for the IDR. Default is (0,0,0).
    CONECT_lines : bool
        Whether to include CONECT lines in the pdb. Default is True.
    graph : bool
        Whether to graph the structure. Default is False.
    num_models : Int
        Number of IDR models to make. 

    Returns
    -------
    if graph==False:
        None. Just saves the PDB.
    if graph==True:
        None. Just graphs the structure.
    '''

    # don't let people waste their time trying to graph multiple models. 
    if graph==True and num_models!=1:
        num_models=1
        print('Cannot graph more than one model. Falling back to one model.')

    # make sure mode is reasonable before getting rolling.
    if mode not in list(parameters.modes.keys()):
        raise dodoException('Invalid mode specified. Please specify super_compact, compact, normal, expanded, super_expanded, or max_expansion, or predicted.')        

    # make dict to hold models
    models={}

    for model in range(1, num_models+1):
    # build the IDR
        models[model]= build_idr_from_sequence(sequence, mode=mode, end_coord=end_coord, 
            attempts_per_idr=attempts_per_idr, attempts_per_coord=attempts_per_res)

    # if user doesn't want CONECT lines, don't return them.
    if CONECT_lines==True:
        CONECT_coords=models[1]['CONECT_coords']
    else:
        CONECT_coords=None

    if graph==False:
        if out_path=='':
            raise dodoException('Please specify an output path.')
        for model_number in models:
            if model_number==1:
                add_mode='w'
            else:
                add_mode='a'
            if model_number==num_models:
                last_model=True
            else:
                last_model=False     
            idr_dict=models[model_number]
            write_pdb(idr_dict['xyz_list'], out_path, atom_indices=idr_dict['atom_indices'],
                atom_names=idr_dict['atom_names'], residue_indices=idr_dict['residue_indices'],
                residue_names=idr_dict['residue_names'], CONECT_LINES=CONECT_coords,
                add_mode=add_mode, last_model=last_model, model_num=model_number)
    else:
        plot_structure(models[1]['xyz_list'], {'idr': [0, len(sequence)]})




