# various tools
import metapredict as meta
from dodo.parameters import AADICT_3_to_1, AADICT
from dodo.dodo_exceptions import dodoException, dodoPDBException
import matplotlib.pyplot as plt


# takes in lines from a pbd file and returns a dict.
def af2_lines_to_structure_dict(af2_pdb_lines):
    '''
    function to return a dictionary of the AF2 structure info.

    parameters
    ----------
    af2_pdb_lines : list
        the lines from the AF2 pdb file.

    Returns
    -------
    dict
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
    '''
    # lists to hold info (and a str for sequence)
    all_atom_coords=[]
    ca_coords=[]
    all_atoms_types=[]
    sequence=''
    ca_three_letter_seq=[]
    all_atom_three_letter_seq=[]
    coords_by_aa={}
    atom_coords_per_aa=[]
    start_resnum=1
    for line in af2_pdb_lines:
        if line != '':
            if line[:4]=='ATOM':
                resname=line[17:20]
                chain_id=line[21]                
                restype=line[13:16]
                resnum=int(line[22:26])
                x=float(line[30:38])
                y=float(line[38:46])
                z=float(line[46:54])
                if resnum not in coords_by_aa:
                    coords_by_aa[resnum]=[[x,y,z]]
                else:
                    cur_coords=coords_by_aa[resnum]
                    cur_coords.append([x,y,z])
                    coords_by_aa[resnum]=cur_coords
                all_atom_coords.append([x,y,z])
                atom_coords_per_aa.append([x,y,z])
                all_atom_three_letter_seq.append(resname)
                all_atoms_types.append(restype)
                if restype=='CA ' or restype==' CA':
                    ca_coords.append([x,y,z])
                    ca_three_letter_seq.append(resname)
                    sequence+=AADICT_3_to_1[resname]



    # return lists
    return {'all_atom_coords':all_atom_coords, 'ca_coords':ca_coords, 
        'all_atoms_types':all_atoms_types, 'sequence':sequence, 
        'ca_three_letter_seq':ca_three_letter_seq, 
        'all_atom_three_letter_seq':all_atom_three_letter_seq,
        'coords_by_aa':coords_by_aa}



def plot_simple(all_coordinates, highlight=[], dims=[10,10]):
    '''
    simplified version of plot_structure. Doesn't require any coords
    info on fd vs. idr vs. loop.

    Parameters
    ----------
    all_coordinates : list
        a list of the coordinates of all atoms in the structure.
    highlight : list, optional
        a list of lists of coordinates to highlight. The default is [].
    dims : list, optional
        the dimensions of the plot. The default is [10,10].

    Returns
    -------
    None
        Just plots the graph and shows it.
    '''

    # Create a 3D plot
    fig = plt.figure(figsize=dims)
    ax = fig.add_subplot(111, projection='3d')
    if highlight != []:
        for coord in highlight:
            all_coordinates.remove(coord)
    x_c, y_c, z_c = zip(*all_coordinates)
    ax.scatter(x_c, y_c, z_c, s=70, depthshade=True)
    if highlight != []:
        hi_xc, hi_yc, hi_zc=zip(*highlight)
        for hi in range(0, len(highlight)):
            curx=hi_xc[hi]
            cury=hi_yc[hi]
            curz=hi_zc[hi]
            ax.scatter(curx, cury, curz, s=70, depthshade=True, color='red')
    
    # Set make sure axes all same size so structure looks good. 
    all_coords=[]
    all_coords.extend(x_c)
    all_coords.extend(y_c)
    all_coords.extend(z_c)
    ax.set_ylim(min(all_coords),max(all_coords))
    ax.set_xlim(min(all_coords),max(all_coords))
    ax.set_zlim(min(all_coords),max(all_coords))
    
    # get rid of annoying stuffs
    all_ax=[ax.xaxis, ax.yaxis, ax.zaxis]
    for curax in all_ax:
        curax.set_pane_color((1.0, 1.0, 1.0, 0.0))
        curax.set_ticks([])
        curax._axinfo["grid"]['color'] =  (1,1,1,0)
        curax._axinfo['axisline']['linewidth']=0
        curax._axinfo['grid']['linewidth']=0
        curax._axinfo['juggled']=(0,0,0)

    plt.show()

def plot_structure(all_coordinates, fd_idr_loop_dict, dims=[10,10], save_path='', 
    fd_color='blue', IDR_color='orange', loop_color='red'):
    '''
    Quick way to view a structure using matplotlib.

    Parameters
    ----------
    all_coordinates : list
        a list of the coordinates of all atoms in the structure.
    fd_idr_loop_dict : dict
        a dictionary of the structure info. keys have 'idr'
        in them for idrs, 'folded' in them for FDs, and
        'loop' in them for loops. Values are lists 
        that contain the coordinates in the sequence for each
        different 'feature' type.
    dims : list, optional
        the dimensions of the plot. The default is [10,10].
    save_path : string, optional
        if you want to save the plot, specify a path here.
        The default is ''. If it is just not specified, just shows plot.
    fd_color : string, optional
        the color to plot the folded domains. The default is 'blue'.
    IDR_color : string, optional
        the color to plot the IDRs. The default is 'orange'.
    loop_color : string, optional
        the color to plot the loops. The default is 'red'.


    Returns
    -------
    None
        Shows plot of the structure.
    '''
    # Create a 3D plot
    fig = plt.figure(figsize=dims)
    ax = fig.add_subplot(111, projection='3d')

    # break up coords into each 'part'
    fd_coords=[]
    idr_coords=[]
    loop_coords=[]
    for fd_idr_loop in fd_idr_loop_dict:
        if 'idr' in fd_idr_loop:
            cur_coords=fd_idr_loop_dict[fd_idr_loop]
            idr_coords.append(all_coordinates[cur_coords[0]:cur_coords[1]])
        elif 'folded' in fd_idr_loop:
            cur_coords=fd_idr_loop_dict[fd_idr_loop]
            fd_coords.extend(all_coordinates[cur_coords[0]:cur_coords[1]])
        elif 'loop' in fd_idr_loop:
            fd_loop_coords=fd_idr_loop_dict[fd_idr_loop]
            cur_fd_coords=fd_loop_coords[0]
            cur_loop_in_fd_coords=fd_loop_coords[1]
            all_fd_coords=all_coordinates[cur_fd_coords[0]:cur_fd_coords[1]]
            for loop in cur_loop_in_fd_coords:

                cur_loop_coords=all_coordinates[loop[0]:loop[1]]
                for coord in cur_loop_coords:
                    all_fd_coords.remove(coord)
                loop_coords.append(cur_loop_coords)
            fd_coords.extend(all_fd_coords)
        else:
            raise dodoException('Error - coordinate input did not seem to be an IDR, FD, or loop.')

    # extract coords for each.
    if fd_coords != []:
        fd_x_coords, fd_y_coords, fd_z_coords = zip(*fd_coords)
        if fd_color != 'blue':
            ax.scatter(fd_x_coords, fd_y_coords, fd_z_coords, s=70, depthshade=True, color=fd_color,)
        else:
            ax.scatter(fd_x_coords, fd_y_coords, fd_z_coords, s=70, depthshade=True)

    if idr_coords != []:
        if len(idr_coords) > 1:
            for cur_coord in idr_coords:
                idr_x_coords, idr_y_coords, idr_z_coords = zip(*cur_coord)
                # plt the IDR
                ax.plot(idr_x_coords, idr_y_coords, idr_z_coords, color=IDR_color, linewidth=2) 
        else:   
            idr_x_coords, idr_y_coords, idr_z_coords = zip(*idr_coords)
            ax.plot(idr_x_coords, idr_y_coords, idr_z_coords, color=IDR_color, linewidth=2) 

    if loop_coords != []:
        if len(loop_coords) > 1:
            for cur_coord in loop_coords:
                loop_x_coords, loop_y_coords, loop_z_coords = zip(*cur_coord)  
                ax.plot(loop_x_coords, loop_y_coords, loop_z_coords, color=loop_color, linewidth=2)              
        else:
            loop_x_coords, loop_y_coords, loop_z_coords = zip(*loop_coords[0])
            # plot the loops
            ax.plot(loop_x_coords, loop_y_coords, loop_z_coords, color=loop_color, linewidth=2)



    # Set make sure axes all same size so structure looks good. 
    all_coords_x, all_coords_y, all_coords_z = zip(*all_coordinates)
    all_coords=[]
    all_coords.extend(all_coords_x)
    all_coords.extend(all_coords_y)
    all_coords.extend(all_coords_z)
    ax.set_ylim(min(all_coords),max(all_coords))
    ax.set_xlim(min(all_coords),max(all_coords))
    ax.set_zlim(min(all_coords),max(all_coords))
    
    # get rid of annoying stuffs
    all_ax=[ax.xaxis, ax.yaxis, ax.zaxis]
    for curax in all_ax:
        curax.set_pane_color((1.0, 1.0, 1.0, 0.0))
        curax.set_ticks([])
        curax._axinfo["grid"]['color'] =  (1,1,1,0)
        curax._axinfo['axisline']['linewidth']=0
        curax._axinfo['grid']['linewidth']=0
        curax._axinfo['juggled']=(0,0,0)

    # save or show
    if save_path=='':
        # Show the plot
        plt.show()
    else:
        plt.savefig(save_path)


#
# ----------------------------------------------------------------------------
#
# function to write PDBs. Shoutout to Alex Holehouse for writing this up.

def write_pdb(xyz_list,
              outname,
              box_dims=[500,500,500],
              atom_start=1,
              atom_indices = None,
              atom_name='CA',
              atom_names=None,
              residue_start=1,
              residue_indices=None,
              residue_name='GLY',
              residue_names=None,
              chain_id='A',
              chain_ids=None,
              beta=None,
              one_atom_per_residue=True, 
              include_seq_res=False,
              sequence=''):
                 
    
    """
    Create a new pdb file based on parsing ATOM lines ONLY. Skips HETROATOMS for now...
    
    This parser follows the PDB specification to the letter.
        
    Example line:
    ATOM   2209  CB  TYR A 299       6.167  22.607  20.046  1.00  8.12           C
    00000000011111111112222222222333333333344444444445555555555666666666677777777778
    12345678901234567890123456789012345678901234567890123456789012345678901234567890
    ATOM line format description from
    http://deposit.rcsb.org/adit/docs/pdb_atom_format.html:
    COLUMNS        DATA TYPE       CONTENTS
    --------------------------------------------------------------------------------
    1 -  6         Record name     "ATOM  "
    7 - 11         Integer         Atom serial number.
    13 - 16        Atom            Atom name.
    17             Character       Alternate location indicator.
    18 - 20        Residue name    Residue name.
    22             Character       Chain identifier.
    23 - 26        Integer         Residue sequence number.
    27             AChar           Code for insertion of residues.
    31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)       Occupancy (Default = 1.0).
    61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
    73 - 76        LString(4)      Segment identifier, left-justified.
    77 - 78        LString(2)      Element symbol, right-justified.
    79 - 80        LString(2)      Charge on the atom.
        
    Simple constructor that takes in a structured data and writes out a valid PDB file.
    
    Right now the implementation doesn't allow a ton of variability but we can and will build on this going forward as
    new features are needed.
    
    Parameters
    --------------
    xyz_list : list, np.array or other iterable type
        This contains a nx3 matrix where each element contains three positions that define the x,y and z positions in 
        the PDB file.

    outname : str
        Name and location of the PDB file to be written

    box_dims : list or tuple of len(3)
        A list or tuple that defines the x, y, and z sizes of the PDB box. Default = [500,500,500]

    atom_start : int
        Defines the index of the first atom to be written. Defaults = 1

    atom_indices : list
        If provided, must be a list that is equal to the length of the xyz_list, and each position
        in the XYZ list will be mapped to the corresponding atom_index here. This over-rides atom_start.

    atom_name : str
        If provide, means all atoms have the same name and this name is defined by atom_name.
        Default = 'CA'.

    atom_names : list or None
        If provided, over-writes atom_name and assigns a name to each atom. Must be a list that
        equals in length the length of the xyz_list. Can be a max of 4 characters. Default = None

    residue_start : int 
        Defines the index of the first residueatom to be written. Defaults to 1

    residue_indices : list
        If provided, must be a list that is equal to the length of the xyz_list, and each position
        in the XYZ list will be mapped to the corresponding residue index in here. This over-rides residue_start.

    residue_name : str
        If provide, means all residues have the same name and this name is defined by residue_name.
        Can be a max of 3 chars. Default = 'GLY'.

    residue_names : list or None
        If provided, over-writes residue_name and assigns a name to each residue. Must be a list that
        equals in length the length of the xyz_list (assuming one_atom_per_residue is True, which for
        now it must be). Default = None.
        
    chain_id : str
        Defies the ID used for defining the chain. Must be a single character. Default = 'A'.

    chain_ids : list 
        If provided, must be a list that is equal to the length of the xyz_list, and each position
        in the XYZ list will be mapped to the corresponding chain index in here. This over-rides chain_id.

    beta : list or None
        If provided, assigns a beta value to every atom in the xyz_list. Default = None


    Returns
    --------------
    None
        Does not return anything, but a valid PDB file is written to the location defined by
        outname. The pdb extension is not added automatically.


    ## TO DO:
    
    (1) Add option to center in the box
    (2) Add support for one-atom-per-residue being false
    
    """

    # ########################################################################################
    # ###
    # START OF LOCAL FUNCTIONS
    # ###
    
    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #
    def numeric_padder(v, length, rounder=5, head_padding=False):
        """
        Local function for building numerical values that conform to the
        PDB specification.
        
        Specifically, this converts v (a value) into a string that's equal
        in length to length.

        Parameters
        -------------
        v : float or int
            Numeric type input

        length : int
            The precise length that a string of v should be

        rounder : int
            A passed value that defines the initial degree of rounding
            to be used. Probably leave this as is unless there's a good
            reason to change --> is used here as a recursive index tracker.
            Default = 5.

        head_padding : bool
            Flag which, if set to true, means whitespace is added at the
            front of the string. If set to False, whitespaceis added at 
            the end of the string. Default = False.

        Returns
        ----------
        str
            Returns a string of length $length that maps to the value 
            passed in v.


        """
        raw = str(v)


        # correct length!
        if len(raw) == length:
            return raw

        # too short - pad with whitespace
        if len(raw) < length:
            if head_padding:                    
                return " "*(length-len(raw)) + raw 
            else:
                return raw + " "*(length-len(raw))
                    

        # too long - recursive
        if rounder > 0:
            new_v = round(v, rounder)
            return numeric_padder(new_v, length, rounder-1)
        elif rounder == 0:
            new_v = int(round(v, 0))
            return numeric_padder(new_v, length, -1)
        else:
            raise dodoException('Could not round this')

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    #
    #
    def string_padder(s, length, head_padding=False):

        """
        Local function for building strings that conform to the
        PDB specification.
        
        Specifically, this converts an input string (s) into a string that's 
        equal in length to the input variable length.

        Parameters
        -------------
        s : str 
            Input string

        length : int
            The precise length that a string of s should be

        head_padding : bool
            Flag which, if set to true, means whitespace is added at the
            front of the string. If set to False, whitespaceis added at 
            the end of the string. Default = False.

        Returns
        ----------
        str
            Returns a string of length $length that maps to the value 
            passed in s.


        """

        # correct length!
        if len(s) == length:
            return s

        # too short - pad with whitespace
        if len(s) < length:
            if head_padding:                    
                return " "*(length-len(s)) + s 
            else:
                return s + " "*(length-len(s))
        else:            
            # truncate. Brutal.
            return s[0:length]
                
    # ###
    # END OF LOCAL FUNCTIONS
    # ###
    # ########################################################################################

                
    if one_atom_per_residue is not True:
        raise dodoPDBException('Right now, one-atom-per-residue is required for PDB writing... sorry')

    # get number of atom-lines the final PDB file will have
    n_atom_lines = len(xyz_list)

    # parse atom names
    if atom_names is not None:
        if len(atom_names) != n_atom_lines:
            raise dodoPDBException('Number of elements in atom_names must equal number of elements in xyz_list')
    else:
        atom_names = [atom_name]*n_atom_lines

    ## parse residue names
    if residue_names is not None:
        if len(residue_names) != n_atom_lines:
            raise dodoPDBException('Number of elements in residue_names must equal number of elements in xyz_list')
    else:
        residue_names = [residue_name]*n_atom_lines

    # parse atom indices
    if atom_indices is not None:
        if len(atom_indices) != n_atom_lines:
            raise dodoPDBException('Number of elements in atom_indices must equal number of elements in xyz_list')
    else:
        atom_indices = range(atom_start, n_atom_lines+1)

    # parse residue indices
    if residue_indices is not None:
        if len(residue_indices) != n_atom_lines:
            raise dodoPDBException('Number of elements in residue_indices must equal number of elements in xyz_list')
    else:
        residue_indices = range(residue_start, n_atom_lines+1)

    # parse chain IDs
    if chain_ids is not None:
        if len(chain_ids) != n_atom_lines:
            raise dodoPDBException('Number of elements in chain_ids must equal number of elements in xyz_list')

        # check...
        for i in chain_ids:
            if len(i) > 1:
                raise dodoPDBException('Sadly all chain_ids must be single-character identifiers, and the identifier {i} has more than one character!')
    else:
        chain_ids = [chain_id]*n_atom_lines
        
    # parse beta
    if beta is not None:
        if len(beta) != n_atom_lines:
            raise dodoPDBException('Number of elements in chain_ids must equal number of elements in xyz_list')
    else:
        beta = [1]*n_atom_lines

    # open the filehandle
    fh = open(outname, 'w')

    # write the CRYSTAL line. The box dimensions need to be 7 chars long each
    box_x = numeric_padder(box_dims[0], 7)
    box_y = numeric_padder(box_dims[1], 7)
    box_z = numeric_padder(box_dims[2], 7)

    fh.write(f'CRYST1  {box_x}  {box_y}  {box_x}  90.00  90.00  90.00 P 1\n')
    fh.write('MODEL    1\n')

    # define widths of PDB columns

    # ATOM_RECORD LINE
    ATOM_SERIAL_NUMBER_WIDTH = 5
    ATOM_NAME_WIDTH = 4
    # ALTERNATE LOCATION INDICATOR (17)
    RESNAME_WIDTH=3
    CHAIN_IDENTIFIER_WITHD=1
    RES_NUMBER_WIDTH=4
    # INSERTION (27)
    XWIDTH=8
    YWIDTH=8
    ZWIDTH=8
    OCCUPANCY_WIDTH=6
    BETA_WIDTH=6
    FINAL_PAD = ' '*8
    SEQRES_WIDTH=4
    RES_PER_SEQRES=13
    SEQRES_PAD = ' '*10
    SEQRES_LEN_PAD=5


    # initialize
    idx = 0
    res_idx =  residue_start
    atom_idx = atom_start
    
    if include_seq_res:
        if sequence=='':
            raise dodoPDBException('Need to provide a sequence if you want to add SEQRES lines. There are other ways to do this but... sorry!')

        # figure out how many lines. 
        n_seqres_lines = int(len(sequence)/RES_PER_SEQRES)
        if n_seqres_lines*RES_PER_SEQRES != len(sequence):
            n_seqres_lines+=1
        for seqresline in range(1, n_seqres_lines+1):
            seq_res_line=''
            cur_aas=sequence[(seqresline-1)*RES_PER_SEQRES:seqresline*RES_PER_SEQRES]
            seq_res_line_num=f'{seqresline}'
            for i in range(0, 4-len(seq_res_line_num)):
                seq_res_line_num=seq_res_line_num+' '
            seq_res_line+=f'SEQRES{seq_res_line_num} A'
            seq_res_seq_len=f'{len(sequence)}'
            for i in range(0, 5-len(seq_res_seq_len)):
                seq_res_seq_len=' '+seq_res_seq_len
            seq_res_line+=seq_res_seq_len
            seq_res_line+='  '
            for aa in cur_aas:
                seq_res_line+=f'{AADICT[aa]} '
            seq_res_line+='         '
            fh.write(seq_res_line+'\n')

    for atom in xyz_list:

        # extract specific data to write
        atom_idx = atom_indices[idx]
        atom_name = atom_names[idx]
        residue_name = residue_names[idx]
        chain_id = chain_ids[idx]
        res_idx = residue_indices[idx]
        b = beta[idx]

        atom_line = 'ATOM  '
        
        atom_line = atom_line + numeric_padder(atom_idx, ATOM_SERIAL_NUMBER_WIDTH, rounder=0, head_padding=True)
        atom_line = atom_line + ' '  # note space here at 12
        atom_line = atom_line + string_padder(atom_name, ATOM_NAME_WIDTH, head_padding=True)
        atom_line = atom_line + ' '  # alternate location indicate, col 17
        atom_line = atom_line + string_padder(residue_name, RESNAME_WIDTH)
        atom_line = atom_line + ' ' + string_padder(chain_id, 1) # note single padding space here (21)
        atom_line = atom_line + numeric_padder(res_idx, RES_NUMBER_WIDTH, rounder=0, head_padding=True)
        atom_line = atom_line + ' '  # insertion for res (col 27)
        atom_line = atom_line + ' '*3  # gap between col 27 and col 31
        
        
        if len(atom_line) != 27:
            dodoPDBException(f'Error writing PDB file. This is a bug [col 27, line so far:{atom_line}]')
            
        atom_line = atom_line + numeric_padder(atom[0], XWIDTH, rounder=4, head_padding=True)
        atom_line = atom_line + numeric_padder(atom[1], YWIDTH, rounder=4, head_padding=True)
        atom_line = atom_line + numeric_padder(atom[2], ZWIDTH, rounder=4, head_padding=True)

        atom_line = atom_line + " " * OCCUPANCY_WIDTH # keep in case we add data here later
        atom_line = atom_line + numeric_padder(b, BETA_WIDTH, rounder=2,  head_padding=True) # keep in case we add data here later
        atom_line = atom_line + FINAL_PAD + '\n'
        fh.write(atom_line)

        idx = idx + 1

    fh.write('ENDMDL\n')
    fh.write('END')

        
                 
        


