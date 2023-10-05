# rewriting some PDB stuff. Admitedly a large amount of it could be done with dependencies
# but where's the fun in that? This is some gluten-free, chemical-free, alkaline,
# all natural good code. None of that chemical-filled dependency nonsense. 

import numpy as np
from dodo.parameters import AADICT_3_to_1, AADICT
from dodo.dodo_exceptions import dodoPDBException

# in case array instead np.array is returned.. ugh. Will fix this problem later but here's a bandaid 
class array(list):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.array = np.array(self)

def format_string(input_string, total_length):
    """
    Format a string to a specified total length by trimming or padding.

    Parameters
    ----------
    input_string : str 
        The input string to be formatted.
    total_length :int 
        The desired total length of the resulting string.

    Returns
    --------
    str :
        The formatted string.
    """
    if len(input_string) > total_length:
        # Truncate the string if it's longer than the specified length
        formatted_string = input_string[:total_length]
    else:
        # Pad the string with spaces if it's shorter than the specified length
        formatted_string = input_string.ljust(total_length)

    return formatted_string


class PDBParser:
    '''
    A PDB file parser because none exist in the world so I made one.
    Not really, but I wanted weird stuff so I made this. 
    '''
    def __init__(self, pdb_lines):
        self.sequence = ""
        self.sequence_3aa_by_index = {}
        self.all_atom_coords_by_index = {}
        self.index_to_3aa = {}
        self.beta_vals_by_index = {}
        self.seqres = []
        self.chain = ''
        self.ter=[]
        self.regions_dict={}
        # separating out now.
        self.FD_loop_coords={}
        self.FD_coords={}
        self.IDR_coords={}
        self.number_atoms=0
        self.parse_pdb_lines(pdb_lines)


    def parse_pdb_lines(self, pdb_lines):
        for line in pdb_lines:
            if line.startswith("ATOM"):
                if self.chain == '':
                    self.chain = line[21]
                beta_val=line[60:66].strip()
                atom_index=int(line[6:11])-1
                aa = line[17:20].strip()
                res_index = int(line[22:26])-1
                atom_type = line[12:16].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                
                if res_index not in self.all_atom_coords_by_index:
                    self.index_to_3aa[res_index] = aa
                    self.sequence+=AADICT_3_to_1[aa]
                    self.all_atom_coords_by_index[res_index] = {}
                    self.sequence_3aa_by_index[res_index]={}
                    self.beta_vals_by_index[res_index] = {}

                self.all_atom_coords_by_index[res_index][atom_type] = (x, y, z)
                self.sequence_3aa_by_index[res_index][atom_type]=aa
                self.beta_vals_by_index[res_index] = beta_val

            if line.startswith("SEQRES"):
                self.seqres.append(line)
            if line.startswith('TER'):
                self.ter.append(line)

    def __str__(self):
        return f"PDBParser: {len(self.sequence)} amino acids"





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
              sequence='', 
              CONECT_LINES=None):
                 
    
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

    # just list of cooridnates for the alpha carbons of the IDRs. Otherwise they don't show up connected in VMD. 
    if CONECT_LINES != None:
        CONECT=[]
        for pair in CONECT_LINES:
            CONECT.append(f'CONECT{str(pair[0]).rjust(5)}{str(pair[1]).rjust(5)}')

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


    if CONECT_LINES!=None:
        fh.write(f'\n')
        for line in CONECT:
            fh.write(f'{line}\n')

    fh.write('ENDMDL\n')
    fh.write('END')


    # close file handle
    fh.close()

        
                 
def save_pdb_from_PDBParserObj(PDBParserObj, out_path, 
    include_FD_atoms, CONECT_lines):
    '''
    Function for saving pdbs from PDBParserObj.
    
    Parameters
    ----------
    PDBParserObj : PDBParser object
        Object with all the info needed to make a pdb.
    out_path : str
        Path to save pdb to.
    include_FD_atoms : bool
        Whether to include all atoms from the AF2 structure.
    CONECT_lines : bool
        Whether to include CONECT lines in the pdb.

    Returns
    -------
    None. Just saves the PDB.
    '''
    # make lists based on what to include...
    xyz_list=[]
    residue_names=[]
    residue_indices=[]
    atom_indices=[]
    beta_vals=[]
    atom_names=[]
    CONECT_COORDS=[]

    # get all the coords.
    all_coords = PDBParserObj.all_atom_coords_by_index

    # let user not use all atom from FDs if they want. 
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
                # only conect N-CA-C because this is for the IDR. 
                if atom_count < PDBParserObj.number_atoms:
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

    # if user doesn't want CONECT lines, nuke them. 
    if CONECT_lines==False:
        CONECT_COORDS=None

    # write pdb
    write_pdb(xyz_list, out_path, atom_indices=atom_indices,
            atom_names=atom_names, residue_indices=residue_indices,
            residue_names=residue_names, beta=beta_vals,CONECT_LINES=CONECT_COORDS)

