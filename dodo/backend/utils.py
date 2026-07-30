'''
random useful utilities
'''
import os
import numpy as np


def write_dumby_pdb(CA_coordinates, save_path, max_res=9999):
    '''
    function to write a dumby pdb so we can 
    easily visualize things we do to CA coordinates.

    CA_coordinates: np.array, shape=(N, 3)
    save_path: str, path to save pdb file
    '''
    # open file
    atom_per_res = int(len(CA_coordinates)/ max_res)+1
    cur_res=1
    with open(save_path, 'w') as f:
        # iterate over CA coordinates
        for i, CA in enumerate(CA_coordinates):
            if i % atom_per_res == 0:
                cur_res+=1
            # write atom line
            f.write(f'ATOM  {i+1:5}  CA  ALA A{cur_res:4}    {CA[0]:8.3f}{CA[1]:8.3f}{CA[2]:8.3f}  1.00  0.00           C\n')
        # write end line
        f.write('END\n')
    # close file
    f.close()

def amino_acid_3_to_1(three_letter_code):
    '''
    converts a 3 letter amino acid code to a 1 letter amino acid code.
    '''
    amino_acid_dict = {
        'ALA': 'A',
        'ARG': 'R',
        'ASN': 'N',
        'ASP': 'D',
        'CYS': 'C',
        'GLU': 'E',
        'GLN': 'Q',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LEU': 'L',
        'LYS': 'K',
        'MET': 'M',
        'PHE': 'F',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'TRP': 'W',
        'TYR': 'Y',
        'VAL': 'V'
    }
    return amino_acid_dict[three_letter_code]


def write_chain_to_pdb(chain_obj, save_path,
                       renumber_atoms=True):
    '''
    function to write a complex to a pdb file.

    Parameters
    ---------
    chain_obj : Chain
        Chain object to write to pdb file
    
    save_path : str
        path to save pdb file

    renumber_atoms : bool
        if True, renumber atoms to start at 1 and increment by 1
        for each atom. if False, use the atom_id from the chain_obj.
    
    Returns
    -------
    None
    '''
    # make sure path to save exists
    if os.path.exists(os.path.dirname(save_path)) == False:
        raise ValueError(f'Path to save pdb file does not exist: {save_path}')
    
    atom_number = 1

    # open file
    with open(save_path, 'w') as f:
        # iterate over domains
        for domain in chain_obj.domains.values():
            # iterate over monomers
            for monomer in domain.monomers.values():
                # iterate over atoms
                for atom in monomer.atoms.values():
                    if renumber_atoms:
                        atom.atom_id = atom_number
                        atom_number += 1
                    # Format each component according to PDB specifications
                    atom_line = (
                        f'ATOM  '                        # cols 1-6
                        f'{atom.atom_id:5} '            # cols 7-11 (+space)
                        f'{atom.atom_name:<4}'          # cols 13-16
                        f' {monomer.monomer_name:3} '    # cols 18-20
                        f'{chain_obj.chain_id}'         # col 22
                        f'{monomer.monomer_number:4}'   # cols 23-26
                        f'    '                         # cols 27-30
                        f'{float(atom.x_coord):8.3f}'   # cols 31-38
                        f'{float(atom.y_coord):8.3f}'   # cols 39-46
                        f'{float(atom.z_coord):8.3f}'   # cols 47-54
                        f'{float(monomer.occupancy):6.2f}'  # cols 55-60
                        f'{float(monomer.b_factor):6.2f}'   # cols 61-66
                        f'          '                   # cols 67-76
                        f'{atom.element:>2}\n'          # cols 77-78
                    )
                    f.write(atom_line)
        f.write('END\n')
    # close file
    f.close()




def write_multiple_conformations_to_pdb(conformations, save_path):
    '''
    function to write multiple conformations to a pdb file.

    Parameters
    ---------
    conformations : Np.array 
        Np.array of shape N,L,3 where N is the number of conformations,
        L is the length of the IDR and 3 is the X, Y, Z coordinates.

    save_path : str
        path to save pdb file

    Returns
    -------
    None
    '''
    pass

