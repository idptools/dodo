'''
random useful utilities
'''
import os
import numpy as np


def write_dumby_pdb(CA_coordinates, save_path):
    '''
    function to write a dumby pdb so we can 
    easily visualize things we do to CA coordinates.

    CA_coordinates: np.array, shape=(N, 3)
    save_path: str, path to save pdb file
    '''
    # open file
    with open(save_path, 'w') as f:
        # iterate over CA coordinates
        for i, CA in enumerate(CA_coordinates):
            # write atom line
            f.write(f'ATOM  {i+1:5}  CA  ALA A{i+1:4}    {CA[0]:8.3f}{CA[1]:8.3f}{CA[2]:8.3f}  1.00  0.00           C\n')
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
