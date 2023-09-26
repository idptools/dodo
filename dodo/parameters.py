"""
Various parameter dictionaries for PDB file creation

"""
          
AADICT = {"A":"ALA", 
          "C":"CYS",
          "D":"ASP",
          "E":"GLU",
          "F":"PHE",
          "G":"GLY",
          "H":"HIS",
          "I":"ILE",
          "K":"LYS",
          "L":"LEU",
          "M":"MET",
          "N":"ASN",
          "P":"PRO",
          "Q":"GLN",
          "R":"ARG",
          "S":"SER",
          "T":"THR",
          "V":"VAL",
          "W":"TRP",
          "Y":"TYR"}

AADICT_3_to_1 = {'ALA': 'A',
                'CYS': 'C',
                'ASP': 'D',
                'GLU': 'E',
                'PHE': 'F',
                'GLY': 'G',
                'HIS': 'H',
                'ILE': 'I',
                'LYS': 'K',
                'LEU': 'L',
                'MET': 'M',
                'ASN': 'N',
                'PRO': 'P',
                'GLN': 'Q',
                'ARG': 'R',
                'SER': 'S',
                'THR': 'T',
                'VAL': 'V',
                'TRP': 'W',
                'TYR': 'Y'}  

# Pre-defined coordinates for atoms of the 20 amino acids shown in AF2 structures
# Warning! These might not be accurate TBH. 
AMINO_ACID_ATOMS_AF2 = {
    'M': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.24, -0.94, 0.0],
        'CB': [-0.32, -1.33, 0.88],
        'O': [2.27, -0.49, 0.0],
        'CG': [0.03, -2.76, 0.0],
        'SD': [-0.41, -4.33, -0.88],
        'CE': [1.78, -5.36, -1.76]
    },
    'K': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.34, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.04, -2.76, -0.89],
        'CD': [-0.51, -3.74, -0.88],
        'CE': [1.06, -4.74, -0.89],
        'NZ': [1.06, -5.98, -0.88]
    },
    'A': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.25, -0.94, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.26, -0.49, 0.0]
    },
    'P': {
        'N': [0.0, 1.19, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.24, -0.94, 0.0],
        'CB': [-0.42, -1.33, 0.0],
        'O': [2.25, -0.49, 0.0],
        'CG': [-0.48, -2.76, 0.0],
        'CD': [0.56, -3.74, 0.0]
    },
    'S': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.49, 0.0],
        'OG': [0.02, -2.67, 0.0]
    },
    'N': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'ND2': [0.02, -3.77, 0.0],
        'OD1': [1.30, -3.77, 0.0]
    },
    'G': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.25, -0.94, 0.0],
        'O': [2.26, -0.49, 0.0]
    },
    'F': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD1': [-1.30, -3.77, 0.0],
        'CD2': [1.34, -3.77, 0.0],
        'CE1': [-2.61, -4.79, 0.0],
        'CE2': [2.65, -4.79, 0.0],
        'CZ': [0.02, -5.80, 0.0]
    },
    'L': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.25, -0.94, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.26, -0.49, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD1': [-1.30, -3.77, 0.0],
        'CD2': [1.34, -3.77, 0.0]
    },
    'E': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD': [1.51, -3.74, 0.0],
        'OE1': [1.36, -5.09, 0.0],
        'OE2': [2.75, -3.40, 0.0]
    },
    'I': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.25, -0.94, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.26, -0.49, 0.0],
        'CG1': [0.02, -2.76, 0.88],
        'CG2': [0.02, -2.76, -0.88],
        'CD1': [0.02, -4.21, 0.88]
    },
    'Q': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD': [1.51, -3.74, 0.0],
        'NE2': [1.36, -5.09, 0.0],
        'OE1': [2.75, -3.40, 0.0]
    },
    'W': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD1': [-1.30, -3.77, -0.88],
        'CD2': [1.34, -3.77, -0.88],
        'CE2': [2.65, -4.79, -0.88],
        'CE3': [0.02, -5.80, -0.88],
        'NE1': [-1.30, -5.15, -1.77],
        'CH2': [1.34, -5.15, -1.77],
        'CZ2': [3.04, -6.19, -1.77],
        'CZ3': [0.42, -7.20, -1.77]
    },
    'H': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD2': [0.02, -4.21, 0.0],
        'ND1': [1.31, -3.40, 0.0],
        'CE1': [1.31, -4.68, 0.0],
        'NE2': [0.02, -5.50, 0.0]
    },
    'C': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.25, -0.94, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.26, -0.49, 0.0],
        'SG': [0.02, -2.76, 0.0]
    },
    'V': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.25, -0.94, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.26, -0.49, 0.0],
        'CG1': [-0.32, -2.76, 0.88],
        'CG2': [-0.32, -2.76, -0.88]
    },
    'Y': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD1': [-1.30, -3.77, 0.0],
        'CD2': [1.34, -3.77, 0.0],
        'CE1': [-2.61, -4.79, 0.0],
        'CE2': [2.65, -4.79, 0.0],
        'OH': [0.02, -6.23, 0.0],
        'CZ': [-1.29, -5.82, 0.0]
    },
    'T': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG2': [-0.32, -2.76, 0.88],
        'OG1': [-0.32, -2.76, -0.88]
    },
    'D': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'OD1': [0.02, -4.21, 0.0],
        'OD2': [1.31, -3.40, 0.0]
    },
    'R': {
        'N': [0.0, 1.22, 0.0],
        'CA': [0.0, 0.0, 0.0],
        'C': [1.26, -0.93, 0.0],
        'CB': [-0.32, -1.33, 0.0],
        'O': [2.27, -0.48, 0.0],
        'CG': [0.02, -2.76, 0.0],
        'CD': [1.51, -3.74, 0.0],
        'NE': [1.36, -5.09, 0.0],
        'NH1': [2.75, -3.40, 0.0],
        'NH2': [0.02, -6.25, 0.0],
        'CZ': [1.31, -6.02, 0.0],
        'OXT': [2.27, -1.44, 0.0]
    }
}