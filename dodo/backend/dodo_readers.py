'''
class to read in .cif or .pdb files to populate the dodo_structures.Complex object
'''
import os
import numpy as np
from itertools import islice
from itertools import islice

from dodo.backend.dodo_structures import Atom, Monomer, Domain, Chain, Complex

class Reader:
    def __init__(self, filepath, complex_id=None):
        self.filepath = filepath
        if complex_id is None:
            complex_id = os.path.basename(filepath).split('.')[0]
        self.complex=Complex(complex_id=complex_id)
        
    def read(self):
        """
        Reads in the lines.
        """
        # make sure file exists. 
        if not os.path.exists(self.filepath):
            raise FileNotFoundError(f'File {self.filepath} not found')
        
        with open(self.filepath, 'r') as file:
            lines = file.read().split('\n')
        file.close()
        return lines
    
    def parse_cif(self):
        """
        Parses the CIF file and extracts relevant information.
        """
        # Read and split lines once
        lines = self.read()

        # Find all 'loop_' occurrences
        loop_ind_ids = [i + 1 for i, line in enumerate(lines) if line.strip() == 'loop_']

        # Locate start and end of atom site loop efficiently
        atom_site_loop_start = next((i for i in loop_ind_ids if '_atom_site.' in lines[i]), None)
        if atom_site_loop_start is None:
            raise ValueError('No atom site loop found')

        atom_site_loop_end = next((i + atom_site_loop_start for i, line in enumerate(lines[atom_site_loop_start:])
                                if line.startswith('ATOM')), None)
        
        if atom_site_loop_end is None:
            raise ValueError('No start of ATOM lines found')

        end_atom_lines = next((i + atom_site_loop_end - 1 for i, line in enumerate(lines[atom_site_loop_end:])
                            if not line.startswith('ATOM')), None)
        
        if end_atom_lines is None:
            raise ValueError('No end of ATOM lines found')

        # Headers processing
        headers = [x.strip() for x in lines[atom_site_loop_start:atom_site_loop_end] if x.strip()]
        header_map = {h: i for i, h in enumerate(headers)}  # Precompute index lookup

        chain_id_ind = header_map['_atom_site.label_asym_id']
        
        # Atom lines processing
        atom_lines = lines[atom_site_loop_end:end_atom_lines + 1]
        split_atom_lines = [l.split() for l in atom_lines]

        # Optimized chain finding
        chain_dict = {}
        prev_chain_id = None
        chain_lines = []

        for line in split_atom_lines:
            chain_id = line[chain_id_ind]
            if chain_id != prev_chain_id:
                if prev_chain_id is not None:
                    chain_dict[prev_chain_id] = chain_lines  # Store previous chain
                prev_chain_id = chain_id
                chain_lines = []
            chain_lines.append(line)
            
        if chain_lines:
            chain_dict[prev_chain_id] = chain_lines  # Store last chain

        return header_map, chain_dict

    def parse_pdb(self):
        '''
        parses a PDB file. similar to parse_cif
        '''

        # Read and split lines once
        lines = self.read()

        # find DBREF lines
        dbref_lines = [line for line in lines if line[:6] == 'DBREF ']
        
        # find lines with UNP as the DB.
        uniprot_lines = [line for line in dbref_lines if 'UNP' in line]

        # if dbref lines are found, extract uniprot id
        chain_to_uniprot = {}
        if uniprot_lines != []:
            for line in uniprot_lines:
                chain_id = line[12]
                uniprot_id = line[33:41].strip()
                chain_to_uniprot[chain_id] = uniprot_id
            

        # find ATOM lines
        atom_lines = [line for line in lines if line[:4] == 'ATOM']

        # Optimized chain finding
        chain_dict = {}
        prev_chain_id = None
        chain_lines = []

        for line in atom_lines:
            chain_id = line[21]
            if chain_id != prev_chain_id:
                if prev_chain_id is not None:
                    chain_dict[prev_chain_id] = chain_lines  # Store previous chain
                prev_chain_id = chain_id
                chain_lines = []
            chain_lines.append(line)

        if chain_lines:
            chain_dict[prev_chain_id] = chain_lines  # Store last chain

        return chain_dict, chain_to_uniprot        
    
    def build_complex_from_pdb(self):
        '''
        constructs a Complex object from a pdb file
        '''
        chain_dict, chain_to_uniprot = self.parse_pdb()

        for chain_id, lines in chain_dict.items():
            # see if we can get the uniprot ID in the PDB. If not, oh well, still is None.
            uniprot_id=None
            if chain_id in chain_to_uniprot:
                uniprot_id = chain_to_uniprot[chain_id]
            chain = Chain(chain_id=chain_id, uniprot_ID=uniprot_id)
            domain = Domain(domain_id=0, domain_type='unknown')

            # Map monomer numbers to their lines in one pass
            monomer_map = {}
            for token in lines:
                mon_num = int(token[22:26])
                monomer_map.setdefault(mon_num, []).append(token)

            # Process each monomer efficiently
            for mon_num, monomer_tokens in monomer_map.items():
                monomer = Monomer(monomer_number=mon_num, 
                                  monomer_name=monomer_tokens[0][17:20].strip(),
                                  occupancy=monomer_tokens[0][54:60].strip(),
                                  b_factor=monomer_tokens[0][60:66].strip(),
                                  charge=monomer_tokens[0][78:80].strip()
                )
                
                monomer.add_atoms([
                    Atom(atom_id=int(t[6:11]), 
                        atom_name=t[12:16].strip(),
                        x_coord=float(t[30:38]), 
                        y_coord=float(t[38:46]), 
                        z_coord=float(t[46:54]), 
                        element=t[76:78].strip())
                    for t in monomer_tokens
                ])
                domain.add_monomer(monomer)

            chain.add_domain(domain)
            self.complex.add_chain(chain)

        return self.complex

    def build_complex_form_cif(self):
        """
        Construct molecular complex from CIF file.
        """
        headers, chain_dict = self.parse_cif()

        # Extract indices once to avoid multiple `index()` calls
        atom_id_ind = headers['_atom_site.id']
        atom_name_ind = headers['_atom_site.auth_atom_id']
        chain_id_ind = headers['_atom_site.label_asym_id']
        monomer_name_ind = headers['_atom_site.auth_comp_id']
        monomer_number_ind = headers['_atom_site.auth_seq_id']
        x_ind = headers['_atom_site.Cartn_x']
        y_ind = headers['_atom_site.Cartn_y']
        z_ind = headers['_atom_site.Cartn_z']
        element_ind = headers['_atom_site.type_symbol']

        for chain_id, split_lines in chain_dict.items():
            chain = Chain(chain_id=chain_id)
            domain = Domain(domain_id=0, domain_type='unknown')

            # Map monomer numbers to their lines in one pass
            monomer_map = {}
            for tokens in split_lines:
                mon_num = int(tokens[monomer_number_ind])
                monomer_map.setdefault(mon_num, []).append(tokens)

            # Process each monomer efficiently
            for mon_num, monomer_tokens in monomer_map.items():
                monomer = Monomer(monomer_number=mon_num, 
                                  monomer_name=monomer_tokens[0][monomer_name_ind],
                                  occupancy=monomer_tokens[0][headers['_atom_site.occupancy']],
                                  b_factor=monomer_tokens[0][headers['_atom_site.B_iso_or_equiv']],
                                  charge=monomer_tokens[0][headers['_atom_site.pdbx_formal_charge']])
                monomer.add_atoms([
                    Atom(atom_id=t[atom_id_ind], atom_name=t[atom_name_ind],
                        x_coord=float(t[x_ind]), y_coord=float(t[y_ind]), 
                        z_coord=float(t[z_ind]), element=t[element_ind])
                    for t in monomer_tokens
                ])
                domain.add_monomer(monomer)

            chain.add_domain(domain)
            self.complex.add_chain(chain)

        return self.complex
    
    def return_complex(self):
        if self.filepath.endswith('.cif'):
            return self.build_complex_form_cif()
        elif self.filepath.endswith('.pdb'):
            return self.build_complex_from_pdb()
        else:
            raise ValueError('File format not supported. Only PDB and CIF files are supported.')
        

