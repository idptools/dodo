'''
holds code that defines the following structures:
    - Atoms
    - Monomers
    - Domains
    - Chains
    - Complexes

    Atoms are the basic unit. This is equivalent to a bead.
    Mers are a collection of atoms. This could make up a residue or a nucleic acid or something else.
    Domains are a collection of mers. This could be a protein domain.
    Chains are a collection of domains. This could be a protein chain.
    Complexes are a collection of chains. This could be a protein complex.
'''

import numpy as np
from scipy.spatial import cKDTree
from utils import amino_acid_3_to_1

class Atom:
    __slots__ = ['atom_id', 'atom_name', 'x_coord', 'y_coord', 'z_coord', 'element']
    
    def __init__(self, atom_id=None, atom_name=None, 
                 x_coord=None, y_coord=None, z_coord=None, 
                 element=None):
        self.atom_id = atom_id
        self.atom_name = atom_name
        self.x_coord = float(x_coord)
        self.y_coord = float(y_coord)
        self.z_coord = float(z_coord)
        self.element = element
    
    def __str__(self):
        return f'Atom: {self.atom_id} {self.atom_name} {self.x_coord} {self.y_coord} {self.z_coord} {self.element}'
    
    def __repr__(self):
        return f'Atom: {self.atom_id} {self.atom_name} {self.x_coord} {self.y_coord} {self.z_coord} {self.element}'
    
    def __eq__(self, other):
        return self.atom_id == other.atom_id and self.atom_name == other.atom_name and  self.x_coord == other.x_coord and self.y_coord == other.y_coord and self.z_coord == other.z_coord and self.element == other.element 
    
    def get_atom_id(self):
        return self.atom_id

    def get_atom_name(self):
        return self.atom_name

    def coordinates(self):
        return (self.x_coord, self.y_coord, self.z_coord)
    
    def get_element(self):
        return self.element
    
    def set_atom_id(self, atom_id):
        self.atom_id = atom_id

    def set_coordinates(self, x_coord, y_coord, z_coord):
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.z_coord = z_coord
    
class Monomer:
    __slots__ = ['monomer_number', 'monomer_name', 'occupancy', 
                 'b_factor', 'charge', 'atoms', '_coord_cache']
    
    def __init__(self, 
                 monomer_number=None, 
                 monomer_name=None,
                 occupancy=None,
                 b_factor=None,
                 charge=None):
        self.monomer_number = monomer_number
        self.monomer_name = monomer_name
        self.occupancy = occupancy
        self.b_factor = b_factor
        self.charge = charge
        self.atoms = []
        self._coord_cache = None
    
    def __str__(self):
        return f'monomer: {self.monomer_number} {self.monomer_name} {self.atoms}'
    
    def __repr__(self):
        return f'monomer: {self.monomer_number} {self.monomer_name} {self.atoms}'
    
    def __eq__(self, other):
        return self.monomer_number == other.monomer_number and self.monomer_name == other.monomer_name and self.atoms == other
    
    def get_monomer_number(self):
        return self.monomer_number
    
    def get_monomer_name(self):
        return self.monomer_name
    
    def get_atoms(self):
        return self.atoms
    
    def get_atoms_dict(self):
        return {atom.get_atom_id(): {'coords':atom.coordinates(), 'name':atom.get_atom_name()} for atom in self.atoms}
    
    def get_chain_id(self):
        return self.chain_id
    
    def add_atom(self, atom):
        self.atoms.append(atom)
        self.invalidate_cache()

    def add_atoms(self, atoms):
        self.atoms.extend(atoms)
        self.invalidate_cache()
    
    def remove_atom(self, atom):
        self.atoms.remove(atom)
        self.invalidate_cache()
    
    def get_atom(self, atom_id):
        for atom in self.atoms:
            if atom.atom_id == atom_id:
                return atom
        return None

    def get_coordinates_array(self):
        """Cache and return numpy array of coordinates"""
        if self._coord_cache is None:
            self._coord_cache = np.array([[atom.x_coord, atom.y_coord, atom.z_coord] 
                                        for atom in self.atoms])
        return self._coord_cache
    
    def invalidate_cache(self):
        """Invalidate coordinate cache when atoms change"""
        self._coord_cache = None

class Domain:
    __slots__ = ['domain_id', 'domain_type', 'monomers', 'loop_indices', '_coord_cache', '_spatial_index']
    
    def __init__(self, 
                 domain_id=None, 
                 domain_type=None, 
                 monomers=[]):
        self.domain_id = domain_id
        self.domain_type = domain_type
        self.monomers = monomers
        self.loop_indices = []
        self._coord_cache = None
        self._spatial_index = None
    
    def __str__(self):
        return f'Domain: {self.domain_id} {self.domain_type} ({len(self.monomers)} monomers)'
    
    def __repr__(self):
        return self.__str__()
    
    def add_monomer(self, monomer):
        self.monomers.append(monomer)
        self.invalidate_cache()
    
    def get_monomers(self):
        return self.monomers
    
    def get_coordinates_array(self):
        """Cache and return numpy array of coordinates"""
        if self._coord_cache is None:
            coords = []
            for mon in self.monomers:
                coords.extend(mon.get_coordinates_array())
            self._coord_cache = np.array(coords)
        return self._coord_cache
    
    def build_spatial_index(self):
        """Build KD-tree for spatial queries"""
        coords = self.get_coordinates_array()
        self._spatial_index = cKDTree(coords)
    
    def rotate(self, rotation_matrix):
        """Rotate domain by rotation matrix"""
        coords = self.get_coordinates_array()
        rotated_coords = np.dot(coords, rotation_matrix.T)
        self._update_coordinates(rotated_coords)
    
    def translate(self, translation_vector):
        """Translate domain by vector"""
        coords = self.get_coordinates_array()
        translated_coords = coords + translation_vector
        self._update_coordinates(translated_coords)
    
    def _update_coordinates(self, new_coords):
        """Update atomic coordinates in domain"""
        idx = 0
        for monomer in self.monomers:
            for atom in monomer.atoms:
                atom.x_coord, atom.y_coord, atom.z_coord = new_coords[idx]
                idx += 1
        self.invalidate_cache()
    
    def invalidate_cache(self):
        """Invalidate coordinate cache"""
        self._coord_cache = None
        self._spatial_index = None

    def assign_domain_type(self, domain_type):
        self.domain_type = domain_type

    def add_loop_indices(self, indices):
        self.loop_indices.append(indices)

    def get_loop_indices(self):
        return self.loop_indices

    def translate_domain_to_origin(self, point=(0, 0, 0)):
        '''
        translates the domain such that
        the center of mass of the domain is 
        located at a specific point. 
        Default is the origin.
        '''
        center_of_mass = self.identify_center_of_mass()
        translation_vector = np.array(point) - center_of_mass
        self.translate(translation_vector)


class Chain:
    __slots__ = ['chain_id', 'domains', '_coord_cache', '_spatial_index', 'uniprot_ID', 'full_sequence']
    
    def __init__(self, chain_id=None, uniprot_ID=None):
        self.chain_id = chain_id
        self.domains = []
        self._coord_cache = None
        self._spatial_index = None
        self.uniprot_ID = uniprot_ID
        self.full_sequence = None
    
    def __str__(self):
        return f'Chain: {self.chain_id} {self.domains}'
    
    def __repr__(self):
        return f'Chain: {self.chain_id} {self.domains}'
    
    def __eq__(self, other):
        return self.chain_id == other.chain_id and self.domains == other.domains
    
    def get_chain_id(self):
        return self.chain_id
    
    def get_domains(self):
        return self.domains
    
    def add_domain(self, domain):
        self.domains.append(domain)
        self.invalidate_cache()
    
    def remove_domain(self, domain):
        self.domains.remove(domain)
        self.invalidate_cache()
    
    def get_domain(self, domain_id):
        for domain in self.domains:
            if domain.domain_id == domain_id:
                return domain
        return None

    def build_spatial_index(self):
        """Build KD-tree for spatial queries"""
        coords = self.get_coordinates_array()
        self._spatial_index = cKDTree(coords)
        
    def get_coordinates_array(self):
        """Cache and return numpy array of all coordinates in chain"""
        if self._coord_cache is None:
            coords = []
            for domain in self.domains:
                coords.extend(domain.get_coordinates_array())
            self._coord_cache = np.array(coords)
        return self._coord_cache
    

    def invalidate_cache(self):
        """Invalidate coordinate cache when domains change"""
        self._coord_cache = None
        self._spatial_index = None


    def get_all_atom_coordinates(self):
        """Cached version of coordinate retrieval"""
        if self._coord_cache is None:
            coordinates = []
            atom_info = []
            for domain in self.domains:
                for monomer in domain.monomers:
                    for atom in monomer.atoms:
                        coordinates.append([atom.x_coord, atom.y_coord, atom.z_coord])
                        atom_info.append((atom, monomer, domain))

            self._coord_cache = (np.array(coordinates), atom_info)
        
        return self._coord_cache
    
    def get_coords_by_monomer(self):
        """Return coordinates by monomer"""
        coords = {}
        for domain in self.domains:
            for monomer in domain.monomers:
                coords[monomer.get_monomer_number()]=monomer.get_coordinates_array()
        return coords
    
    def rotate(self, rotation_matrix):
        """Batch rotate all coordinates"""
        coords, atom_info = self.get_all_atom_coordinates()
        rotated_coords = np.dot(coords, rotation_matrix.T)
        self._update_coordinates(rotated_coords)
        
    def translate(self, translation_vector):
        """Batch translate all coordinates"""
        coords, atom_info = self.get_all_atom_coordinates()
        translated_coords = coords + translation_vector
        self._update_coordinates(translated_coords)

    def identify_center_of_mass(self):
        """Identify the center of mass of the chain"""
        coords = self.get_coordinates_array()
        return np.mean(coords, axis=0)
    
    def get_uniprot_id(self):
        return self.uniprot_ID
    
    def set_full_sequence(self, sequence):
        self.full_sequence = sequence
    
    def get_full_sequence(self):
        return self.full_sequence
    
    def get_current_sequence(self):
        sequence = ''
        for domain in self.domains:
            for monomer in domain.monomers:
                sequence+=amino_acid_3_to_1(monomer.get_monomer_name())
        return sequence

    def find_missing_residues(self):
        """Find missing residues in the chain"""
        sequence = self.get_current_sequence()
        missing_residues = []
        for i, aa in enumerate(self.full_sequence):
            if aa != sequence[i]:
                missing_residues.append((i, aa))
        return missing_residues




    

class Complex:
    __slots__ = ['complex_id', 'chains', '_coord_cache', '_spatial_index', 'filepath']
    
    def __init__(self, complex_id=None, filepath=None):
        self.complex_id = complex_id
        self.chains = []
        self._coord_cache = None
        self._spatial_index = None
        self.filepath = filepath
    
    def __str__(self):
        return f'Complex: {self.complex_id} {self.chains}'
    
    def __repr__(self):
        return f'Complex: {self.complex_id} {self.chains}'
    
    def __eq__(self, other):
        return self.complex_id == other.complex_id and self.chains == other.chains
    
    def get_complex_id(self):
        return self.complex_id
    
    def get_chains(self):
        return self.chains
    
    def add_chain(self, chain):
        self.chains.append(chain)
        self.invalidate_cache()
    
    def remove_chain(self, chain):
        self.chains.remove(chain)
        self.invalidate_cache()
    
    def get_chain(self, chain_id):
        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain
        return None

    def build_spatial_index(self):
        """Build KD-tree for spatial queries"""
        coords = self.get_coordinates_array()
        self._spatial_index = cKDTree(coords)
    
    def get_coordinates_array(self):
        """Cache and return numpy array of all coordinates in complex"""
        if self._coord_cache is None:
            coords = []
            for chain in self.chains:
                coords.extend(chain.get_coordinates_array())
            self._coord_cache = np.array(coords)
        return self._coord_cache
    
    def get_atoms_within_radius(self, point, radius):
        """Efficiently find atoms within radius of point"""
        if self._spatial_index is None:
            self.build_spatial_index()
        indices = self._spatial_index.query_ball_point(point, radius)
        return indices

    def invalidate_cache(self):
        """Invalidate coordinate cache when chains change"""
        self._coord_cache = None
        self._spatial_index = None

    def rotate(self, rotation_matrix):
        """Batch rotate all coordinates"""
        coords = self.get_coordinates_array()
        rotated_coords = np.dot(coords, rotation_matrix.T)
        self._update_coordinates(rotated_coords)
        
    def translate(self, translation_vector):
        """Batch translate all coordinates"""
        coords = self.get_coordinates_array()
        translated_coords = coords + translation_vector
        self._update_coordinates(translated_coords)

    def _update_coordinates(self, new_coords):
        """Update atomic coordinates and invalidate caches"""
        idx = 0
        for chain in self.chains:
            for domain in chain.domains:
                for monomer in domain.monomers:
                    for atom in monomer.atoms:
                        atom.x_coord, atom.y_coord, atom.z_coord = new_coords[idx]
                        idx += 1
        self._coord_cache = None
        self._spatial_index = None
    
    def get_all_atom_coordinates(self):
        """Cached version of coordinate retrieval"""
        if self._coord_cache is None:
            coordinates = []
            atom_info = []
            
            for chain in self.chains:
                for domain in chain.domains:
                    for monomer in domain.monomers:
                        for atom in monomer.atoms:
                            coordinates.append([atom.x_coord, atom.y_coord, atom.z_coord])
                            atom_info.append((atom, monomer, domain, chain))
            
            self._coord_cache = (np.array(coordinates), atom_info)
        
        return self._coord_cache
    
    def calculate_all_atom_distances(self):
        """Calculate distances between all atoms using vectorized operations."""
        coords, atom_info = self.get_all_atom_coordinates()
        
        # Calculate pairwise distances using broadcasting
        # This creates a matrix of shape (n_atoms, n_atoms)
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        distances = np.sqrt(np.sum(diff * diff, axis=-1))
        
        return distances, atom_info
    
    def check_atom_clashes(self, distance_threshold=2.0):
        """Optimized clash detection using spatial partitioning"""
        if self._spatial_index is None:
            coords, atom_info = self.get_all_atom_coordinates()
            self._spatial_index = cKDTree(coords)
        
        # Find all pairs of points within distance_threshold
        pairs = self._spatial_index.query_pairs(distance_threshold)
        
        # Convert to clash information
        clashes = []
        coords, atom_info = self._coord_cache
        for i, j in pairs:
            dist = np.sqrt(np.sum((coords[i] - coords[j])**2))
            clashes.append((atom_info[i], atom_info[j], dist))
        
        return clashes
    
    def get_chain_id_seq_dict(self):
        chain_id_seq_dict = {}
        for chain in self.chains:
            chain_id_seq_dict[chain.chain_id] = ''
            for domain in chain.domains:
                for monomer in domain.monomers:
                    chain_id_seq_dict[chain.chain_id]+=amino_acid_3_to_1(monomer.get_monomer_name())
        return chain_id_seq_dict
    
    def get_chain_to_uniprot_id_dict(self):
        chain_to_uniprot_id_dict = {}
        for chain in self.chains:
            chain_to_uniprot_id_dict[chain.chain_id] = chain.uniprot_ID
        return chain_to_uniprot_id_dict

    
    def return_filepath(self):
        return self.filepath
    
 