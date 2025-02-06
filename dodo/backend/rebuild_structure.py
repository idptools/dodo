'''
Final backend module for the DODO project.
This module contains everything necessary for actually reconstructing
the structure of a protein. This module will ultimately 
be imported and used for all of the frontend functionality.
'''

import numpy as np
import os

from dodo.backend.generate_alpha_carbon_points import generate_ca_points
import idr_construction_utils as idr_utils
import position_fds
import utils
from dodo_readers import Reader
from dodo_structures import Atom, Monomer, Domain, Chain, Complex
from identify_domains import get_fds_loops_idrs

test=Reader('/Users/ryanemenecker/Desktop/lab_packages/dodo/dodo/data/6kn7.cif')
test_complex=test.return_complex()
print(test_complex.get_chain_to_uniprot_id_dict())