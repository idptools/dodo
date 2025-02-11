'''
various functions we want to use to
construct the IDR. Specifically, this functionality
will focus on taking the candidate coordinates from
the generate_alpha_carbon_points.py function and determine
which one is the best to use for the next alpha carbon
position
'''

import random

import numpy as np

from dodo.backend.dodo_structures import Atom, Monomer, Domain, Chain, Complex
from dodo.backend.generate_alpha_carbon_points import generate_ca_points
import dodo_math



