#!/usr/bin/env python

# CLI for making a PDB from a protein name

import os
import argparse

from dodo.build import pdb_from_name as build_pbd_from_name

def pdb_from_name():
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description='Generate a structure from a protein name that is based on an AF2 structure having redesigned IDRs.')
    # add args
    parser.add_argument('name', nargs='+', help='Name of the protein.')
    parser.add_argument('-o', '--out_path', help='Path to the output file.', required=True)
    parser.add_argument('-m', '--mode', default='predicted', help='Mode to use for generating the structure.')
    parser.add_argument('-n', '--num_models', default=1, help='Number of models to generate.', type=int)
    parser.add_argument('-l', '--linear_placement', action='store_true', default=False, help='Place the disordered regions linearly.')
    parser.add_argument('-c', '--no_CONECT_lines', action='store_false', default=True, help='Include CONECT lines in the output PDB.')
    parser.add_argument('-f', '--no_FD_atoms', action='store_false', default=True, help='Include the full disordered atoms in the output PDB.')
    parser.add_argument('-u', '--use_metapredict', action='store_true', default=False, help='Use Metapredict to predict the disordered regions.')
    parser.add_argument('-s', '--silent', action='store_true', default=False, help='Silence the output.')
    parser.add_argument('-apr', '--attempts_per_region', default=40, help='Number of attempts to make per region.', type=int)
    parser.add_argument('-apc', '--attempts_per_coord', default=2000, help='Number of attempts to make per coordinate.', type=int)
    
    # parser args
    args = parser.parse_args()

    # change name from list to string
    prot_name=''
    for i in args.name:
        prot_name+=f'{i} '
    protein_name=prot_name[:-1]

    if args.silent:
        verbose=False
    else:
        verbose=True

    # build pdb
    build_pbd_from_name(protein_name, out_path=args.out_path, mode=args.mode, 
        linear_placement=args.linear_placement, CONECT_lines=args.no_CONECT_lines, 
        include_FD_atoms=args.no_FD_atoms, use_metapredict=args.use_metapredict, 
        verbose=verbose, attempts_per_region=args.attempts_per_region, 
        attempts_per_coord=args.attempts_per_coord, num_models=args.num_models)

