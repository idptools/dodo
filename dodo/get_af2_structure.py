# code for getting af2 structures from the internets
import requests
import io
from getSequence import getseq
import os
from dodo.dodo_exceptions import dodoAF2Exception
from dodo.pdb_tools import array, PDBParser

def get_af2_pdb_lines(protein_name, verbose=True, return_lines=False):
    '''
    function to get the af2 pdb from uniprot using the protein name.
    Uses a pypi package last updated in 2015, so if functionality breaks
    check wget package

    protein_name : string
        the protein name as a string. Can be multiple words
        like 'human p53' or just p53. Can also be tune uniprot ID. 
        Will take the *top* hit on Uniprot.

    verbose : bool
        whether to print the protein name recieved from uniprot.

    return_lines : bool
        whether to return the lines of the pdb file as a list

    '''
    # get the uniprot ID
    try:
        prot_id_and_seq = getseq(protein_name)
    except:
        raise Exception('Unable to find protein!')
    prot_id=prot_id_and_seq[0].split('|')[1].upper()

    #generate the URL
    url=f'https://alphafold.ebi.ac.uk/files/AF-{prot_id}-F1-model_v4.pdb'
    # print output
    if verbose==True:
        print(f'Downloading AF2 structure for {prot_id_and_seq[0]}')

    # download the af2 structure
    response=requests.get(url)
    in_memory_file = io.StringIO(response.text)
    lines=[]
    for line in in_memory_file:
        lines.append(line.strip())
    PDBParserObj=PDBParser(lines)
    sequence=PDBParserObj.sequence
    if sequence != prot_id_and_seq[1]:
        raise dodoAF2Exception('AF2 sequence does not match sequence retrieved from Uniprot!')
    if return_lines==True:
        return lines
    return PDBParserObj


