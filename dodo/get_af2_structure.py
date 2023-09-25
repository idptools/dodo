# code for getting af2 structures from the internets
import requests
import io
from getSequence import getseq
import os
from dodo.dodo_exceptions import dodoAF2Exception
from dodo.dodo_tools import af2_lines_to_structure_dict


def get_af2_pdb(protein_name, outpath='', silent=True, save=False):
    '''
    function to get the af2 pdb from uniprot using the protein name.
    Uses a pypi package last updated in 2015, so if functionality breaks
    check wget package

    protein_name : string
        the protein name as a string. Can be multiple words
        like 'human p53' or just p53. Can also be tune uniprot ID. 
        Will take the *top* hit on Uniprot.

    outpath : string
        Optional. Can specify the output location to save the pdb. Otherwise,
        just saves to cwd as the af2 name.

    silent : bool
        whether to print the protein name recieved from uniprot.

    save : bool 
        whether to save a pdb. Otherwise returns the lines.

    '''
    # get the uniprot ID
    try:
        prot_id_and_seq = getseq(protein_name)
    except:
        raise Exception('Unable to find protein!')
    prot_id=prot_id_and_seq[0].split('|')[1].upper()

    #generate the URL
    url=f'https://alphafold.ebi.ac.uk/files/AF-{prot_id}-F1-model_v4.pdb'
    if outpath=='':
        outpath=f'{os.getcwd()}/data/AF-{prot_id}-F1-model_v4.pdb'
    # make sure we are getting the absolute path to where 
    # this saves so we can nuke it later.
    abspath=os.path.abspath(outpath)
    # print output
    if silent==False:
        print(f'Downloading AF2 structure for {prot_id_and_seq[0]}')

    # download the af2 structure
    response=requests.get(url)
    in_memory_file = io.StringIO(response.text)
    lines=[]
    for line in in_memory_file:
        lines.append(line.strip())

    af2_seq = af2_lines_to_structure_dict(lines)['sequence']
    if af2_seq != prot_id_and_seq[1]:
        raise dodoAF2Exception('AF2 sequence does not match sequence retrieved from Uniprot!')
    return lines




