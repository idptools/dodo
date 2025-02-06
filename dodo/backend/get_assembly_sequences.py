'''
For assemblies, we don't have full chain info. 
Even worse, we don't have uniprot IDs in .cif files because someone decided they want
to watch the world burn.
SO. For .cif files, the workaround is to query the structure ID at 
RCSB, get the scop ID for each unique chain, use that to get the full sequence
from the scop API, and then use that to get the full sequence for each chain.

For PDBs, we have luxury living because the full chain info is in the PDB.
'''
import os
import requests
import logging
import time
from typing import Dict, Optional
from requests.exceptions import RequestException

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SequenceRetrievalError(Exception):
    """Custom exception for sequence retrieval errors"""
    pass

def retry_request(url: str, max_retries: int = 3, delay: float = 1.0) -> Optional[dict]:
    """Helper function to retry failed requests"""
    for attempt in range(max_retries):
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.json()
        except RequestException as e:
            if attempt == max_retries - 1:
                logger.error(f"Failed to retrieve data from {url}: {str(e)}")
                return None
            time.sleep(delay * (attempt + 1))
    return None

def get_structure_uniprot_ids(structure_id: str) -> Dict[str, str]:
    """
    Retrieves UniProt IDs for chains in a given structure.

    Parameters
    ----------
    structure_id : str
        PDB/mmCIF structure identifier

    Returns
    -------
    Dict[str, str]
        Dictionary mapping chain IDs to UniProt IDs

    Raises
    ------
    SequenceRetrievalError
        If there is an error retrieving data from RCSB
    ValueError
        If structure_id is invalid
    """
    if not isinstance(structure_id, str) or not structure_id.strip():
        raise ValueError("Structure ID must be a non-empty string")
    
    structure_id = structure_id.strip().upper()
    chain_id_to_uniprot_id = {}

    try:
        # Get structure information
        url = f'https://data.rcsb.org/rest/v1/core/entry/{structure_id}'
        data = retry_request(url)
        
        if not data:
            raise SequenceRetrievalError(f"Failed to retrieve data for structure {structure_id}")
            
        try:
            ids = data['rcsb_entry_container_identifiers']['polymer_entity_ids']
        except KeyError as e:
            logger.error(f"Missing polymer entity IDs for structure {structure_id}")
            raise SequenceRetrievalError(f"Invalid structure data: {str(e)}")

        # Process each polymer entity
        for entity_id in ids:
            url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{structure_id}/{entity_id}'
            entity_data = retry_request(url)
            
            if not entity_data:
                logger.warning(f"Skipping entity {entity_id} - failed to retrieve data")
                continue
                
            try:
                chains_asym = entity_data['rcsb_polymer_entity_container_identifiers']['asym_ids']
                chains_auth_asym = entity_data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']
                uniprot_ids = entity_data['rcsb_polymer_entity_container_identifiers']['uniprot_ids']
                
                if not uniprot_ids:
                    logger.warning(f"No UniProt IDs found for entity {entity_id}")
                    continue
                    
                uniprot_id = uniprot_ids[0]  # Take first UniProt ID if multiple exist
                
                for chain in chains_auth_asym:    
                    chain_id_to_uniprot_id[chain] = uniprot_id
                # we can add additional things as long as we aren't overwriting previous values with non auth values. 
                for chain in chains_asym:
                    if chain not in chain_id_to_uniprot_id:
                        chain_id_to_uniprot_id[chain] = uniprot_id
                    
            except KeyError as e:
                logger.warning(f"Missing data for entity {entity_id}: {str(e)}")
                continue

        if not chain_id_to_uniprot_id:
            raise SequenceRetrievalError(f"No valid chain/UniProt mappings found for {structure_id}")

    except Exception as e:
        logger.error(f"Error processing structure {structure_id}: {str(e)}")
        raise SequenceRetrievalError(f"Failed to process structure {structure_id}: {str(e)}")

    return chain_id_to_uniprot_id

def get_sequences_form_uniprot_ids(chain_to_uniprot_ids_dict: Dict[str, str]) -> Dict[str, str]:
    """
    takes in chain to uniprot ID dict, returns chain to full length seq dict. 
    This function takes a dictionary mapping PDB chain IDs to UniProt IDs and returns
    a dictionary mapping chain IDs to their corresponding amino acid sequences by 
    querying the UniProt REST API.
    Parameters
    ----------
    chain_to_uniprot_ids_dict : Dict[str, str]
        Dictionary mapping PDB chain IDs (keys) to UniProt IDs (values)
    Returns
    -------
    Dict[str, str]
        Dictionary mapping PDB chain IDs (keys) to amino acid sequences (values)
    Raises
    ------
    ValueError
        If input is not a dictionary
    SequenceRetrievalError
        If there is an error retrieving sequences from UniProt
    Notes
    -----
    - Uses the UniProt REST API to fetch sequences
    - Skips chains where sequence retrieval fails
    - Logs warnings for chains where no sequence is found
    """
    
    if not isinstance(chain_to_uniprot_ids_dict, dict):
        raise ValueError("Input must be a dictionary")

    uniprot_seqs_to_fetch = list(set(chain_to_uniprot_ids_dict.values()))
    uniprot_id_to_seq = {}
    chain_seq_dict = {}

    try:
        for uniprot_id in uniprot_seqs_to_fetch:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
            data = retry_request(url)
            
            if not data:
                logger.error(f"Failed to retrieve sequence for UniProt ID: {uniprot_id}")
                continue

            try:
                sequence = data['sequence']['value']
                uniprot_id_to_seq[uniprot_id] = sequence
            except KeyError as e:
                logger.error(f"Missing sequence data for {uniprot_id}: {str(e)}")
                continue

        for chain_id, uniprot_id in chain_to_uniprot_ids_dict.items():
            if uniprot_id in uniprot_id_to_seq:
                chain_seq_dict[chain_id] = uniprot_id_to_seq[uniprot_id]
            else:
                logger.warning(f"No sequence found for chain {chain_id} (UniProt ID: {uniprot_id})")

    except Exception as e:
        logger.error(f"Unexpected error retrieving sequences: {str(e)}")
        raise SequenceRetrievalError("Failed to retrieve sequences from UniProt")

    return chain_seq_dict
