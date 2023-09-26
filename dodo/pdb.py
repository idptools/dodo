# rewriting some PDB stuff.
import numpy as np

# in case array instead np.array is returned.. ugh
class array(list):
    def __init__(self, *args):
        list.__init__(self, *args)
        self.array = np.array(self)

def format_string(input_string, total_length):
    """
    Format a string to a specified total length by trimming or padding.

    Args:
    input_string (str): The input string to be formatted.
    total_length (int): The desired total length of the resulting string.

    Returns:
    str: The formatted string.
    """
    if len(input_string) > total_length:
        # Truncate the string if it's longer than the specified length
        formatted_string = input_string[:total_length]
    else:
        # Pad the string with spaces if it's shorter than the specified length
        formatted_string = input_string.ljust(total_length)

    return formatted_string


def write_single_chain_pdb():
    pass










