
# From BioPython: https://github.com/biopython/biopython/blob/c5a6b1374267d769b19c1022b4b45472316e78b4/Bio/Seq.py#L36
def _maketrans(complement_mapping):
    """Make a python string translation table (PRIVATE).

    Arguments:
     - complement_mapping - a dictionary.

    Returns a translation table (a bytes object of length 256) for use with
    the python string's translate method.

    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    keys = "".join(complement_mapping.keys()).encode("ASCII")
    values = "".join(complement_mapping.values()).encode("ASCII")
    return bytes.maketrans(keys + keys.lower(), values + values.lower())