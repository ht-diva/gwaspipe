"""
Custom allele sorting functionality.
"""

from functools import cmp_to_key


def custom_alleles_sort(strings):
    """
    Sort alleles using custom comparison:
    - Single-length alleles come after multi-length
    - Alphabetical ordering within same length

    Parameters
    ----------
    strings : list
        List of alleles to sort

    Returns
    -------
    list
        Sorted alleles
    """

    def compare(a, b):
        # If both strings have the same length = 1
        if len(a) == len(b) == 1:
            return ord(a) - ord(b)
        # If both strings have length > 1
        elif len(a) == len(b) > 1:
            for char_a, char_b in zip(a, b):
                if char_a != char_b:
                    return ord(char_a) - ord(char_b)
            return len(a) - len(b)
        # If one string is longer than the other
        else:
            return len(b) - len(a)

    return sorted(strings, key=cmp_to_key(compare))
