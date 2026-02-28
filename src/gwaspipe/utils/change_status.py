"""
Status change utility functions.

This module provides functions for changing STATUS codes in GWAS data.
"""

import pandas as pd

STATUS_CATEGORIES = [
    str(j + i)
    for j in [1300000, 1800000, 1900000, 3800000, 3900000, 9700000, 9800000, 9900000]
    for i in range(0, 100000)
]


def vchange_status_from_version_3_6_16(status, digit, before, after):
    """
    Change STATUS codes at specific digit positions.

    Parameters
    ----------
    status : pd.Series
        Series containing STATUS codes
    digit : int
        Digit position to change (1-indexed)
    before : str
        Characters to replace
    after : str
        Replacement characters

    Returns
    -------
    pd.Categorical
        Updated STATUS codes as categorical
    """
    dic = dict(zip(str(before), str(after)))
    left = status.str.slice(0, digit - 1)
    mid_orig = status.str.get(digit - 1)
    mid = mid_orig.map(dic).fillna(mid_orig)
    right = status.str.slice(digit, None)
    return pd.Categorical(left + mid + right, categories=STATUS_CATEGORIES)
