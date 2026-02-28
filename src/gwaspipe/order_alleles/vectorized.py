"""
Vectorized allele ordering functionality.
"""

import numpy as np
import pandas as pd
from gwaslab.g_Log import Log

from gwaspipe.utils.change_status import vchange_status_from_version_3_6_16 as vchange_status
from gwaslab.qc_fix_sumstats import start_to, finished

from .constants import TRANSLATE_TABLE_ORDER
from .sorting import custom_alleles_sort


def _orderalleles_status_vec(sumstats, nea="NEA", ea="EA", status="STATUS", verbose=True, log=Log()):
    """
    Vectorized status ordering for fast processing of large datasets.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe
    nea : str, default='NEA'
        Column name for non-effect allele
    ea : str, default='EA'
        Column name for effect allele
    status : str, default='STATUS'
        Column name for status
    verbose : bool, default=True
        Whether to print verbose output
    log : Log
        Log object for recording operations

    Returns
    -------
    pd.DataFrame
        Updated dataframe with reordered status codes
    """
    if sumstats.empty:
        return sumstats

    # Translate the strings to integer numpy arrays in a very fast way
    _ea = sumstats[ea]
    max_len_ea = _ea.str.len().max()
    _ea = _ea.str.translate(TRANSLATE_TABLE_ORDER).to_numpy().astype(f"<U{max_len_ea}")
    _ea = _ea.view("<u4").reshape(-1, max_len_ea).astype(np.uint8)

    _nea = sumstats[nea]
    max_len_nea = _nea.str.len().max()
    _nea = _nea.str.translate(TRANSLATE_TABLE_ORDER).to_numpy().astype(f"<U{max_len_nea}")
    _nea = _nea.view("<u4").reshape(-1, max_len_nea).astype(np.uint8)

    ea_non_zero = np.sum(_ea != 0, axis=1)
    nea_non_zero = np.sum(_nea != 0, axis=1)

    # First condition: swap if NEA is longer than EA
    should_swap = nea_non_zero > ea_non_zero

    # When NEA and EA have the same length, check if the first different value is smaller
    equal_length = nea_non_zero == ea_non_zero
    ea_equal = _ea[equal_length]
    nea_equal = _nea[equal_length]
    max_shape = max(ea_equal.shape[1], nea_equal.shape[1]) if ea_equal.shape[0] > 0 else 0

    if max_shape > 0:
        if ea_equal.shape[1] > nea_equal.shape[1]:
            nea_equal = np.pad(nea_equal, ((0, 0), (0, max_shape - _nea.shape[1])))
        elif ea_equal.shape[1] < nea_equal.shape[1]:
            ea_equal = np.pad(ea_equal, ((0, 0), (0, max_shape - _ea.shape[1])))

        first_difference = np.argmax(nea_equal != ea_equal, axis=1)
        ea_different_val = np.take_along_axis(ea_equal, first_difference[:, None], axis=1)
        nea_different_val = np.take_along_axis(nea_equal, first_difference[:, None], axis=1)
        cond_ordering = (nea_different_val < ea_different_val).flatten()
        should_swap[equal_length] = cond_ordering

    log.write(
        f"  -For Flipped match ({sum(should_swap)} matches): convert STATUS xxxxx[0123456789]x to xxxxx3x...",
        verbose=verbose,
    )
    sumstats.loc[should_swap, status] = vchange_status(sumstats.loc[should_swap, status], 6, "0123456789", "3" * 10)

    return sumstats


def vectorizedorderalleles_status(sumstats, nea="NEA", ea="EA", status="STATUS", verbose=True, log=Log()):
    """
    Vectorized allele ordering wrapper with threshold for large/small alleles.

    Uses fast vectorized operations for processing speed.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe
    nea : str, default='NEA'
        Column name for non-effect allele
    ea : str, default='EA'
        Column name for effect allele
    status : str, default='STATUS'
        Column name for status
    verbose : bool, default=True
        Whether to print verbose output
    log : Log
        Log object for recording operations

    Returns
    -------
    pd.DataFrame
        Updated dataframe with reordered status codes
    """

    ##start function with col checking##########################################################
    _start_line = "change status based on custom allele order"
    _end_line = "changing status based on custom allele order"
    _start_cols = [nea, ea, status]
    _start_function = ".order_alleles()"
    _must_args = {}

    is_enough_info = start_to(
        sumstats=sumstats,
        log=log,
        verbose=verbose,
        start_line=_start_line,
        end_line=_end_line,
        start_cols=_start_cols,
        start_function=_start_function,
        **_must_args,
    )
    if not is_enough_info:
        return sumstats
    ############################################################################################

    max_len = 4  # chosen threshold for vectorized processing
    condition = (sumstats[nea].str.len() <= max_len) & (sumstats[ea].str.len() <= max_len)

    log.write(f" -Changing status for records with ( len(NEA) <= {max_len} and len(EA) <= {max_len} )", verbose=verbose)
    sumstats_cond = sumstats[condition]
    out = _orderalleles_status_vec(sumstats_cond, nea=nea, ea=ea, status=status, verbose=verbose, log=log)
    sumstats.loc[condition, status] = out[status]

    log.write(f" -Changing status for records with ( len(NEA) > {max_len} or len(EA) > {max_len} )", verbose=verbose)
    sumstats_not_cond = sumstats[~condition]
    out = _orderalleles_status_vec(sumstats_not_cond, nea=nea, ea=ea, status=status, verbose=verbose, log=log)
    sumstats.loc[~condition, status] = out[status]

    finished(log, verbose, _end_line)

    return sumstats


def orderalleles_status(sumstats, nea="NEA", ea="EA", status="STATUS", verbose=True, log=Log()):
    """
    Non-vectorized allele ordering using custom sorting.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe
    nea : str, default='NEA'
        Column name for non-effect allele
    ea : str, default='EA'
        Column name for effect allele
    status : str, default='STATUS'
        Column name for status
    verbose : bool, default=True
        Whether to print verbose output
    log : Log
        Log object for recording operations

    Returns
    -------
    pd.DataFrame
        Updated dataframe with reordered status codes
    """

    def status_ordering(ea, nea, status):
        to_sort = [ea, nea]
        allele1, allele2 = custom_alleles_sort(to_sort)
        if allele1 != ea:
            status = status[:5] + "3" + status[6:]
        return status

    out = sumstats[[ea, nea, status]].apply(lambda x: status_ordering(x.iloc[0], x.iloc[1], x.iloc[2]), axis=1)

    if sumstats[status].dtype.name == "category":
        sumstats[status] = pd.Categorical(out.values, categories=sumstats[status].cat.categories)
    else:
        sumstats[status] = out.values

    return sumstats


def parallelorderalleles_status(sumstats, nea="NEA", ea="EA", status="STATUS", n_cores=1, verbose=True, log=Log()):
    """
    Parallel allele ordering for distributed processing across multiple cores.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe
    nea : str, default='NEA'
        Column name for non-effect allele
    ea : str, default='EA'
        Column name for effect allele
    status : str, default='STATUS'
        Column name for status
    n_cores : int, default=1
        Number of cores for parallel processing
    verbose : bool, default=True
        Whether to print verbose output
    log : Log
        Log object for recording operations

    Returns
    -------
    pd.DataFrame
        Updated dataframe with reordered status codes
    """
    from functools import partial
    from multiprocessing import Pool

    ##start function with col checking##########################################################
    _start_line = "change status based on custom allele order"
    _end_line = "changing status based on custom allele order"
    _start_cols = [nea, ea, status]
    _start_function = ".order_alleles()"
    _must_args = {}

    is_enough_info = start_to(
        sumstats=sumstats,
        log=log,
        verbose=verbose,
        start_line=_start_line,
        end_line=_end_line,
        start_cols=_start_cols,
        start_function=_start_function,
        n_cores=n_cores,
        **_must_args,
    )
    if not is_enough_info:
        return sumstats
    ############################################################################################

    if n_cores > 1:
        from gwaslab.qc_fix_sumstats import _df_split

        df_split = _df_split(sumstats, n_cores)
        pool = Pool(n_cores)
        map_func = partial(orderalleles_status, ea=ea, nea=nea, status=status, verbose=verbose, log=log)
        sumstats = pd.concat(pool.map(map_func, df_split))
        pool.close()
        pool.join()
    else:
        sumstats = orderalleles_status(sumstats, ea=ea, nea=nea, status=status, verbose=verbose, log=log)

    finished(log, verbose, _end_line)
    return sumstats
