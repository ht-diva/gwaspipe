"""
External module for order_alleles functionality.
Other dependencies are imported from gwaslab modules.
"""

from functools import partial, cmp_to_key
from multiprocessing import Pool

import numpy as np
import pandas as pd
from gwaslab.bd_common_data import _maketrans
from gwaslab.g_Log import Log

from gwaslab.g_version import _get_version
from gwaslab.qc_fix_sumstats import (
    start_to,
    get_reverse_complementary_allele,
    flip_by_swap,
    flip_by_sign,
    flip_by_subtract,
    flip_by_inverse,
    _df_split,
    finished,
)

from gwaspipe.util_change_status import vchange_status_from_version_3_6_16 as vchange_status


# ===== Custom allele sorting =====


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


# ===== Vectorized allele ordering =====

# Create translation table for fast vectorized ordering
ORDER_MAPPING = {el: chr(i + 1) for i, el in enumerate(sorted(["A", "T", "C", "G"], reverse=False))}
assert all(value != chr(0) for value in ORDER_MAPPING.values()), "Mapping should not equal chr(0)"
TRANSLATE_TABLE_ORDER = _maketrans(ORDER_MAPPING)


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


# ===== Allele ordering functions =====


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


# ===== SNP ID building function =====
def build_snpids(sumstats, chrom="CHR", pos="POS", nea="NEA", ea="EA", snpid="SNPID"):
    return (
        sumstats[chrom].astype("string")
        + ":"
        + sumstats[pos].astype("string")
        + ":"
        + sumstats[ea].astype("string")
        + ":"
        + sumstats[nea].astype("string")
    )


def parallelbuildsnpid(
    sumstats, chrom="CHR", pos="POS", nea="NEA", ea="EA", snpid="SNPID", n_cores=1, verbose=True, log=Log()
):
    """
    Build SNPID column in format CHR:POS:EA:NEA in parallel.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe
    chrom : str, default='CHR'
        Column name for chromosome
    pos : str, default='POS'
        Column name for position
    nea : str, default='NEA'
        Column name for non-effect allele
    ea : str, default='EA'
        Column name for effect allele
    snpid : str, default='SNPID'
        Column name for SNP ID
    n_cores : int, default=1
        Number of cores for parallel processing
    verbose : bool, default=True
        Whether to print verbose output
    log : Log
        Log object for recording operations

    Returns
    -------
    pd.DataFrame
        Updated dataframe with built SNPID column
    """

    if snpid in sumstats.columns:
        log.write(f"Start to build SNPID column on {n_cores} cores...", verbose=verbose)

        if n_cores > 1:
            df_split = _df_split(sumstats[[chrom, pos, ea, nea]], n_cores)
            pool = Pool(n_cores)
            map_func = partial(build_snpids, chrom=chrom, pos=pos, nea=nea, ea=ea, snpid=snpid)
            snpids = pd.concat(pool.map(map_func, df_split))
            pool.close()
            pool.join()
        else:
            snpids = build_snpids(sumstats, chrom=chrom, pos=pos, nea=nea, ea=ea, snpid=snpid)

        sumstats[snpid] = snpids
        log.write("Finished building SNPID column.", verbose=verbose)
    else:
        log.warning(f"'{snpid}' column is not found in the DataFrame. Skipping the build of SNPID.", verbose=verbose)

    return sumstats


def _flip_allele_statistics(
    sumstats,
    status="STATUS",
    reverse_compl=True,
    flip_ref=True,
    flip_ref_und=True,
    flip_rev_strand=True,
    verbose=True,
    log=Log(),
):
    ##start function with col checking##########################################################
    _start_line = "adjust statistics based on STATUS code"
    _end_line = "adjusting statistics based on STATUS code"
    _start_cols = []
    _start_function = ".flip_allele_stats()"
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

    if_stats_flipped = False
    ###################get reverse complementary####################
    if reverse_compl:
        # pattern = r"\w\w\w\w\w[45]\w"
        # matched_index = status_match(sumstats[status],6,[4,5]) #
        matched_index = sumstats[status].str[5].str.match(r"4|5")
        if sum(matched_index) > 0:
            log.write(
                "Start to convert alleles to reverse complement for SNPs with status xxxxx[45]x...{}".format(
                    _get_version()
                ),
                verbose=verbose,
            )
            log.write(" -Flipping " + str(sum(matched_index)) + " variants...", verbose=verbose)
            if ("NEA" in sumstats.columns) and ("EA" in sumstats.columns):
                log.write(" -Converting to reverse complement : EA and NEA...", verbose=verbose)
                reverse_complement_nea = sumstats.loc[matched_index, "NEA"].apply(
                    lambda x: get_reverse_complementary_allele(x)
                )
                reverse_complement_ea = sumstats.loc[matched_index, "EA"].apply(
                    lambda x: get_reverse_complementary_allele(x)
                )
                categories = (
                    set(sumstats["EA"])
                    | set(sumstats["NEA"])
                    | set(reverse_complement_ea)
                    | set(reverse_complement_nea)
                )
                sumstats["EA"] = pd.Categorical(sumstats["EA"], categories=categories)
                sumstats["NEA"] = pd.Categorical(sumstats["NEA"], categories=categories)
                sumstats.loc[matched_index, ["NEA"]] = reverse_complement_nea
                sumstats.loc[matched_index, ["EA"]] = reverse_complement_ea
                sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 6, "4", "2")
                log.write(" -Changed the status for flipped variants : xxxxx4x -> xxxxx2x", verbose=verbose)
            if_stats_flipped = True
    ###################flip ref####################
    if flip_ref:
        # pattern = r"\w\w\w\w\w[35]\w"
        # matched_index = status_match(sumstats[status],6,[3,5]) #sumstats[status].str.match(pattern)
        matched_index = sumstats[status].str[5].str.match(r"3|5")
        if sum(matched_index) > 0:
            log.write(
                "Start to flip allele-specific stats for SNPs with status xxxxx[35]x: ALT->EA , REF->NEA ...{}".format(
                    _get_version()
                ),
                verbose=verbose,
            )
            log.write(" -Flipping " + str(sum(matched_index)) + " variants...", verbose=verbose)

            flip_by_swap(sumstats, matched_index, log, verbose)
            flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
            flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
            flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)

            # change status
            log.write(" -Changed the status for flipped variants : xxxxx[35]x -> xxxxx[12]x", verbose=verbose)
            sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 6, "35", "12")
            if_stats_flipped = True

    ###################flip ref for undistingushable indels####################
    if flip_ref_und:
        # pattern = r"\w\w\w\w[123][67]6"
        # matched_index = status_match(sumstats[status],6,[1,2,3])|status_match(sumstats[status],6,[6,7])|status_match(sumstats[status],7,6) #sumstats[status].str.match(pattern)
        matched_index = sumstats[status].str[4:].str.match(r"[123][67]6")
        if sum(matched_index) > 0:
            log.write(
                "Start to flip allele-specific stats for standardized indels with status xxxx[123][67][6]: ALT->EA , REF->NEA...{}".format(
                    _get_version()
                ),
                verbose=verbose,
            )
            log.write(" -Flipping " + str(sum(matched_index)) + " variants...", verbose=verbose)

            flip_by_swap(sumstats, matched_index, log, verbose)
            flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
            flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
            flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)

            # change status
            log.write(" -Changed the status for flipped variants xxxx[123][67]6 -> xxxx[123][67]4", verbose=verbose)
            sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 7, "6", "4")
            if_stats_flipped = True
            # flip ref
    ###################flip statistics for reverse strand panlindromic variants####################
    if flip_rev_strand:
        # pattern = r"\w\w\w\w\w[012]5"
        # matched_index = status_match(sumstats[status],6,[0,1,2]) | status_match(sumstats[status],7,[5])#sumstats[status].str.match(pattern)
        matched_index = sumstats[status].str[5:].str.match(r"05|15|25")
        if sum(matched_index) > 0:
            log.write(
                "Start to flip allele-specific stats for palindromic SNPs with status xxxxx[12]5: (-)strand <=> (+)strand...{}".format(
                    _get_version()
                ),
                verbose=verbose,
            )
            log.write(" -Flipping " + str(sum(matched_index)) + " variants...", verbose=verbose)

            flip_by_sign(sumstats, matched_index, log, verbose, cols=None)
            flip_by_subtract(sumstats, matched_index, log, verbose, cols=None, factor=1)
            flip_by_inverse(sumstats, matched_index, log, verbose, cols=None, factor=1)

            # change status
            log.write(" -Changed the status for flipped variants:  xxxxx[012]5: ->  xxxxx[012]2", verbose=verbose)
            sumstats.loc[matched_index, status] = vchange_status(sumstats.loc[matched_index, status], 7, "5", "2")
            if_stats_flipped = True

    if not if_stats_flipped:
        log.write(" -No statistics have been changed.")

    finished(log, verbose, _end_line)
    return sumstats


# ===== Main order_alleles function =====


def order_alleles(
    sumstats_data,
    log=None,
    ea="EA",
    nea="NEA",
    status="STATUS",
    chrom="CHR",
    pos="POS",
    snpid="SNPID",
    format_snpid=True,
    n_cores=1,
    mode="v",
    flipallelestats_args=None,
    verbose=True,
):
    """
    Order alleles based on custom ordering and update statistics accordingly.

    This function reorders effect and non-effect alleles based on a defined
    order, flips allele-specific statistics accordingly, and optionally
    rebuilds SNP IDs with the new allele order.

    **Three-step process:**
    1. Set status to appropriate value if ea and nea should be flipped based on custom ordering
    2. Fix stats to match the new allele order (ea and nea are swapped where needed and stats are fixed)
    3. Build snpid based on the new allele order -> chr:pos:allele1:allele2

    Parameters
    ----------
    sumstats_data : pd.DataFrame
        Summary statistics dataframe containing allele and status columns
    log : Log, optional
        Log object for recording operations. If None, creates new Log()
    ea : str, default='EA'
        Column name for effect allele
    nea : str, default='NEA'
        Column name for non-effect allele
    status : str, default='STATUS'
        Column name for status information
    chrom : str, default='CHR'
        Column name for chromosome
    pos : str, default='POS'
        Column name for position
    snpid : str, default='SNPID'
        Column name for SNP ID
    format_snpid : bool, default=True
        Whether to rebuild SNP ID based on new allele order
    n_cores : int, default=1
        Number of cores for parallel processing
    mode : str, default='v'
        Processing mode: 'v' for vectorized or 'p' for parallel
    flipallelestats_args : dict, optional
        Additional arguments to pass to flipallelestats
    verbose : bool, default=True
        Whether to print verbose output

    Returns
    -------
    pd.DataFrame
        Updated dataframe with reordered alleles and adjusted statistics

    Examples
    --------
    >>> # Vectorized mode (faster, recommended)
    >>> sumstats_ordered = order_alleles(sumstats_df, mode='v')

    >>> # Parallel mode (multi-threaded)
    >>> sumstats_ordered = order_alleles(sumstats_df, mode='p', n_cores=4)

    >>> # Custom parameters
    >>> sumstats_ordered = order_alleles(
    ...     sumstats_df,
    ...     ea='A1',
    ...     nea='A2',
    ...     format_snpid=True,
    ...     n_cores=2,
    ...     mode='p'
    ... )
    """
    # Handle empty DataFrame case
    if sumstats_data.empty:
        if verbose:
            log.log("Empty DataFrame provided to order_alleles")
        return sumstats_data

    if log is None:
        log = Log()

    if flipallelestats_args is None:
        flipallelestats_args = {}

    # Step 1: set status to appropriate value if ea and nea should be flipped
    # based on custom ordering
    if mode == "v":
        sumstats_data = vectorizedorderalleles_status(
            sumstats_data, ea=ea, nea=nea, status=status, log=log, verbose=verbose
        )
    elif mode == "p":
        sumstats_data = parallelorderalleles_status(
            sumstats_data, ea=ea, nea=nea, status=status, n_cores=n_cores, log=log, verbose=verbose
        )
    else:
        if verbose:
            log.log("Invalid mode. Using vectorized mode.")
        sumstats_data = vectorizedorderalleles_status(
            sumstats_data, ea=ea, nea=nea, status=status, log=log, verbose=verbose
        )

    # Step 2: fix stats to match the new allele order
    # (ea and nea are swapped where needed and stats are fixed accordingly)
    categories = set()
    if sumstats_data[ea].dtype.name == "category":
        categories = categories | set(sumstats_data[ea].cat.categories.tolist())
    if sumstats_data[nea].dtype.name == "category":
        categories = categories | set(sumstats_data[nea].cat.categories.tolist())

    sumstats_data[ea] = pd.Categorical(sumstats_data[ea], categories=categories)
    sumstats_data[nea] = pd.Categorical(sumstats_data[nea], categories=categories)

    # Use default flipallelestats args optimized for speed
    base_flipallelestats_args = dict(reverse_compl=False, flip_ref=True, flip_ref_und=False, flip_rev_strand=False)
    base_flipallelestats_args.update(flipallelestats_args)
    sumstats_data = _flip_allele_statistics(sumstats_data, log=log, **base_flipallelestats_args)

    # Step 3: build snpid based on the new allele order -> chr:pos:ea:nea
    if format_snpid:
        sumstats_data = parallelbuildsnpid(
            sumstats_data, chrom=chrom, pos=pos, ea=ea, nea=nea, snpid=snpid, n_cores=n_cores, log=log, verbose=verbose
        )

    return sumstats_data
