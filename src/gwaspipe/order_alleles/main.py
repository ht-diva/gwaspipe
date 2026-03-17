"""
Main order_alleles function.
"""

import pandas as pd
from gwaslab.g_Sumstats import _flip_allele_stats
from gwaslab.info.g_Log import Log

from .snpid import parallelbuildsnpid
from .vectorized import parallelorderalleles_status, vectorizedorderalleles_status


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
    3. Build snpid based on the new allele order -> chr:pos:ea:nea

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

    sumstats_data = _flip_allele_stats(sumstats_data, log=log)

    # Step 3: build snpid based on the new allele order -> chr:pos:ea:nea
    if format_snpid:
        sumstats_data = parallelbuildsnpid(
            sumstats_data, chrom=chrom, pos=pos, ea=ea, nea=nea, snpid=snpid, n_cores=n_cores, log=log, verbose=verbose
        )

    return sumstats_data
