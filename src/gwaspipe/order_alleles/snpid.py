"""
SNP ID building functionality.
"""

import pandas as pd
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_fix_sumstats import _df_split
from functools import partial
from multiprocessing import Pool


def build_snpids(sumstats, chrom="CHR", pos="POS", nea="NEA", ea="EA", snpid="SNPID"):
    """
    Build SNPID column in format CHR:POS:EA:NEA.

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

    Returns
    -------
    pd.Series
        Series containing built SNPID strings
    """
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
