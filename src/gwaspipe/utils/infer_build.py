"""
Genome build inference utility.

This module provides functionality to infer genome build version using HapMap3 SNPs.
"""

import pandas as pd
from gwaslab.g_Log import Log
from gwaslab.qc_fix_sumstats import finished
from gwaslab.qc_fix_sumstats import start_to

from gwaspipe.utils.change_status import vchange_status_from_version_3_6_16 as vchange_status
import importlib.resources


def infergenomebuild(
    sumstats, status="STATUS", chrom="CHR", pos="POS", ea="EA", nea="NEA", build="19", verbose=True, log=Log()
):
    """
    Infer genome build version using HapMap3 SNPs.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe
    status : str, default='STATUS'
        Column name for status
    chrom : str, default='CHR'
        Column name for chromosome
    pos : str, default='POS'
        Column name for position
    ea : str, default='EA'
        Column name for effect allele
    nea : str, default='NEA'
        Column name for non-effect allele
    build : str, default='19'
        Default build to use
    verbose : bool, default=True
        Whether to print verbose output
    log : Log
        Log object for recording operations

    Returns
    -------
    tuple
        (sumstats, inferred_build) where sumstats is the updated dataframe and
        inferred_build is the inferred genome build version
    """
    ##start function with col checking##########################################################
    _start_line = "infer genome build version using hapmap3 SNPs"
    _end_line = "inferring genome build version using hapmap3 SNPs"
    _start_cols = [chrom, pos]
    _start_function = ".infer_build()"
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
        return sumstats, "Unknown"
    ############################################################################################

    inferred_build = "Unknown"
    log.write("Start to infer genome build version using hapmap3 SNPs...", verbose=verbose)
    # Use importlib.resources to get the path from gwaslab package
    with importlib.resources.path("gwaslab.data.hapmap3_SNPs", "hapmap3_db150_hg19.snplist.gz") as data_path_19:
        data_path_19 = str(data_path_19)

    with importlib.resources.path("gwaslab.data.hapmap3_SNPs", "hapmap3_db151_hg38.snplist.gz") as data_path_38:
        data_path_38 = str(data_path_38)

    log.write(" -Loading Hapmap3 variants data...", verbose=verbose)
    hapmap3_ref_19 = pd.read_csv(
        data_path_19, sep=r"\s+", usecols=["#CHROM", "POS"], dtype={"#CHROM": "string", "POS": "string"}
    )
    hapmap3_ref_38 = pd.read_csv(
        data_path_38, sep=r"\s+", usecols=["#CHROM", "POS"], dtype={"#CHROM": "string", "POS": "string"}
    )

    log.write(" -CHR:POS will be used for matching...", verbose=verbose)
    raw_chrpos = sumstats[chrom].astype("string") + ":" + sumstats[pos].astype("string")

    hapmap3_ref_19["chr:pos"] = hapmap3_ref_19["#CHROM"] + ":" + hapmap3_ref_19["POS"]
    hapmap3_ref_38["chr:pos"] = hapmap3_ref_38["#CHROM"] + ":" + hapmap3_ref_38["POS"]

    match_count_for_19 = sum(raw_chrpos.isin(hapmap3_ref_19["chr:pos"].values))
    match_count_for_38 = sum(raw_chrpos.isin(hapmap3_ref_38["chr:pos"].values))

    log.write(" -Matching variants for hg19: num_hg19 = ", match_count_for_19, verbose=verbose)
    log.write(" -Matching variants for hg38: num_hg38 = ", match_count_for_38, verbose=verbose)

    if max(match_count_for_19, match_count_for_38) < 10000:
        log.warning("Please be cautious due to the limited number of variants.", verbose=verbose)

    if match_count_for_19 > match_count_for_38:
        log.write(" -Since num_hg19 >> num_hg38, assigning genome build hg19...", verbose=verbose)
        sumstats[status] = vchange_status(sumstats[status], 1, "9", "1")
        sumstats[status] = vchange_status(sumstats[status], 2, "9", "9")
        inferred_build = "19"
    elif match_count_for_19 < match_count_for_38:
        log.write(" -Since num_hg19 << num_hg38, assigning genome build hg38...", verbose=verbose)
        sumstats[status] = vchange_status(sumstats[status], 1, "9", "3")
        sumstats[status] = vchange_status(sumstats[status], 2, "9", "8")
        inferred_build = "38"
    else:
        log.write(" -Since num_hg19 = num_hg38, unable to infer...", verbose=verbose)

    finished(log, verbose, _end_line)
    return sumstats, inferred_build
