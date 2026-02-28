"""
Allele statistics flipping functionality.
"""

import pandas as pd
from gwaslab.g_Log import Log
from gwaslab.g_version import _get_version
from gwaslab.qc_fix_sumstats import (
    start_to,
    get_reverse_complementary_allele,
    flip_by_swap,
    flip_by_sign,
    flip_by_subtract,
    flip_by_inverse,
    finished,
)

from gwaspipe.utils.change_status import vchange_status_from_version_3_6_16 as vchange_status


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
    """
    Flip allele statistics based on STATUS codes.

    Parameters
    ----------
    sumstats : pd.DataFrame
        Summary statistics dataframe
    status : str, default='STATUS'
        Column name for status
    reverse_compl : bool, default=True
        Whether to reverse complement alleles
    flip_ref : bool, default=True
        Whether to flip reference alleles
    flip_ref_und : bool, default=True
        Whether to flip reference for indistinguishable indels
    flip_rev_strand : bool, default=True
        Whether to flip reverse strand palindromic variants
    verbose : bool, default=True
        Whether to print verbose output
    log : Log
        Log object for recording operations

    Returns
    -------
    pd.DataFrame
        Updated dataframe with flipped statistics
    """
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
