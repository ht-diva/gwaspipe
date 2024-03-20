import numpy as np
import pandas as pd
from typing import Dict

import gwaslab.hm_harmonize_sumstats
import gwaslab.bd_common_data as bd_common_data
import gwaslab.qc_fix_sumstats as qc_fix_sumstats
import gwaslab.g_Log as g_Log

from Bio import SeqIO

from .utils import _maketrans


### CONSTANTS AND MAPPINGS ###

PADDING_VALUE = 100

# chr(0) should not be used in the mapping dict, because it's a reserved value.
# Instead of starting from chr(1), we start from chr(2) because this could be useful in the future
# to compute the complementary allele with a simple XOR operation (e.g. 2 ^ 1 = 3, 3 ^ 1 = 2, 4 ^ 1 = 5, 5 ^ 1 = 4, ...)
MAPPING = {
    "A": chr(2),
    "T": chr(3),
    "C": chr(4),
    "G": chr(5),
    "N": chr(6),
}
assert all(value != chr(0) for value in MAPPING.values()), "Mapping in the dictionary should not be equal to chr(0). This is a reserved value"

_COMPLEMENTARY_MAPPING = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "N": "N",
}
COMPLEMENTARY_MAPPING = {k: MAPPING[v] for k,v in _COMPLEMENTARY_MAPPING.items()}

TRANSLATE_TABLE = _maketrans(MAPPING)
TRANSLATE_TABLE_COMPL = _maketrans(COMPLEMENTARY_MAPPING)


### FUNCTIONS ###

# From https://github.com/Cloufield/gwaslab/blob/f6b4c4e58a26e5d67d6587141cde27acf9ce2a11/src/gwaslab/hm_harmonize_sumstats.py#L293
# It's basically the original function but with the addition of the fast_check_status function call (the try-except block)
def fast_checkref(sumstats,ref_path,chrom="CHR",pos="POS",ea="EA",nea="NEA",status="STATUS",chr_dict=bd_common_data.get_chr_to_number(),remove=False,verbose=True,log=g_Log.Log()):
    #log.write("++++++++++ CUSTOM FAST CHECKREF ++++++++++",verbose=verbose)  
    ##start function with col checking##########################################################
    _start_line = "check if NEA is aligned with reference sequence"
    _end_line = "checking if NEA is aligned with reference sequence"
    _start_cols = [chrom,pos,ea,nea,status]
    _start_function = ".check_ref()"
    _must_args ={}

    is_enough_info = qc_fix_sumstats.start_to(sumstats=sumstats,
                              log=log,
                              verbose=verbose,
                              start_line=_start_line,
                              end_line=_end_line,
                              start_cols=_start_cols,
                              start_function=_start_function,
                              **_must_args)
    if is_enough_info == False: return sumstats
    ############################################################################################
    log.write(" -Reference genome FASTA file: "+ ref_path,verbose=verbose)  
    log.write(" -Loading fasta records", verbose=verbose)
    chromlist = bd_common_data.get_chr_list(add_number=True)
    records = SeqIO.parse(ref_path, "fasta")

    all_records_dict = {}
    chroms_in_sumstats = sumstats[chrom].unique() # load records from Fasta file only for the chromosomes present in the sumstats
    for record in records:
        #record = next(records)
        if record is not None:
            record_chr = str(record.id).strip("chrCHR").upper()
            if record_chr in chr_dict.keys():
                i = chr_dict[record_chr]
            else:
                i = record_chr
            if (i in chromlist) and (i in chroms_in_sumstats):
                all_records_dict.update({i: record})

    # Try to apply the fast implemention, if it fails, fall back to the original implementation 
    if len(all_records_dict) > 0:
        try:
            log.write(" -Checking records with custom fast_check_status()", verbose=verbose)
            all_records_dict = dict(sorted(all_records_dict.items())) # sort by key in case the fasta records are not already ordered by chromosome
            to_check_ref = (sumstats[chrom].isin(list(all_records_dict.keys()))) & (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna())
            sumstats_to_check = sumstats.loc[to_check_ref,[chrom,pos,ea,nea,status]]
            sumstats.loc[to_check_ref,status] = fast_check_status(sumstats_to_check, all_records_dict, log=log, verbose=verbose)
            log.write(" -Finished checking records", verbose=verbose)
        except:
            log.write(f"Error in checking records with the fast_check_status function. Falling back to original gwaslab check_status implementation.",verbose=verbose)
            log.write(" -Checking records: ", end="",verbose=verbose) 
            for record_chr, record in all_records_dict.items():
                log.write(record_chr," ", end="",show_time=False,verbose=verbose) 
                to_check_ref = (sumstats[chrom]==i) & (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna())
                sumstats.loc[to_check_ref,status] = sumstats.loc[to_check_ref,[pos,ea,nea,status]].apply(lambda x:gwaslab.hm_harmonize_sumstats.check_status(x,record),axis=1)
            log.write("\n",end="",show_time=False,verbose=verbose) 
 
    sumstats[status] = sumstats[status].astype("string")

    # These checks are quite expensive and are used only for logging purposes.
    # We could remove/disable them if we want to speed up the function
    available_to_check =sum( (~sumstats[pos].isna()) & (~sumstats[nea].isna()) & (~sumstats[ea].isna()))
    status_0=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[0]\w", case=False, flags=0, na=False))
    status_3=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[3]\w", case=False, flags=0, na=False))
    status_4=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[4]\w", case=False, flags=0, na=False))
    status_5=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[5]\w", case=False, flags=0, na=False))
    status_6=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[6]\w", case=False, flags=0, na=False))
    #status_7=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[7]\w", case=False, flags=0, na=False))
    status_8=sum(sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w", case=False, flags=0, na=False))
    
    log.write(" -Variants allele on given reference sequence : ",status_0,verbose=verbose)
    log.write(" -Variants flipped : ",status_3,verbose=verbose)
    raw_matching_rate = (status_3+status_0)/available_to_check
    flip_rate = status_3/available_to_check
    log.write("  -Raw Matching rate : ","{:.2f}%".format(raw_matching_rate*100),verbose=verbose)
    if raw_matching_rate <0.8:
        log.warning("Matching rate is low, please check if the right reference genome is used.")
    if flip_rate > 0.85 :
        log.write("  -Flipping variants rate > 0.85, it is likely that the EA is aligned with REF in the original dataset.",verbose=verbose)
    
    log.write(" -Variants inferred reverse_complement : ",status_4,verbose=verbose)
    log.write(" -Variants inferred reverse_complement_flipped : ",status_5,verbose=verbose)
    log.write(" -Both allele on genome + unable to distinguish : ",status_6,verbose=verbose)
    #log.write(" -Reverse_complementary + both allele on genome + unable to distinguish: ",status_7)
    log.write(" -Variants not on given reference sequence : ",status_8,verbose=verbose)
    
    if remove is True:
        sumstats = sumstats.loc[~sumstats["STATUS"].str.match("\w\w\w\w\w[8]\w"),:]
        log.write(" -Variants not on given reference sequence were removed.",verbose=verbose)
    
    qc_fix_sumstats.finished(log, verbose, _end_line)
    return sumstats


def build_fasta_records(fasta_records_dict: Dict, pos_as_dict=True, log=g_Log.Log(), verbose=True):
    log.write("   -Building numpy fasta records from dict", verbose=verbose)

    # Let's do some magic to convert the fasta record to a numpy array of integers in a very fast way.
    # fasta_record.seq._data is a byte-string, so we can use the bytes.maketrans to apply a translation.
    # Here we map the bytes to the unicode character representing the desired integer as defined in the mapping dict
    # (i.e. b'A' -> '\x02', b'T' -> '\x03', b'C' -> '\x04', b'G' -> '\x05', b'N' -> '\x06')
    # Then, using np.array(... dtype=<U..) we convert the string to a numpy array of unicode characters.
    # Then, we do a magic with view('<u4') to convert the unicode characters to 4-byte integers, so we obtain the actual integer representation of the characters
    # Lastly, we cast the array to np.uint8 to convert the 4-byte integers to 1-byte integers to save memory
    # Full example:
    # fasta_record.seq._data = b'ACTGN' -> b'\x02\x04\x03\x05\x06' -> np.array(['\x02\x04\x03\x05\x06'], dtype='<U5') -> np.array([2, 4, 3, 5, 6], dtype=uint32) -> np.array([2, 4, 3, 5, 6], dtype=uint8)
    all_r = []
    for r in fasta_records_dict.values():
        r = r.seq._data.translate(TRANSLATE_TABLE)
        r = np.array([r], dtype=f'<U{len(r)}').view('<u4').astype(np.uint8)
        all_r.append(r)

    # We've just created a list of numpy arrays, so we can concatenate them to obtain a single numpy array
    # Then we keep track of the starting position of each record in the concatenated array. This will be useful later
    # to index the record array depending on the position of the variant and the chromosome
    records_len = np.array([len(r) for r in all_r])
    starting_positions = np.cumsum(records_len) - records_len
    if pos_as_dict:
        starting_positions = {k: v for k, v in zip(fasta_records_dict.keys(), starting_positions)}
    record = np.concatenate(all_r)
    del all_r # free memory

    return record, starting_positions


def fast_check_status(sumstats: pd.DataFrame, fasta_records_dict: Dict, log=g_Log.Log(), verbose=True):

    chrom,pos,ea,nea,status = sumstats.columns

    # First, convert the fasta records to a single numpy array of integers
    record, starting_positions_dict = build_fasta_records(fasta_records_dict, pos_as_dict=True, log=log, verbose=verbose)

    # In _fast_check_status(), several 2D numpy arrays are created and they are padded to have shape[1] == max_len_nea or max_len_ea
    # Since most of the NEA and EA strings are short, we perform the check first on the records having short NEA and EA strings,
    # and then we perform the check on the records having long NEA and EA strings. In this way we can speed up the process (since the 
    # arrays are smaller) and save memory.
    max_len = 4 # this is a chosen value, we could compute it using some stats about the length and count of NEA and EA strings
    condition = (sumstats[nea].str.len() <= max_len) * (sumstats[ea].str.len() <= max_len)

    log.write(f"   -Checking records for ( len(NEA) <= {max_len} and len(EA) <= {max_len} )", verbose=verbose)
    sumstats_cond = sumstats[condition]
    starting_pos_cond = np.array([starting_positions_dict[k] for k in sumstats_cond[chrom].unique()])
    sumstats.loc[condition, status] = _fast_check_status(sumstats_cond, record=record, starting_positions=starting_pos_cond)

    log.write(f"   -Checking records for ( len(NEA) > {max_len} or len(EA) > {max_len} )", verbose=verbose)
    sumstats_not_cond = sumstats[~condition]
    starting_not_pos_cond = np.array([starting_positions_dict[k] for k in sumstats_not_cond[chrom].unique()])
    sumstats.loc[~condition, status] = _fast_check_status(sumstats_not_cond, record=record, starting_positions=starting_not_pos_cond)

    return sumstats[status].values


def _fast_check_status(x: pd.DataFrame, record: np.array, starting_positions: np.array):
    # status 
    #0 /  ----->  match
    #1 /  ----->  Flipped Fixed
    #2 /  ----->  Reverse_complementary Fixed
    #3 /  ----->  flipped
    #4 /  ----->  reverse_complementary 
    #5 / ------>  reverse_complementary + flipped
    #6 /  ----->  both allele on genome + unable to distinguish
    #7 /  ----> reverse_complementary + both allele on genome + unable to distinguish
    #8 / -----> not on ref genome
    #9 / ------> unchecked
    if x.empty:
        return np.array([])
    
    # x is expected to be a DataFrame with these columns in that order: ['CHR', 'POS', 'EA', 'NEA', 'STATUS']
    # In this way, we don't need to specify the columns names
    _chrom = x.iloc[:, 0]
    _pos = x.iloc[:, 1]
    _ea = x.iloc[:, 2]
    _nea = x.iloc[:, 3]
    _status = x.iloc[:, 4]

    # position of the status (i.e. x['STATUS']) that will be modified
    status_flip_idx = 5 

    pos = _pos.values.astype(np.int64) # convert to int64 because they could be of type 'object'

    # Rebase the chromosome numbers to 0-based indexing
    # e.g. ['1', '2', '4', '2'] -> [0, 1, 2, 1]
    # This is needed because record is a single 1D array containing all the records for all the selected chromosomes,
    # so for instance if record contains the records for chr1, chr2, chr4 ([...chr1...chr2...chr4...]), we need to 
    # rebase the chromosome numbers to 0-based indexing to index the correct record portion when we do starting_positions[chrom]
    # Note that in x there are only the rows for the same chromosomes for which we have the records in record
    # (i.e. we don't have rows for chr3 if we don't have the record for chr3). This filtering is done in the caller function
    _chrom = _chrom.values
    unique_values, _ = np.unique(_chrom, return_inverse=True) # Get the sorted unique values and their indices
    chrom = np.searchsorted(unique_values, _chrom) # Replace each value in '_chrom' with its corresponding index in the sorted unique values

    max_len_nea = _nea.str.len().max()
    max_len_ea = _ea.str.len().max()


    # Let's apply the same magic used for the fasta records (check build_fasta_records() for details) to convert the NEA and EA to
    # a numpy array of integers in a very fast way.
    # In that case we start from a pd.Series to we can apply some built-in methods.
    # Also, when doing nea.view('<u4'), each row will be automatically right-padded with zeros to reach the max_len_nea.
    # For this reason, we then replace the zeros with out padding value
    # (and that's why the mapping dict can't have chr(0) as a value, otherwise we would have zeros for both padding and a character)
    # Reshaping is needed because .view('<u4') will create a flattened array    
    nea = _nea.str.translate(TRANSLATE_TABLE).to_numpy().astype(f'<U{max_len_nea}')
    nea = nea.view('<u4').reshape(-1, max_len_nea).astype(np.uint8)
    nea[nea == 0] = PADDING_VALUE # padding value

    # Create a mask holding True at the position of non-padding values
    mask_nea = nea != PADDING_VALUE

    # Create the reverse complement of NEA
    # In this case, we manually left-pad the translated string with the padding value, since the padding done by view('<u4') would be right-padded
    # and that will make hard the reverse operation (because we would have e.g. [2, 2, 4, 100, ..., 100] which will be hard to convert into [4, 2, 2, 100, ..., 100])
    rev_nea = _nea.str.translate(TRANSLATE_TABLE_COMPL).str.pad(max_len_nea, 'left', chr(PADDING_VALUE)).to_numpy().astype(f'<U{max_len_nea}')
    rev_nea = rev_nea.view('<u4').reshape(-1, max_len_nea).astype(np.uint8)
    rev_nea = rev_nea[:, ::-1]


    # Let's do everything again for EA
    ea = _ea.str.translate(TRANSLATE_TABLE).to_numpy().astype(f'<U{max_len_ea}')
    ea = ea.view('<u4').reshape(-1, max_len_ea).astype(np.uint8)
    ea[ea == 0] = PADDING_VALUE # padding value

    mask_ea = ea != PADDING_VALUE

    rev_ea = _ea.str.translate(TRANSLATE_TABLE_COMPL).str.pad(max_len_ea, 'left', chr(PADDING_VALUE)).to_numpy().astype(f'<U{max_len_ea}')
    rev_ea = rev_ea.view('<u4').reshape(-1, max_len_ea).astype(np.uint8)
    rev_ea = rev_ea[:, ::-1]


    # Convert the status (which are integers represented as strings) to a numpy array of integers.
    # Again, use the same concept as before to do this in a very fast way.
    # e.g. ["9999999", "9939999", "9929999"] -> [[9, 9, 9, 9, 9, 9, 9], [9, 9, 3, 9, 9, 9, 9], [9, 9, 2, 9, 9, 9, 9]]
    assert _status.str.len().value_counts().nunique() == 1 # all the status strings should have the same length, let's be sure of that.
    status_len = len(_status.iloc[0])
    mapping_status = {str(v): chr(v) for v in range(10)}
    table_stats = _maketrans(mapping_status)
    status = _status.str.translate(table_stats).to_numpy().astype(f'<U{status_len}')
    status = status.view('<u4').reshape(-1, status_len).astype(np.uint8)


    # Expand the position to a 2D array and subtract 1 to convert to 0-based indexing
    # e.g. [2, 21, 46] -> [[1], [20], [45]]
    pos = np.expand_dims(pos, axis=-1) - 1

    # Create a modified indices array specifying the starting position of each chromosome in the concatenated record array
    modified_indices = starting_positions[chrom]
    modified_indices = modified_indices[:, np.newaxis] # Add a new axis to modified_indices to align with the dimensions of pos

    # Create the range of indices: [0, ..., max_len_nea-1]
    indices_range = np.arange(max_len_nea)

    # Add the range of indices to the starting indices
    # e.g. pos = [[1], [20], [45]], indices_range = [0, 1, 2], indices = [[1, 2, 3], [20, 21, 22], [45, 46, 47]]
    indices = pos + indices_range

    # Modify indices to select the correct absolute position in the concatenated record array
    indices = indices + modified_indices

    # Let's pad the fasta records array because if there is a (pos, chrom) for which (pos+starting_position[chrom]+max_len_nea > len(record) we get out of bounds error.
    # This basically happens if there is a pos for the last chromosome for which pos+max_len_nea > len(record for that chrom).
    # This is very unlikely to happen but we should handle this case.
    record = np.pad(record, (0, max_len_nea), constant_values=PADDING_VALUE)
    
    # Index the record array using the computed indices.
    # Since we use np.take, indices must all have the same length, and this is why we added the padding to NEA
    # and we create the indices using max_len_nea (long story short, we can't obtain a scattered/ragged array)
    output_nea = np.take(record, indices)

    # Check if the NEA is equal to the reference sequence at the given position
    # In a non-matrix way, this is equivalent (for one single element) to:
    # nea == record[pos-1: pos+len(nea)-1]
    # where for example:
    #  a) nea = "AC", record = "ACTG", pos = 1 -> True
    #  b) nea = "T", record = "ACTG", pos = 3 -> True
    #  c) nea = "AG", record = "ACTG", pos = 1 -> False
    # Since we want to do everything in a vectorized way, we will compare the padded NEA with the output 
    # and then we use the mask to focus only on the non-padded elements
    # Pseudo example (X represents the padding value):
    #  nea = ['AC', 'T'], record = 'ACTGAAG', pos = [1, 3]
    #  -> nea = ['AC', 'TX'], indices = [[1, 2], [3, 4]], mask = [[True, True], [True, False]], output_nea = [['A', 'C'], ['T', 'G']]
    #  -> nea == output_nea: [[True, True], [True, False]], mask: [[True, True], [True, False]]
    #  -> nea == output_nea + ~mask: [[True, True], [True, True]]
    #  -> np.all(nea == output_nea + ~mask, 1): [True, True]
    nea_eq_ref = np.all((nea == output_nea) + ~mask_nea, 1)
    rev_nea_eq_ref = np.all((rev_nea == output_nea) + ~mask_nea, 1)

    # Let's do everything again for EA
    indices_range = np.arange(max_len_ea)
    indices = pos + indices_range
    indices = indices + modified_indices
    output_ea = np.take(record, indices)

    ea_eq_ref = np.all((ea == output_ea) + ~mask_ea, 1)
    rev_ea_eq_ref = np.all((rev_ea == output_ea) + ~mask_ea, 1)

    masks_max_len = max(mask_nea.shape[1], mask_ea.shape[1])

    len_nea_eq_len_ea = np.all(
        np.pad(mask_nea, ((0,0),(0, masks_max_len-mask_nea.shape[1])), constant_values=False) == 
        np.pad(mask_ea, ((0,0),(0, masks_max_len-mask_ea.shape[1])), constant_values=False)
        , axis=1) # pad masks with False to reach same shape
    len_rev_nea_eq_rev_len_ea = len_nea_eq_len_ea

    # The following conditions replicates the if-else statements of the original check_status function:
    # https://github.com/Cloufield/gwaslab/blob/f6b4c4e58a26e5d67d6587141cde27acf9ce2a11/src/gwaslab/hm_harmonize_sumstats.py#L238

    # nea == ref && ea == ref && len(nea) != len(ea)
    status[nea_eq_ref * ea_eq_ref * ~len_nea_eq_len_ea, status_flip_idx] = 6

    # nea == ref && ea != ref
    status[nea_eq_ref * ~ea_eq_ref, status_flip_idx] = 0

    # nea != ref && ea == ref
    status[~nea_eq_ref * ea_eq_ref, status_flip_idx] = 3

    # nea != ref && ea != ref && rev_nea == ref && rev_ea == ref && len(rev_nea) != len(rev_ea)
    status[~nea_eq_ref * ~ea_eq_ref * rev_nea_eq_ref * rev_ea_eq_ref * ~len_rev_nea_eq_rev_len_ea, status_flip_idx] = 8

    # nea != ref && ea != ref && rev_nea == ref && rev_ea != ref
    status[~nea_eq_ref * ~ea_eq_ref * rev_nea_eq_ref * ~rev_ea_eq_ref, status_flip_idx] = 4

    # nea != ref && ea != ref && rev_nea != ref && rev_ea == ref
    status[~nea_eq_ref * ~ea_eq_ref * ~rev_nea_eq_ref * rev_ea_eq_ref, status_flip_idx] = 5

    # nea != ref && ea != ref && rev_nea != ref && rev_ea != ref
    status[~nea_eq_ref * ~ea_eq_ref * ~rev_nea_eq_ref * ~rev_ea_eq_ref, status_flip_idx] = 8

    # Convert back the (now modified) 2D status array to a numpy array of strings in a very fast way.
    # Since 'status' is a 2D array of integers ranging from 0 to 9, we can build the integer representation
    # of each row using the efficent operation below (e.g. [1, 2, 3, 4, 5] -> [12345]).
    # Then we convert this integer to a string using the f'<U{status.shape[1]}' dtype (e.g. 12345 -> '12345')
    # The "naive" way would be:
    #   status_str = [''.join(map(str, l)) for l in status]
    #   status_arr = np.array(status_str)
    status_flat = np.sum(status * 10**np.arange(status.shape[1]-1, -1, -1), axis=1)
    status_arr = status_flat.astype(f'<U{status.shape[1]}')

    return status_arr