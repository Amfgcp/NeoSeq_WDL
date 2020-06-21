"""
Helper functions to deal with string sequences.
Developed to work with amino acid sequences.
"""

from math import floor
import logging

"""
Returns array of some peptides of provided 'size' created using a sliding window.
The last index where a sliding windows starts is 'size' - 1.
Note that for a given input 'short_pep' different sliding windows may result in
the same short peptide sequence (e.g, in repetitive sequences).
"""
def calc_mers(short_pep, size):
    results = []
    for i in range(size):
        slide_peptide = short_pep[i:i + size]
        if len(slide_peptide) < size:
            # depends on how len() works for arrays that are "indexed" outside of
            # bounds: only outputs the real length of the array
            break
        results.append(slide_peptide)
    return results

"""
Assumes 'var_id' format of: chrX_Y_W/Z. E.g., chr3_112563516_AGTATTCTGCCAAT/A
(frameshift variant) located in chromosome 3 position 112563516, where
in the reference there was a 'AGTATTCTGCCAAT' and in the tumor an 'A'.
"""
def check_if_frameshift(var_id):
    split_var_id = var_id.split("/")
    ref = split_var_id[0].split("_")[-1]
    mut = split_var_id[-1]
    indel_len = abs(len(ref) - len(mut))
    if indel_len % 3 != 0:
        logging.debug('Found frameshift var %s', var_id)
        return True
    return False

"""
Returns a tuple with the subsequence of specified 'size' centered (when possible)
around the mutation and wild type sequences, 'mut_pep' and 'wt_pep', respectively.

NOTE: when 'size' is even and pos is in the right boundary of getting cut
("elif mut_len - pos - 1 < half_size" condition), it still gets to this
condition even though the size is just right to accommodate an even sized
sequence, i.e. it's technically a "center" case (last "else" condition). Can be
solved with even/odd check but this is unnecessary unless there is a need to
know, e.g., if the mutation is not centered in the return value of this function.
Note that an even sized sequence will never actually be centered around the
mutation as there is center.
"""
def take_sub_peptide(mut_pep, wt_pep, pos, size, is_frameshift_or_empty_WT, mut_offset, wt_offset):
    mut_len = len(mut_pep)
    if pos >= mut_len or pos < 0:
        incorrect_mut_pos = ("Mutation position is not valid, value is "
            "negative or bigger than MUT sequence length. The value of pos "
            "was: {} and length: {}").format(pos, mut_len - 1)
        raise Exception(incorrect_mut_pos)
    if mut_len <= size:
        logging.debug("Sequence length too small (%i) for provided size (%s), \
                       returning entire sequence", mut_len, size)
        return (mut_pep, wt_pep[: size])
    else:
        half_size = floor(size/2)
        if pos <= half_size:
            logging.debug("Sequence \"cut\" on left")
            if is_frameshift_or_empty_WT:
                return (mut_pep, wt_pep[: size])
            return (mut_pep[: size], wt_pep[: size])
        elif mut_len - pos - 1 < half_size:
            logging.debug("Sequence \"cut\" on right")
            mut_pep_start = mut_len - size
            wt_sub_pep = wt_pep[mut_pep_start - mut_offset + wt_offset: \
                          mut_pep_start - mut_offset + wt_offset + size]
            return (mut_pep[mut_len - size:], wt_sub_pep)
        else:
            logging.debug("Mutation in range")
            mut_pep_start = pos - half_size
            wt_sub_pep = wt_pep[mut_pep_start - mut_offset + wt_offset: \
                          mut_pep_start - mut_offset + wt_offset + size]
            if is_frameshift_or_empty_WT:
                return (mut_pep[pos - half_size:], wt_sub_pep)
            if size % 2 == 0:
                return (mut_pep[pos - half_size:pos + half_size], wt_sub_pep)
            return (mut_pep[pos - half_size:pos + half_size + 1], wt_sub_pep)

"""
Returns a sub-peptide of 'mut_pep' centered around 'mut_pos' and of size
2 * 'size' - 1.
"""
def calc_short_peptide(mut_pep, mut_pos, size):
    return mut_pep[max(0, mut_pos - size + 1):mut_pos + size]