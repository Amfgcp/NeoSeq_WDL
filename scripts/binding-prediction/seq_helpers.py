"""
Helper functions to deal with string sequences.
Developed to work with amino acid sequences.
"""

import logging

"""
Returns array of some peptides of provided 'size' created using a sliding window.
The last index where a sliding windows starts is `size` - 1.
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
Assumes varID format of: chrX_Y_W/Z. E.g., chr3_112563516_AGTATTCTGCCAAT/A
(frameshift variant) located in chromosome 3 position 112563516, where
in the reference there was a 'AGTATTCTGCCAAT' and in the tumor an 'A'.
"""
def check_if_frameshift(varID):
    split_varID = varID.split("/")
    ref = split_varID[0].split("_")[-1]
    alt = split_varID[-1]
    indel_len = abs(len(ref) - len(alt))
    if indel_len % 3 != 0:
        logging.debug('Found frameshift var %s', varID)
        return True
    return False

"""
Returns a tuple with the subsequence of specified 'size' centered (when possible)
around the mutation and its starting index on the input alt sequence, 'mut_seq'.
This index is useful for obtaining the Wild Type subsequence aligned at the start
of the alt subsequence.

NOTE: when 'size' is even and pos is in the right boundary of getting cut
("elif mut_len - pos - 1 < half_size" condition), it still gets to this
condition even though the size is just right to accommodate an even sized
sequence, i.e. it's technically a "center" case (last "else" condition). Can be
solved with even/odd check but this is unnecessary unless there is a need to
know, e.g., if the mutation is not centered in the return value of this function.
Note that an even sized sequence will never actually be centered around the
mutation as there is center.
"""
def take_sub_peptide(mut_seq, pos, size, is_frameshift):
    mut_len = len(mut_seq)
    if pos >= mut_len or pos < 0:
        incorrect_mut_pos = ("Mutation position is not valid, value is "
                             "negative or bigger than ALT sequence length. "
                             "The value of pos was: {} and length: {}").format(pos, mut_len - 1)
        raise Exception(incorrect_mut_pos)
    if mut_len <= size:
        logging.debug("Sequence length too small (%i) for provided size (%s), returning entire sequence", mut_len, size)
        return (mut_seq, 0)
    else:
        half_size = floor(size/2)
        if pos <= half_size:
            logging.debug("Sequence \"cut\" on left")
            if is_frameshift:
                return (mut_seq, 0)
            return (mut_seq[: size], 0)
        elif mut_len - pos - 1 < half_size:
            logging.debug("Sequence \"cut\" on right")
            return (mut_seq[mut_len - size:], mut_len - size)
        else:
            logging.debug("Mutation in range")
            if is_frameshift:
                return (mut_seq[pos - half_size:], pos - half_size)
            if size % 2 == 0:
                return (mut_seq[pos - half_size:pos + half_size], pos - half_size)
            return (mut_seq[pos - half_size:pos + half_size + 1], pos - half_size)

"""
Returns a sub-peptide of 'mut_seq' centered around 'mut_pos' and of size
2 * 'size' - 1.
"""
def calc_short_peptide(mut_seq, mut_pos, size):
    return mut_seq[max(0, mut_pos - size + 1):mut_pos + size]