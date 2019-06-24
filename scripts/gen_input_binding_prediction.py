"""
This script is meant to generate the inputs for the binding prediction phase (e.g netMHC)

Usage:
  gen_input_binding_prediction.py -f NIC4.25peptide.txt -p 8 -d /path/to/db

Options:
  -h --help        Show this screen.
  --version        Show version.
  -p=8,9,10,11,12  Peptide length to be generated
  -f=<file>        Three column file containing varID, WTpeptide and ALTpeptide. ALT peptide has been currated for phases variants (using ex. ISOVAR)
  -d=<path/to/db> Path to database
"""

import fileinput
import sys
import re
import os

from docopt import docopt

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
from Bio.Blast import NCBIXML

from math import log, floor

'''
Returns a sub-peptide of 'ALTsequence' centered around 'mutpos'. If
'ALTsequence' is long enough, sub-peptide is of size 2 * 'size' - 1 for odd
'size' and of 2 * 'size' if it's even.
'''
def get_short_peptide(ALTsequence, mutpos, size):
## subtract -1 from int(mutpos) because python is 0 based
    start = max(0, mutpos + 1 - size)
    end = mutpos + size
    return ALTsequence[start:end]

# create peptides of size 'size' using a sliding window
# NOTE: can occur that for a given input peptide  different sliding windows
# contain the same short peptide sequence
def get_mers(short_pep, size):
    results = []
    for i in range(size):
        slide_peptide = short_pep[i:i + size]
        if len(slide_peptide) < size:
            # depends on python len for arrays that are indexed ouside of bounds
            # to still only output the length of the content
            break
        results.append(slide_peptide)
    return results

def gap_function(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Gap open penalty
        return -0.6 + (x/51) / 2
    return -0.6 + (x/51) / 2 - 0.1 # gap extension penalty

# assumes varID format: chrX_X_X/X, e.g., chr3_112563516_AGTATTCTGCCAAT/A
def check_if_frameshift(varID):
    split_varID = varID.split("/")
    ref = split_varID[0].split("_")[-1]
    alt = split_varID[-1]
    indel_len = abs(len(ref) - len(alt))
    if indel_len % 3 != 0:
        print(varID, "is frameshift")
        return True
    return False

'''
Returns a tuple with the subsequence of specified 'size' centered (when possible)
around the mutation and its starting index on the input alt sequence, 'ALTseq'.
This index is useful for obtaining the Wild Type subsequence aligned at the start
of the alt subsequence.
'''
def take_sub_peptide(ALTseq, pos, size, is_frameshift):
    alt_len = len(ALTseq)
    if pos >= alt_len or pos < 0:
        incorrect_mut_pos = ("Mutation position is not valid, value is "
                             "negative or bigger than ALT sequence length. "
                             "The value of pos was: {} and length: {}").format(pos, alt_len - 1)
        raise Exception(incorrect_mut_pos)
    if alt_len <= size:
        print("too small of seq")
        return (ALTseq, 0)
    if alt_len > size:
        half_size = floor(size/2)
        if pos <= half_size:
            print("cut on left")
            if is_frameshift:
                return (ALTseq, 0)
            return (ALTseq[: size], 0)
        elif alt_len - pos - 1 <= half_size:
            print("cut on right")
            return (ALTseq[alt_len - size:], alt_len - size)
        else:
            print("centered")
            if is_frameshift:
                return (ALTseq[pos - half_size:], pos - half_size)
            return (ALTseq[pos - half_size:pos + half_size + 1], pos - half_size)

### MAIN ###
def run(sizes):

    mass_spec_folder_name = "51aa/"
    if not os.path.exists(mass_spec_folder_name):
        os.mkdir(mass_spec_folder_name)
    mass_spec_file_name = mass_spec_folder_name + sample_name + "_2massSpec" + "_" + database + ".txt"
    mass_spec_file = open(mass_spec_file_name, "w+")
    mass_spec_file.write("varID\tALTpeptide_all_aa\n")
    mass_spec_written = False

    reactivity_25aa_folder_name = "25aa/"
    if not os.path.exists(reactivity_25aa_folder_name):
        os.mkdir(reactivity_25aa_folder_name)
    reactivity_25aa_file_name = reactivity_25aa_folder_name + sample_name + "_25aa_" + database + ".txt"
    reactivity_25aa_file = open(reactivity_25aa_file_name, "w+")
    reactivity_25aa_file.write("varID\tWTpeptide\tALTpeptide_25aa\n")
    reactivity_25aa_written = False


    for size in sizes:
        blast_peps_file_name = sample_name + "_blasted_peptides_" + str(size) + "mers" + "_" + database + ".fsa"
        blast_peps_file = open(blast_peps_file_name, "w+")
        blast_peps = []

        # NOTE: "If you just want to read or write one file see open()"
        for line in fileinput.input(input_file):
            words = line.split()
            if not words[0].startswith("chr"):
                continue
            varID      = words[0]
            WTpeptide  = words[1]
            ALTpeptide = words[2]
            # Global alignment with gap penalty of -0.5 (and -0.1 for extending it)
            # The 'x' default means the mismastch score is 0 and match is 1
            alignments = pairwise2.align.globalxc(WTpeptide, ALTpeptide, gap_function, gap_function, penalize_end_gaps=(False, True))
            print("---", varID)
            for ali in alignments:
                print(format_alignment(*ali))
            mutpos = []
            start_WT = False
            start_ALT = False
            ALT_offset = 0
            WT_offset = 0
            is_frameshift = check_if_frameshift(varID)
            found_first_mut = False
            for i, c in enumerate(zip(alignments[0][0][0:len(ALTpeptide)], alignments[0][1])):
                if found_first_mut and is_frameshift:
                    for x in range(i + ALT_offset, len(ALTpeptide) + ALT_offset):
                        mutpos.append(x)
                    break
                if c[0] != "-":
                    start_WT = True
                if c[1] != "-":
                    start_ALT = True
                if not start_WT and c[0] == "-":
                    WT_offset -= 1
                if not start_ALT and c[1] == "-":
                    ALT_offset -= 1
                if start_WT and start_ALT and c[0] != c[1]:
                    found_first_mut = True
                    mutpos.append((i + ALT_offset, i + WT_offset))
                    print("Appended new mutpos", mutpos)

            mass_spec_suffix = []
            # only slide window over an ALTpeptide if it had differences with WTpeptide
            for e in mutpos:
                if isinstance(e, tuple):
                    i = e
                elif isinstance(e, int):
                    i = (e, None)
                else:
                    raise Exception("mutpos can only be of type int or tuple: {}".format(i))
                print("inside for loop", i, mutpos)
                short_pep = get_short_peptide(ALTpeptide, i[0], size)
                mers = get_mers(short_pep, size)
                # print ("%s\t%s\t%s\t%s\t%s" % (varID, WTpeptide, ALTpeptide, mutpos, mers))
                # NOTE: allows repeated peptides, thus later blasted more than once
                mass_spec_suffix.append(str(i[0] + 1)) # mass spec starts counting at 1
                for mer in mers:
                    blast_peps_file.write(">" + varID + "\n" + mer + "\n")
                    blast_peps.append(mer)

                if not reactivity_25aa_written:
                    alt_pep_25aa, alt_pep_25aa_start = take_sub_peptide(ALTpeptide, i[0], 25, is_frameshift)
                    # NOTE: 'take_sub_peptide' can be made to return 'wt_pep_25aa' directly
                    wt_pep_25aa = WTpeptide[alt_pep_25aa_start - ALT_offset + WT_offset:\
                                        alt_pep_25aa_start - ALT_offset + WT_offset + 25]

                    reactivity_25aa_file.write(varID + "\t" + wt_pep_25aa \
                                                     + "\t" + alt_pep_25aa \
                                                     + "\t" + str(i) + " " + str(is_frameshift) + "\n")

            if not mass_spec_written:
                mass_spec_file.write(varID + "_M" + ",".join(mass_spec_suffix) + "\t" + ALTpeptide + "\n")

        blast_peps_file.close()
        mass_spec_file.close()
        reactivity_25aa_file.close()
        mass_spec_written = True
        reactivity_25aa_written = True

        ### BLAST ###
        if database == "RS":
            db_path = "/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/refseq_taxid_9606/GRCh38_latest_protein"
        elif database == "SP":
            db_path = "/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606"
        else:
            raise Exception("Unknown database: {}".format(database))
        # TODO: create user param for DB and create file name for XML file accordingly
        xml_file_name = sample_name + "_blast_output_" + str(size) + "mers_" + database + ".xml"
        blastp_cline = NcbiblastpCommandline( \
                        query = blast_peps_file_name, \
                        task = "blastp-short", \
                        db = db_path, \
                        outfmt = 5, \
                        out = xml_file_name, \
                        remote = False)
        print("*****    BLAST   *****")
        print("Executing:\n\t", blastp_cline)
        stdout, stderr = blastp_cline()
        print("stdout: ", stdout)
        print("stderr: ", stderr)
        print("***** END OF BLAST *****")

        result_handle = open(xml_file_name)
        blast_records = NCBIXML.parse(result_handle)

        neoantigen_candidates = {}

        decoding_file_name = sample_name + "_peptides_" + str(size) + "mers_decoding" + "_" + database + ".txt"
        decoding_file = open(decoding_file_name, "w+")

        perfect_hits_file_name = sample_name + "_perfect_hits_" + str(size) + "mers" + "_" + database + ".txt"
        perfect_hits_file = open(perfect_hits_file_name, "w+")

        for i, blast_record in enumerate(blast_records):
            print("/////////////////////////////", blast_record.query)
            print(i, " ", blast_peps[i])
            candidate = True
            for alignment in blast_record.alignments:
                print("****Alignment****")
                print("sequence:", alignment.title)
                for hsp in alignment.hsps:
                    print(hsp.query)
                    print(hsp.match)
                    print(hsp.sbjct)
                    print("Identities: ", hsp.identities)
                    print("Alignment Lenght: ", hsp.align_length)
                    print("Size: ", size)
                    # Disregard perfect hits. `align_length` is also checked because of
                    # possible gaps in the alignment.
                    if (size == hsp.identities == hsp.align_length):
                        print("FOUND PERFECT HIT:\t" + blast_peps[i])
                        candidate = False
                        perfect_hits_file.write(blast_peps[i] + "\t" + blast_record.query + "\t" + alignment.title + "\n")
                        break
                if not candidate:
                    break
            if candidate:
                decoding_file.write(str(i) + "\t" + blast_record.query + "\n")
                try:
                    neoantigen_candidates[blast_peps[i]].add(i)
                except KeyError: # key not in dictionary
                    neoantigen_candidates[blast_peps[i]] = {i}

        decoding_file.close()
        perfect_hits_file.close()

        output_file_name = sample_name + "_peptides_" + str(size) + "mers" + "_" + database + ".fsa"
        # Note that the same peptide generated from two different variants will be repeated
        with open(output_file_name, "w+") as output_file:
            for pep, varIDs in neoantigen_candidates.items():
                for _id in varIDs:
                    output_file.write(">" + str(_id) + "\n" + pep + "\n")


if __name__ == '__main__':
    arguments = docopt(__doc__, version='Gen Inputs Bind Pred 1.0')
    sizes = [int(i) for i in arguments["-p"].split(",")]
    input_file = arguments["-f"]
    sample_name = os.path.splitext(os.path.basename(input_file))[0]
    database = arguments["-d"]
    run(sizes)
