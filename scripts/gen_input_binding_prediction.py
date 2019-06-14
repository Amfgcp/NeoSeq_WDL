"""
This script is meant to generate the inputs for the binding prediction phase (e.g netMHC)

Usage:
  gen_input_binding_prediction.py -f NIC4.25peptide.txt -p 8

Options:
  -h --help        Show this screen.
  --version        Show version.
  -p=8,9,10,11,12  Peptide length to be generated
  -f=<file>        Three column file containing varID, WTpeptide and ALTpeptide. ALT peptide has been currated for phases variants (using ex. ISOVAR)
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

from math import log

def getShortPeptide(ALTsequence, mutpos, size):
## subtract -1 from int(mutpos) because python is 0 based
    start = max(0, (int(mutpos) + 1  - size))
    end   = (int(mutpos) ) + size
    result = ALTsequence[start:end]
    return result

# create peptides of size 'size' using a sliding window
# NOTE: can occur that for a given input peptide  different sliding windows
# contain the same short peptide sequence
def getMer(shortPep, size):
    results = []
    for i in range(size):
        slidePeptide = shortPep[i:i+size]
        if len(slidePeptide) < size:
            break
        results.append(slidePeptide)
    return results

def gap_function(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Gap open penalty
        return -0.6 + (x/51) / 2
    return -0.6 + (x/51) / 2 - 0.1 # gap extension penalty

### MAIN ###
def run(sizes):

    mass_spec_folder_name = "51aa/"
    if not os.path.exists(mass_spec_folder_name):
        os.mkdir(mass_spec_folder_name)
    mass_spec_file_name = mass_spec_folder_name + sample_name + "_2massSpec.txt"
    mass_spec_file = open(mass_spec_file_name, "w+")
    mass_spec_file.write("varID\tALTpeptide\n")
    mass_spec_written = False

    for size in sizes:
        blast_peps_file_name = sample_name + "_blasted_peptides_" + str(size) + "mers.fsa"
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
            # for ali in alignments:
            #     print(format_alignment(*ali))
            mutpos = []
            start_WT = False
            start_ALT = False
            ALT_offset = 0
            # NOTE: Check if 1st alignment has the best score
            for i, c in enumerate(zip(alignments[0][0][0:len(ALTpeptide)], alignments[0][1])):
                if c[0] != "-":
                    start_WT = True
                if c[1] != "-":
                    start_ALT = True
                if start_WT and start_ALT and c[0] != c[1]:
                    mutpos.append(i + ALT_offset)
                if not start_ALT and c[1] == "-":
                    ALT_offset -= 1

            mass_spec_suffix = []
            # only slide window over a ALTpeptide if this one had differences with WTpeptide
            for i in mutpos:
                shortPep = getShortPeptide(ALTpeptide, i, size)
                mers = getMer(shortPep, size)
                # print ("%s\t%s\t%s\t%s\t%s" % (varID, WTpeptide, ALTpeptide, mutpos, mers))
                # NOTE: allows repeated peptides, thus later blasted more than once
                mass_spec_suffix.append(str(i + 1)) # mass spec starts counting at 1
                for mer in mers:
                    blast_peps_file.write(">" + varID + "\n" + mer + "\n")
                    blast_peps.append(mer)

            if not mass_spec_written:
                mass_spec_file.write(varID + "_M" + ",".join(mass_spec_suffix) + "\t" + ALTpeptide + "\n")

        blast_peps_file.close()
        mass_spec_file.close()
        mass_spec_written = True

        ### BLAST ###
        xml_file_name = sample_name + "_blast_output_" + str(size) + "mers_swissprot9606.xml"
        # TODO: create param for DB and create file name for XML file accordingly
        blastp_cline = NcbiblastpCommandline( \
                        query=blast_peps_file_name, \
                        task="blastp-short", \
                        db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606", \
                        outfmt=5, \
                        out=xml_file_name, \
                        remote=False)
        print("*****    BLAST   *****")
        print("Executing:\n\t", blastp_cline)
        stdout, stderr = blastp_cline()
        print("stdout: ", stdout)
        print("stderr: ", stderr)
        print("***** END OF BLAST *****")

        result_handle = open(xml_file_name)
        blast_records = NCBIXML.parse(result_handle)

        neoantigen_candidates = {}

        decoding_file_name = sample_name + "_peptides_" + str(size) + "mers_decoding.txt"
        decoding_file = open(decoding_file_name, "w+")

        perfect_hits_file_name = sample_name + "_perfect_hits_" + str(size) + "mer.txt"
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

        output_file_name = sample_name + "_peptides_" + str(size) + "mers.fsa"
        # Note that the same peptide generated from two different variants will be repeated
        with open(output_file_name, "w+") as output_file:
            for pep, varIDs in neoantigen_candidates.items():
                for _id in varIDs:
                    output_file.write(">" + str(_id) + "\n" + pep + "\n")


if __name__ == '__main__':
    arguments = docopt(__doc__, version='Gen Inputs Bind Pred 1.0')
    sizes = [int(i) for i in arguments["-p"].split(",")]
    input_file = arguments["-f"]
    sample_name = input_file.split(".")[0]

    run(sizes)
