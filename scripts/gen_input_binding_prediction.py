"""
This script is meant to generate the inputs for the binding prediction phase

Usage:
  gen_input_binding_prediction.py -f NIC4.25peptide.txt -p 8 -b blast_peps

Options:
  -h --help             Show this screen.
  --version             Show version.
  -p=8,9,10,11,12       Peptide length to be generated
  -f=<file>             Three column file containing varID, WTpeptide and ALTpeptide. ALT peptide has been currated for phases variants (using ex. ISOVAR)
  -b=file_name          File name to write the peptides that are blasted (FASTA format)
"""

import fileinput
import sys
import re

from docopt import docopt

from Bio import pairwise2
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Gen Inputs Bind Pred 1.0')
    size = int(arguments["-p"])
    file = arguments["-f"]
    blast_peps_file_name = arguments["-b"] + ".fasta"

#size = 8

def getShortPeptide(ALTsequence, mutpos):
## subtract -1 from int(mutpos) because python is 0 based
    start = max(0, (int(mutpos) + 1  - size))
    end   = (int(mutpos) ) + size
    result = ALTsequence[start:end]
    return result

# create peptides of size 'size' using a sliding window
def getMer(shortPep, size):
    results = []
    for i in range(size):
        slidePeptide = shortPep[i:i+size]
        if len(slidePeptide) < size:
            break
        results.append(slidePeptide)
    return results


blast_peps_file = open(blast_peps_file_name, "w+")

for line in fileinput.input(file):
    words = line.split()
    if not words[0].startswith("chr"):
        continue
    varID      = words[0]
    WTpeptide  = words[1]
    ALTpeptide = words[2]
    # Global alignment with gap penalty of -0.5 (and -0.1 for extending it)
    # By default (given by the 'x') the mismastch score is 0 and match is 1    
    alignments = pairwise2.align.globalxs(WTpeptide, ALTpeptide, -0.5, -0.1)    
    mutpos = []
    start_WT = False
    start_ALT = False
    for i, c in enumerate(zip(alignments[0][0][0:len(ALTpeptide)], alignments[0][1])):
        if c[0] != "-":
            start_WT = True
        if c[1] != "-":
            start_ALT = True
        if start_WT and start_ALT and c[0] != c[1]:
            mutpos.append(i)

    # only slide window over a ALTpeptide if this one had differences with WTpeptide
    for i in mutpos:
        shortPep = getShortPeptide(ALTpeptide, i)
        mers = getMer(shortPep, size)
        # print ("%s\t%s\t%s\t%s\t%s" % (varID, WTpeptide, ALTpeptide, mutpos, mers))
        for mer in mers:
            blast_peps_file.write(">\n" + mer + "\n")

blast_peps_file.close()

# blastp_cline = NcbiblastpCommandline(query=blast_peps_file_name, \
#                 task="blastp-short", \
#                 db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606", \
#                 outfmt=5, \
#                 out="NOW_swissprot9606_1.xml", \
#                 remote=False)
# print("Executing:\n\t", blastp_cline)
# blastp_cline()
# break # NOTE: remove, testing
