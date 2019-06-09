"""
This script is meant to generate the inputs for the binding prediction phase (e.g netMHC)

Usage:
  gen_input_binding_prediction.py -f NIC4.25peptide.txt -p 8 -b blast_peps

Options:
  -h --help        Show this screen.
  --version        Show version.
  -p=8,9,10,11,12  Peptide length to be generated
  -f=<file>        Three column file containing varID, WTpeptide and ALTpeptide. ALT peptide has been currated for phases variants (using ex. ISOVAR)
  -b=file_name     File name to write the peptides that are blasted (FASTA format)
"""

import fileinput
import sys
import re

from docopt import docopt

from Bio import pairwise2
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
from Bio.Blast import NCBIXML

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

### MAIN ###
blast_peps_file = open(blast_peps_file_name, "w+")
blast_peps = []

# NOTE: "If you just want to read or write one file see open()"
for line in fileinput.input(file):
    words = line.split()
    if not words[0].startswith("chr"):
        continue
    varID      = words[0]
    WTpeptide  = words[1]
    ALTpeptide = words[2]
    # Global alignment with gap penalty of -0.5 (and -0.1 for extending it)
    # The 'x' default means the mismastch score is 0 and match is 1
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
        # NOTE: allows repeated peptides, thus later blasted more than once
        for mer in mers:
            blast_peps_file.write(">" + varID + "\n" + mer + "\n")
            blast_peps.append(mer)

blast_peps_file.close()

### BLAST ###
xml_file_name = "swissprot9606_1.xml" # NOTE: turn this non-static
blastp_cline = NcbiblastpCommandline( \
                query=blast_peps_file_name, \
                task="blastp-short", \
                db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606", \
                outfmt=5, \
                out=xml_file_name, \
                remote=False)
print("Executing:\n\t", blastp_cline)
stdout, stderr = blastp_cline()
print("stdout: ", stdout)
print("stderr: ", stderr)

result_handle = open(xml_file_name)
blast_records = NCBIXML.parse(result_handle)

neoantigen_candidates = {}

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
                break
        if not candidate:
            break
    if candidate:
        try:
            neoantigen_candidates[blast_peps[i]].add(blast_record.query)
        except KeyError: # key not in dictionary
            neoantigen_candidates[blast_peps[i]] = {blast_record.query}


output_file_name = "inputs_bind_pred.fasta" # NOTE: turn this non-static
# writing to file is done after having discarded duplicate peptides
with open(output_file_name, "w+") as output_file:
    for pep, varIDs in neoantigen_candidates.items():
        output_file.write(">" + ", ".join(varIDs) + "\n" + pep + "\n")
