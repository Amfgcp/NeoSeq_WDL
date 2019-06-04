"""My script

Usage:
  ~/Dropbox/NeoSeq/scripts/generate_SHORTpeptides.py -f /home/druano/NeoSeq/NIC4.25peptide.txt -p 8

Options:
  -h --help             Show this screen.
  --version             Show version.
  -p=8,9,10,11,12       Peptide length to be generated
  -f=<file>             Three column file containing varID, WTpeptide and ALTpeptide. ALT peptide has been currated for phases variants (using ex. ISOVAR)
"""

import fileinput
import sys
import re
from docopt import docopt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

if __name__ == '__main__':
    arguments = docopt(__doc__, version='myscript 1.0')
    size = int(arguments["-p"])
    file = arguments["-f"]


#size = 8
#file = "/exports/path-demiranda/usr/druano/scripts/1207.corrected.txt"

def getShortPeptide(ALTsequence, mutpos):
## subtract -1 from int(mutpos) because python is 0 based
    start = max(0, (int(mutpos) + 1  - size))
    end   = (int(mutpos) ) + size
    result = ALTsequence[start:end]
    return result

# create a peptides of x size using a sliding window
def getMer(shortPep, size):
    results = []
    for i in range(size):
        slidePeptide = shortPep[i:i+size]
        if len(slidePeptide) < size:
            break
        results.append(slidePeptide)
    return results

for line in fileinput.input(file):
    words      = line.split()
    varID      = words[0]
    WTpeptide  = words[1]
    ALTpeptide = words[2]
    mutpos = [i for i, c in enumerate(zip(WTpeptide, ALTpeptide)) if c[0] != c[1]]

#   only slide window over a ALTpeptide if this one had differences with WTpeptide
    for i in mutpos:
        shortPep = getShortPeptide(ALTpeptide, i)
        mers = getMer(shortPep, size)
        print ("%s\t%s\t%s\t%s\t%s" % (varID, WTpeptide, ALTpeptide, mutpos, mers))
        # print ("%s\t%s\t%s\t%s\t%s\t%s" % (varID, WTpeptide, ALTpeptide, mutpos, shortPep, mers))
        # print ("%s" % (mutpos))

proteinSeqQuery="ATEEGDQEQDPEPEEEAAVEGEEEEEGAATAAAAPGHSAVPPPPPQLPPLP"
proteinSeqseqFound="ATEEGDQEQDPEPEEEAAVEGEEEEGAATAAAAPGHSAVPPPPPQLPPLPP"

alignments = pairwise2.align.localxx(proteinSeqQuery,proteinSeqseqFound)
for a, b in zip(alignments[0][0], alignments[0][1]):
    print (a, b)
