"""
Predict peptide binding to MHC

Usage:
  binding_prediction.py -f long_peps.txt -s 8,9 -d /path/to/db -l run_1.log

Options:
  -h --help           Show this screen.
  --version           Show version.
  -s 8,9,10,11,12     Peptide length(s) to be generated.
  -f <file>           Three column file containing varID, WTpeptide and ALTpeptide (mutant peptide).
                      Mutant peptide has been curated for phased variants (using e.g. isovar).
  -d <path/to/db>     Path to database.
  -l <log_file_name>  Name to give to the log file.
"""

import fileinput
import sys
import re
import os
from math import log, floor
import logging

from docopt import docopt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
from Bio.Blast import NCBIXML

import seq_helpers as helper

"""
Returns scores for no gap, gap open and gap extension used in aligning 2
sequences. Here we penalize more heavily opening gaps more to the right of the
sequence rather than as soon as possible.
'x' is the gap position in the sequence and 'y' is gap length.
"""
def gap_function(x, y):
    if y == 0:  # no gap
        return 0
    elif y == 1:  # gap open penalty
        return -0.6 + (x/51) / 2
    return -0.6 + (x/51) / 2 - 0.1 # gap extension penalty

def compute_short_peptides():
    pass

def blast_peptides():
    pass

def predict_binding():
    pass

def main():
    arguments = docopt(__doc__, version='Binding Prediction a.0.1')
    logging.basicConfig(filename=arguments["-l"], \
                        filemode='w', \
                        format="%(asctime)s %(levelname)s: %(message)s", \
                        level=logging.DEBUG) # INFO, WARNING, ERROR, CRITICAL
    sizes = [int(i) for i in arguments["-s"].split(",")]
    input_file = arguments["-f"]
    sample_name = os.path.splitext(os.path.basename(input_file))[0]
    database = arguments["-d"]
    for size in sizes:
        compute_short_peptides()
        blast_peptides()
        predict_binding()

if __name__ == '__main__':
    main()
