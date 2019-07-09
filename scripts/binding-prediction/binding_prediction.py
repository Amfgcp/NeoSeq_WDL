"""
Predict peptide binding affinity to MHC

Usage:
  binding_prediction.py -f long_peps.txt -s 8,9 -d /path/to/db -l run_1.log

Options:
  -h --help           Show this screen.
  --version           Show version.
  -s 8,9,10,11,12     Peptide length(s) to be generated.
  -f <file>           Three column file containing var ID, WT peptide (wild type) and MUT peptide (mutant).
                      Mutant peptide has been curated for phased variants (using e.g. isovar).
  -d <path/to/db>     Path to database.
  -l <log_file_name>  Name to give to the log file.
"""

import fileinput
import sys
import re
import os
from math import floor
import logging

from docopt import docopt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
from Bio.Blast import NCBIXML
from mhctools import NetMHC
from mhctools import NetMHCpan

import seq_helpers as helper

"""
Takes the input file with the WT and MUT peptide sequences, finds the mutation
positions, breaks the sequences into smaller ones according to 'size' and
outputs them: WT as a dictionary indexed by a generated id, 'fasta_id', and
MUT also as the latter but additionally as a file name (useful for blast) to a
file containing these smaller peptides as well.
"""
def compute_short_peptides_from_file(pep_file, size):
    logging.info("Analyzing: %s", pep_file)
    MUT_peps_to_blast_file_name = \
         SAMPLE + "_peptides_to_blast_" + str(size) + "mers" + "_" + DB + ".fsa"
    MUT_peps_to_blast_file = open(MUT_peps_to_blast_file_name, "w+")
    WT_peps = dict()
    MUT_peps = dict()
    line_num = 1
    data_25aa = []
    data_mass_spec = []
    debug_count = 0
    for line in fileinput.input(pep_file):
        words = line.split()
        if not words[0].startswith("chr"):
            logging.warning("Skipping line in input file:\n%s", line)
            continue
        var_id  = words[0]
        WT_pep  = words[1]
        MUT_pep = words[2]
        # using default scores for match and mismatch, 1 and 0, respectively
        alignments = pairwise2.align.globalxc(WT_pep, MUT_pep, \
           gap_function, gap_function, penalize_end_gaps=(False, True))
        for ali in alignments:
            logging.debug("Alignment(s):\n%s", format_alignment(*ali))

        is_frameshift = helper.check_if_frameshift(var_id)
        MUT_WT_pos, WT_offset, MUT_offset = find_mutation_positions(var_id, \
                                     alignments, WT_pep, MUT_pep, is_frameshift)

        mass_spec_suffix = []
        recorded_25aa_peps_once = False
        for e in MUT_WT_pos:
            MUT_short_pep, WT_short_pep = helper.take_sub_peptide(MUT_pep, \
                       WT_pep, e[0], size * 2 - 1, False, MUT_offset, WT_offset)
            WT_mers = helper.calc_mers(WT_short_pep, size)
            MUT_mers = helper.calc_mers(MUT_short_pep, size)
            logging.debug(len(WT_mers) == len(MUT_mers))
            # logging.debug("%s\t%s\t%s\t%s\t%s" % (var_id, WT_pep, MUT_pep, MUT_WT_pos, MUT_mers))
            mass_spec_suffix.append(str(e[0] + 1)) # mass spec starts count at 1
            mer_num = 0
            for mer in zip(WT_mers, MUT_mers):
                fasta_id = var_id + "-L" + str(line_num) + "-WT" +  str(e[1]) \
                                  + "-MUT" +  str(e[0]) + "-M" + str(mer_num)
                WT_peps[fasta_id] = mer[0]
                MUT_peps[fasta_id] = mer[1]
                MUT_peps_to_blast_file.write(">" + fasta_id + "\n" + mer[1] + "\n")
                mer_num += 1
                debug_count += 1

            if recorded_25aa_peps_once and is_frameshift:
                continue # such that these kind of repeats aren't recorded
            if not WRITTEN_25AA_REACTIVITY:
                MUT_pep_25aa, WT_pep_25aa = helper.take_sub_peptide(MUT_pep, \
                         WT_pep, e[0], 25, is_frameshift, MUT_offset, WT_offset)
                data_25aa.append((var_id, WT_pep_25aa, MUT_pep_25aa))
                recorded_25aa_peps_once = True

        if not WRITTEN_MASS_SPEC:
            data_mass_spec.append((var_id, mass_spec_suffix, MUT_pep))
        line_num += 1

    MUT_peps_to_blast_file.close()

    if not WRITTEN_25AA_REACTIVITY:
        write_25aa_reactivity_file(data_25aa)
    if not WRITTEN_MASS_SPEC:
        write_mass_spec_file(data_mass_spec)

    logging.info("Short peptides calculated from: %s", pep_file)
    return (WT_peps, MUT_peps, MUT_peps_to_blast_file_name)

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

"""
Creates a folder where the file for mass spec is written. The input consists of
an array of tuples with the data needed. Example of line in file:
"chr22_50508443_G/A_M23,41    PIRDLTKNWEVDVAAQLGEYLEKLDQICISFDEGKTTMNFLEAALLIQ"
"""
def write_mass_spec_file(data_mass_spec):
    logging.info("Writing Mass Spec file")
    global WRITTEN_MASS_SPEC
    folder_name = "51aa/"
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    file_name = folder_name + SAMPLE + "_2massSpec" + "_" + DB + ".txt"
    file = open(file_name, "w+")
    file.write("var_id\tMUT_peptide\n")

    for var_id, mass_spec_suffix, MUT_pep in data_mass_spec:
        file.write(var_id + "_M" + ",".join(mass_spec_suffix) + "\t" + MUT_pep + "\n")

    WRITTEN_25AA_REACTIVITY = True
    file.close()
    logging.info("Finished writing Mass Spec file")

"""
Creates a folder where the file for peptide reactivity is written. The input
consists of an array of tuples with the data needed. Example of line in file:
"chr9_87651644_C/A    VDYQDRHGNTPLHVACKDGNMPIVV     VDYQDRHGNTPLSVACKDGNMPIVV"
"""
def write_25aa_reactivity_file(data_25aa):
    logging.info("Writing 25aa reactivity file")
    global WRITTEN_25AA_REACTIVITY
    folder_name = "25aa/"
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    file_name = folder_name + SAMPLE + "_25aa_" + DB + ".txt"
    file = open(file_name, "w+")
    file.write("var_id\tWT_peptide\tMUT_peptide\n")

    for var_id, WT_pep_25aa, MUT_pep_25aa in data_25aa:
        file.write(var_id + "\t" + WT_pep_25aa + "\t" + MUT_pep_25aa + "\n")

    WRITTEN_25AA_REACTIVITY = True
    file.close()
    logging.info("Finished writing 25aa reactivity file")

"""
Returns a tuple of: a tuple array with the adjusted mutation positions in
'MUT_pep' and 'WT_pep' (in this order) and the 2 offsets for each of these
sequences that were used to make the adjustment.
"""
def find_mutation_positions(var_id, alignments, WT_pep, MUT_pep, is_frameshift):
    MUT_WT_pos = []
    WT_started = False
    MUT_started = False
    WT_offset = 0
    MUT_offset = 0
    found_first_mut = False
    for e, c in enumerate(zip(alignments[0][0][0:len(MUT_pep)], alignments[0][1])):
        if found_first_mut and is_frameshift:
            logging.debug("Frameshift case so append remaining positions in \
                     MUT: %i-%i", e + MUT_offset, len(MUT_pep) + MUT_offset - 1)
            for pos in range(e + MUT_offset, len(MUT_pep) + MUT_offset):
                MUT_WT_pos.append((pos, None))
            break
        if c[0] != "-":
            WT_started = True
        if c[1] != "-":
            MUT_started = True
        if not WT_started and c[0] == "-":
            WT_offset -= 1
        if not MUT_started and c[1] == "-":
            MUT_offset -= 1
        if WT_started and MUT_started and c[0] != c[1]:
            found_first_mut = True
            MUT_WT_pos.append((e + MUT_offset, e + WT_offset))
            logging.debug("Identified mutation positions, WT: %i MUT: %i", \
                                              e + WT_offset, e + MUT_offset)
    return (MUT_WT_pos, WT_offset, MUT_offset)


"""
Blast peptides contained in the file with name 'peps_file_name'.
Returns file name of blast output (XML).
"""
def blast_peptides(peps_file_name, size):
    logging.info("Blasting peptides from file: %s", peps_file_name)
    # NOTE: turn 'db_path' to user input
    if DB == "RS":
        db_path = "/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/refseq_taxid_9606/GRCh38_latest_protein"
    elif DB == "SP":
        db_path = "/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606"
    else:
        raise Exception("Could not use database: {}".format(DB))

    xml_file_name = SAMPLE + "_blast_output_" + str(size) + "mers_" + DB + ".xml"
    blastp_cline = NcbiblastpCommandline( \
                    query = peps_file_name, \
                    task = "blastp-short", \
                    db = db_path, \
                    outfmt = 5, \
                    out = xml_file_name, \
                    remote = False)

    logging.debug("Executing: %s", blastp_cline)
    stdout, stderr = blastp_cline()
    logging.info("Finished blasting file: %s", peps_file_name)
    logging.info("stdout: %s", stdout)
    logging.info("stderr: %s", stderr)

    return xml_file_name

"""
Disregard peptides found to match perfectly in terms of alignment and 'size'.
Returns dictionary of non-perfect matching peptides.
"""
def filter_db_perfect_matches(MUT_peps, blast_xml_out_file_name, size):
    results_handle = open(blast_xml_out_file_name)
    blast_records = NCBIXML.parse(results_handle)
    neoantigen_candidates = dict()

    for blast_record in blast_records:
            logging.debug("////// Blasting peptide with id: %s //////", blast_record.query)
            logging.debug("////// Peptide: %s", MUT_peps[blast_record.query])
            candidate = True
            for alignment in blast_record.alignments:
                logging.debug("*** Alignment ***")
                logging.debug("Aligned (in db) to: %s", alignment.title)
                for hsp in alignment.hsps:
                    logging.debug(hsp.query)
                    logging.debug(hsp.match)
                    logging.debug(hsp.sbjct)
                    logging.debug("Identities: %s", hsp.identities)
                    logging.debug("Alignment Length: %s", hsp.align_length)
                    logging.debug("Size: %i", size)
                    if (size == hsp.identities == hsp.align_length):
                        logging.debug("Found perfect match: " + MUT_peps[blast_record.query])
                        candidate = False
                        break
                if not candidate:
                    break
            if candidate:
                neoantigen_candidates[blast_record.query] = MUT_peps[blast_record.query]

    results_handle.close()
    return neoantigen_candidates

"""
"""
def predict_binding():
    pass

SAMPLE = "unset-sample"
DB = "unset-db"
WRITTEN_25AA_REACTIVITY = False
WRITTEN_MASS_SPEC = False
def main():
    arguments = docopt(__doc__, version='Binding Prediction a.0.1')
    global SAMPLE, DB
    logging.basicConfig(filename=arguments["-l"], \
                        filemode='w', \
                        format="%(asctime)s %(levelname)s: %(message)s", \
                        level=logging.DEBUG) # INFO, WARNING, ERROR, CRITICAL
    sizes = [int(e) for e in arguments["-s"].split(",")]
    input_file = arguments["-f"]
    SAMPLE = os.path.splitext(os.path.basename(input_file))[0] # NOTE: needed?
    DB = arguments["-d"]
    for size in sizes:
        WT_peps, MUT_peps, MUT_peps_file_name = \
            compute_short_peptides_from_file(input_file, size)
        blast_out_xml_file_name = blast_peptides(MUT_peps_file_name, size)
        neoantigen_candidates = \
            filter_db_perfect_matches(MUT_peps, blast_out_xml_file_name, size)
        predict_binding()

if __name__ == '__main__':
    main()
