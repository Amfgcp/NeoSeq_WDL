import fileinput
import sys
import re
from docopt import docopt
from Bio import pairwise2
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO

### Testing Area ###
proteinSeqQuery="ATEEGDQEQDPEPEEEAAVEGEEEEEGAATAAAAPGHSAVPPPPPQLPPLP"
proteinSeqseqFound="ATEEGDQEQDPEPEEEAAVEGEEEEGAATAAAAPGHSAVPPPPPQLPPLPP"

proteinSeqQuery="ASLEQRVQMQEDDIQLLKSALADVVRRLNITEEQQAVLNRKGPTKARPLMQ"
proteinSeqseqFound="ASLEQRVQMQEDDIQLLKSALADVVWRLNITEEQQAVLNRKGPTKARPLMQ"

proteinSeqQuery="LLQPIRDLTKNWEVDVAAQLGEYLEELDQICISFDEGKTTMNFIEAALLIQ"
proteinSeqseqFound="PIRDLTKNWEVDVAAQLGEYLEKLDQICISFDEGKTTMNFLEAALLIQ"

alignments = pairwise2.align.globalxs(proteinSeqQuery, proteinSeqseqFound, -0.5, -0.1)
pos = []
start = False
for i, c in enumerate(zip(alignments[0][0][0:len(proteinSeqseqFound)], alignments[0][1])):
    if c[1] != "-":
        start = True
    if start and c[0] != c[1]:
        pos.append(i)

print(pos)
print(alignments)

'''
Example query in cmd:
    $ echo KILDAVVAQKQ | /exports/path-demiranda/usr/amfgcp/tools/ncbi-blast-2.9.0+/bin/blastp -query - -task blastp-short -db /exports/path-demiranda/usr/amfgcp/databases/blast-test/human.bdb -out test_exp_out.csv -num_threads 1 -outfmt '10 qseqidsseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'
Equivalent in BioPython:
    blastp_cline = NcbiblastpCommandline(query="exp_string.txt", task="blastp-short", db="/exports/path-demiranda/usr/amfgcp/databases/blast-test/human.bdb", out="test_exp_out.csv", outfmt="'10 qseqidsseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'",  remote=False)
'''
blastp_cline = NcbiblastpCommandline(query="exp_string.txt", task="blastp-short", db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/swissprot/swissprot_v5", outfmt=5, out="test_exp_swissprot_output.xml", remote=False)

blastp_cline = NcbiblastpCommandline(query="exp_string.txt", task="blastp-short", db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606", outfmt=5, out="test_exp_swissprot_output2.xml", remote=False)

blastp_cline = NcbiblastpCommandline(query="exp_string.txt", task="blastp-short", db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/refseq_taxid_9606/GRCh38_latest_protein", outfmt=5, out="test_exp_refseq_output2.xml", remote=False)

print("Executing:\n", blastp_cline)
stdout, stderr = blastp_cline()

for query_results in SearchIO.read("test_exp_swissprot_output2.xml", "blast-xml"):
    print('Results for query')
    print(query_results)
    # print('Results for query', query_results.id)
    # if not query_results: # no hits found
    #     print('No results found for', query_results.id)
    #     continue
    # for hit in query_results:
    #     print('\tHit', hit.id)
    #     for hsp in hit:
    #         if hsp.aln_span == hsp.ident_num:   # perfect match if alignment span is equal to identity number
    #             print ('\tPerfect match found!')