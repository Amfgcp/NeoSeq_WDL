import fileinput
import sys
import re
from docopt import docopt
from Bio import pairwise2
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
from Bio.Blast import NCBIXML

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

blastp_cline = NcbiblastpCommandline(query="exp_string.txt", task="blastp-short", db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606", outfmt="'10 qseqidsseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'", out="test_exp_swissprot_output2.csv", remote=False)

blastp_cline = NcbiblastpCommandline(query="exp_string.txt", task="blastp-short", db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606", outfmt=5, out="test_swissprot9606_1.xml", remote=False)
# blastp_cline = NcbiblastpCommandline(query="exp_peps.txt", task="blastp-short", db="/exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606", outfmt="'10 qseqidsseqid qseq qstart qend sseq sstart send length mismatch pident evalue bitscore'", out="test_exp_peps.csv", remote=False)

print("Executing:\n\t", blastp_cline)
stdout, stderr = blastp_cline()

result_handle = open("test_swissprot9606_1.xml")
blast_records = NCBIXML.parse(result_handle) # returns an iterator
# if we wanted to save the records to process more than once we must:
# blast_records = list(blast_records)

for blast_record in blast_records:
    print(blast_record)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
            print("Identities:", hsp.identities)


# for query_results in SearchIO.read("test_exp_swissprot_output2.xml", "blast-xml"):
   # print('Results for query')
   # print(query_results)
    # print('Results for query', query_results.id)
    # if not query_results: # no hits found
    #     print('No results found for', query_results.id)
    #     continue
    # for hit in query_results:
    #     print('\tHit', hit.id)
    #     for hsp in hit:
    #         if hsp.aln_span == hsp.ident_num:   # perfect match if alignment span is equal to identity number
    #             print ('\tPerfect match found!')
