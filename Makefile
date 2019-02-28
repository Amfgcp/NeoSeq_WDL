PICARD = /home/druano/tools/picard/build/libs/picard.jar

REF = /home/druano/tools/gatk-bundle-b38/Homo_sapiens_assembly38.fasta
data = /exports/path-demiranda/usr/amfgcp/synthetic-data/
TMP_DIR = `mktemp`
dream_data = /mnt/patharchief/ImmunoGenomics/ngsdata/synthetic/DREAMchallenge

samples = set1.normal.v2 set1.tumor.v2 set2.normal set2.tumor set3.normal set3.tumor
library = IDT

#=============================================================
# Realign bam files with hg38
#http://crazyhottommy.blogspot.com/2015/10/convert-bam-to-fastq-and-remap-it-with.html
#=============================================================
createBAMs: $(patsubst %, $(WORKDIR)/bam/%.bam, $(samples))


$(WORKDIR)/bam/%.bam:
	 $(eval TMP_DIR := $(shell mktemp -d))
	java -Xmx32G -jar $(PICARD) SamToFastq \
		I=$(dream_data)/synthetic.challenge.$*.bam \
		INTERLEAVE=true NON_PF=true \
		FASTQ=/dev/stdout TMP_DIR=$(TMP_DIR) | \
	/home/druano/tools/bwa-0.7.16a/bwa mem -t 10 -M -p \
	-R '@RG\tID:$(shell echo $* | sed "s/_R1//")\tLB:$(library)\tSM:$(shell echo $* | \
		sed "s/_L.*_R1//")\tPL:ILLUMINA' $(REF) /dev/stdin - | \
	java -Xmx3g -jar $(PICARD) SortSam I=/dev/stdin \
		SO=coordinate CREATE_INDEX=true QUIET=TRUE \
		O=$(data)/bam/$*.bam

#samtools flagstat /mnt/patharchief/ImmunoGenomics/ngsdata/synthetic/DREAMchallenge/synthetic.challenge.set1.normal.v2.bam
#915782532 + 146511378 in total (QC-passed reads + QC-failed reads)
#0 + 0 secondary
#0 + 0 supplementary
#61595565 + 12843134 duplicates
#877143430 + 106023832 mapped (95.78%:72.37%)
#915782532 + 146511378 paired in sequencing
#457891266 + 73255689 read1
#457891266 + 73255689 read2
#858891201 + 98934868 properly paired (93.79%:67.53%)
#864556748 + 100183220 with itself and mate mapped
#12586682 + 5840612 singletons (1.37%:3.99%)
#4540490 + 1132704 with mate mapped to a different chr
#3324847 + 593857 with mate mapped to a different chr (mapQ>=5)
