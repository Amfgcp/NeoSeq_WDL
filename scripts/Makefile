SHELL := /bin/bash # define which shell to use!
PICARD = /home/druano/tools/picard/build/libs/picard.jar

REF = /home/druano/tools/gatk-bundle-b38/Homo_sapiens_assembly38.fasta
data = /exports/path-demiranda/usr/amfgcp/synthetic-data
TMP_DIR = `mktemp`
dream_data = /mnt/patharchief/ImmunoGenomics/ngsdata/synthetic/DREAMchallenge
WORKDIR = /home/druano/NeoSeq

samples = set1.normal.v2 set1.tumor.v2 set2.normal set2.tumor set3.normal set3.tumor
tumors = $(shell tr " " "\n"  < <(echo $(samples)) | sed "s/N$$//g" | sed "s/T$$//g" | sort -u -V | tr "\n" " ")
library = IDT
	
#======================================================================================
# Realign bam files with hg38
#http://crazyhottommy.blogspot.com/2015/10/convert-bam-to-fastq-and-remap-it-with.html
#======================================================================================
FASTQ: $(patsubst %, $(dream_data)/fastq/%_R1.fq, $(samples))

$(dream_data)/fastq/%_R1.fq:
	java -jar -Xmx8g $(PICARD) SamToFastq QUIET=true VALIDATION_STRINGENCY=LENIENT \
		INPUT=$(dream_data)/synthetic.challenge.$*.bam \
		FASTQ=$(dream_data)/fastq/$*_R1.fq SECOND_END_FASTQ=$(dream_data)/fastq/$*_R2.fq

createSAMs: $(patsubst %, $(WORKDIR)/bam/%.sam, $(samples))
$(WORKDIR)/bam/%.sam: #$(dream_data)/fastq/%_R1.fq
	/home/druano/tools/bwa-0.7.16a/bwa mem -t 10 -M \
		-R '@RG\tID:$(shell echo $* | sed "s/_R1//")\tLB:$(library)\tSM:$(shell echo $* | sed "s/_L.*_R1//")\tPL:ILLUMINA' \
		-o $@ $(REF) $(dream_data)/fastq/$*_R1.fq $(dream_data)/fastq/$*_R2.fq

createBAMs: $(patsubst %, $(data)/bam/%.bam, $(samples))

#	$(eval TMP_DIR := $(shell mktemp -d))
$(data)/bam/%.bam: $(WORKDIR)/bam/%.sam
	java -Xmx8g -jar $(PICARD) SortSam I=$(WORKDIR)/bam/$*.sam \
		SO=coordinate CREATE_INDEX=true QUIET=TRUE TMP_DIR=/home/druano/NeoSeq/tmp_picard \
		O=$(data)/bam/$*.bam

#generate a track IGV that can be used to plot coverage
IGVtrack: $(patsubst %, $(data)/bam/%.bam.tdf, $(samples))
$(data)/bam/%.bam.tdf:
	java -jar ~/tools/IGVTools/igvtools.jar count -f min,max,mean $(data)/bam/$*.bam $(data)/bam/$*.bam.tdf /home/druano/tools/gatk-bundle-b38/Homo_sapiens_assembly38.fasta



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

#=============================================================
# Mark Duplicate Reads per sample
#=============================================================
dedup: $(patsubst %, %.dedup.bam, $(samples))
%.dedup.bam:
	@echo " --------- remove duplicated reads per sample --------- "
	java -jar $(PICARD) MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		QUIET=TRUE \
		ASSUME_SORTED=true \
		CREATE_INDEX=true \
		REMOVE_DUPLICATES=false \
		METRICS_FILE=$*.dupl.metrics \
		INPUT=$(WORKDIR)/bam/$*.bam \
		OUTPUT=$(WORKDIR)/bam/$*.dedup.bam

#======================================================================================
# Calling variants with MuTect
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php
#======================================================================================
MUT: $(patsubst %, $(WORKDIR)/muTect/%_muTect.vcf, $(tumors))	


$(WORKDIR)/muTect/%_raw.vcf: 
	/home/druano/tools/gatk-4.beta.3-SNAPSHOT/gatk-launch Mutect2 \
		-R $(REF) \
		-I $(data)/bam/$*T.dedup.recal.bam \
		-I $(data)/bam/$*N.dedup.recal.bam \
		-tumor $*T \
		-normal $*N \
		--dbsnp $(dbSNP) \
		-O $(WORKDIR)/muTect/$*_raw.vcf

$(WORKDIR)/muTect/%_muTect.vcf: $(WORKDIR)/muTect/%_raw.vcf
	/home/druano/tools/gatk-4.beta.3-SNAPSHOT/gatk-launch FilterMutectCalls \
		--variant $(WORKDIR)/muTect/$*_raw.vcf \
		--uniqueAltReadCount 4 \
		--showHidden \
		--output $(WORKDIR)/muTect/$*_muTect.vcf


#=============================================================
# Filter variants
# - 8 reads minimum in tumor, 
# - at least 4 ALT reads
# - 10% variant reads 
# - at least one ALT read in each direction
#=============================================================
#minALTfreq = 0.1
minALTfreq = 0.02
minALTfreq_inNormal = $(shell awk 'BEGIN{print $(minALTfreq)/2} ')
filterMUT: $(patsubst %, $(WORKDIR)/tmp/%_muTect.filter.vcf.gz, $(tumors))

$(WORKDIR)/tmp/%_muTect.filter.vcf.gz: $(WORKDIR)/muTect/%_muTect.vcf
	@echo " --------- Filter variants --------- "
	cat $(WORKDIR)/muTect/$*_muTect.vcf | awk -F"\t" '$$7~"PASS$$" || $$1~"^#" { if ($$1!~"^#") {$$3=$$1"_"$$2"_"$$4"/"$$5} ; print }' OFS="\t" | \
		sed "s/\tFORMAT.*/\tFORMAT\tNORMAL\tTUMOR/" | \
		$(vcftools)/vcf-sort > $(WORKDIR)/muTect/$*_muTect.filter.vcf
	@echo " --------- compress and index the file --------- "
	bgzip -c $(WORKDIR)/muTect/$*_muTect.filter.vcf > $(WORKDIR)/tmp/$*_muTect.filter.vcf.gz
	tabix -p vcf -f $(WORKDIR)/tmp/$*_muTect.filter.vcf.gz

# - Combine different vcf files
#=============================================================
mergeVar: $(patsubst %, $(WORKDIR)/vcf/%_CombineVariants.vcf, $(tumors)) 

$(WORKDIR)/vcf/%_CombineVariants.vcf:
	java -Xmx4g -jar /home/druano/tools/GenomeAnalysisTK-3.8/GenomeAnalysisTK.jar \
		-T CombineVariants \
		-R $(REF) \
		--genotypemergeoption UNIQUIFY \
		-V:strelka $(WORKDIR)/tmp/$*_strelka.vcf.gz \
		-V:mutect $(WORKDIR)/tmp/$*_muTect.filter.vcf.gz \
		-V:varscan $(WORKDIR)/tmp/$*_varScan.filter.vcf.gz\
		--out $(WORKDIR)/vcf/$*_CombineVariants.vcf

# - change strelka output
#=============================================================
$(WORKDIR)/tmp/%_strelka.vcf:
	@echo " --------- add GT Information - just add 0/1 for every call both in normals and tumors --------- "
#see: https://sites.google.com/site/strelkasomaticvariantcaller/home/faq
	/home/druano/tools/somaticseq-2.3.1/utilities/modify_Strelka.py -infile $(WORKDIR)/strelka/$*/results/passed.somatic.snvs.vcf -outfile $(WORKDIR)/strelka/$*.snvs.vcf
	/home/druano/tools/somaticseq-2.3.1/utilities/modify_Strelka.py -infile $(WORKDIR)/strelka/$*/results/passed.somatic.indels.vcf -outfile $(WORKDIR)/strelka/$*.indels.vcf
	@echo " --------- merge SNVs and INDELs in a single file --------- "
	@echo " --------- fix variant id --------- "
	$(vcftools)/vcf-concat $(WORKDIR)/strelka/$*.indels.vcf $(WORKDIR)/strelka/$*.snvs.vcf | \
	awk -F"\t" '$$7~"PASS" || $$1~"^#CHROM" { if ($$1!~"^#") {$$3=$$1"_"$$2"_"$$4"/"$$5} ; print }' OFS="\t"	| \
	cat ~/Dropbox/NeoSeq/STRELKAheader.txt - | sed "s/OVERLAP;//g" | $(vcftools)/vcf-sort > $(WORKDIR)/strelka/$*_strelka.vcf
	@echo " --------- compress and index the file --------- "
	bgzip -c $(WORKDIR)/strelka/$*_strelka.vcf > $(WORKDIR)/tmp/$*_strelka.vcf.gz
	tabix -p vcf -f $(WORKDIR)/tmp/$*_strelka.vcf.gz


#filter varDict to have less false positives @TONI, can you do this in PYTHON? No, is not a valid answer :-)
#=============================================================
varDict_filter: 
# Note that we decided not to put a filter on var frequency but instead output it to the result file, so that the user can actually filter them at a later stage (&& (S$$NIC7T$$AF) >= 0.05)
# ToDo: Keep only variants that have at least 1 forward ALT read and 1 reverse ALT read 
	zcat /exports/path-demiranda/usr/amfgcp/aligning-somatic/results-NIC7-alt-pon/samples/NIC7T/somatic-variantcalling/vDict-gatk4.1.2.0-pon-with-gnomAD-extended-list/NIC7T-NIC7N.vcf.gz | \
		awk -F"\t" '$$7~"PASS" || $$1~"^#CHROM" { if ($$1!~"^#") {$$3=$$1"_"$$2"_"$$4"/"$$5} ; print }' OFS="\t"	| \
		~/tools/vawk/vawk --header '{ if (I$$STATUS == "StrongSomatic" || I$$STATUS == "LikelySomatic") print}' | \
		~/tools/vawk/vawk '{ if ((S$$NIC7T$$DP) >= 8 && (S$$NIC7T$$VD) >= 4) print}' | \
		cut -f3 | sort -k1,1 | \
		join -v2 - <(cat /mnt/patharchief/ImmunoGenomics/ngsdata/VisualInspection/NIC7.txt |  awk '$$2~/T$$/ {print $$1}' | sort -k1,1) | less  