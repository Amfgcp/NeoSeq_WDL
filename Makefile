PICARD = /home/druano/tools/picard/build/libs/picard.jar

REF = /home/druano/tools/gatk-bundle-b38/Homo_sapiens_assembly38.fasta
data = /exports/path-demiranda/usr/amfgcp/synthetic-data/
TMP_DIR = `mktemp`
dream_data = /mnt/patharchief/ImmunoGenomics/ngsdata/synthetic/DREAMchallenge

samples = set1.normal.v2 set1.tumor.v2 set2.normal set2.tumor set3.normal set3.tumor
library = IDT

#======================================================================================
# Realign bam files with hg38
#http://crazyhottommy.blogspot.com/2015/10/convert-bam-to-fastq-and-remap-it-with.html
#======================================================================================
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


#======================================================================================
# Calling variants with MuTect
# https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php
#======================================================================================
MUT: $(patsubst %, $(WORKDIR)/muTect/%_muTect.vcf, $(tumors))	

doMuTect: 
	$(MAKE) -f $(makefile) $(patsubst %, $(WORKDIR)/tmp/%_muTect.filter.vcf.gz, $(tumors))



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
# previously results were filtered as below:
#		~/tools/vawk/vawk --header '{if ( (S$$TUMOR$$ALT_F1R2 + S$$TUMOR$$ALT_F2R1 + S$$TUMOR$$REF_F1R2 + S$$TUMOR$$REF_F2R1) >= 8 && \
#		(S$$TUMOR$$ALT_F1R2 + S$$TUMOR$$ALT_F2R1 >= 4)  && \
#		S$$TUMOR$$AF>= $(minALTfreq) && \
#		S$$TUMOR$$ALT_F1R2>=1 && S$$TUMOR$$ALT_F2R1>=1 ) print }' | \
#		$(vcftools)/vcf-sort > $(WORKDIR)/muTect/$*_muTect.filter.vcf
	@echo " --------- compress and index the file --------- "
	bgzip -c $(WORKDIR)/muTect/$*_muTect.filter.vcf > $(WORKDIR)/tmp/$*_muTect.filter.vcf.gz
	tabix -p vcf -f $(WORKDIR)/tmp/$*_muTect.filter.vcf.gz