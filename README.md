# HLA Binding Prediction Module

This standalone module was developed in the context of a neoantigen identification pipeline for the Immunogenomics group from Leiden University Medical Center.

## Options
```
  -h --help           Show options.
  --version           Show version.
  -s <size(s)>        List of length(s) for the sub-peptides to be generated and have binding affinity predicted.
  -f <file>           Three column file containing var ID, WT peptide (wild type) and MUT peptide (mutant).
  -d <path/to/db>     Path to database.
  -b <HLA(s)>         HLA alleles for the binding prediction step.
  -o <output_dir>     Output folder where files and sub-folders will be created.
  -n <sample_name>    Sample name. Will be used in the file names.
```

## Command Example
binding_prediction.py -f long_peps.txt -s 8,9 -d /path/to/db -b A*02:01,hla-a0101 -o /path/out/dir -n Sample_X

## Setup and Testing
1- (Optional) Create & activate conda environment with python 3.6 and docopt:
    
`conda create -n <ENV_NAME_HERE> python=3.6 docopt`

`conda activate <ENV_NAME_HERE> OR source activate <ENV_NAME_HERE>`

2- Setup modified mhctools:

`pip install git+https://github.com/Amfgcp/mhctools`

NOTE: mhctools requires most of the binding predictors to be installed locally (see: https://pypi.org/project/mhctools/).

3- Clone this repository:

`git clone https://github.com/Amfgcp/NeoSeq_WDL.git`

4- Create output folder, e.g.:

`mkdir test-ouput`

5- Run example

`python neoseq_wdl/scripts/binding-prediction/binding_prediction.py -n 1207.test -f neoseq_wdl/scripts/binding-prediction/1207.test.txt -s 8 -d /exports/path-demiranda/usr/amfgcp/databases/ncbi/v5/generated/swissprot_taxid_9606/swissprot_taxid_9606 -b HLA-C02:02 -o test-output`

6- Expected output folder structure
```
test-output/
  ./25aa
    1207.test_25aa_swissprot_taxid_9606.fsa
  ./51aa
    1207.test_2massSpec_swissprot_taxid_9606.fsa
  ./bind-pred
    1207.test_binding_prediction_swissprot_taxid_9606.txt
  ./blast
    1207.test_blast_output_8mers_swissprot_taxid_9606.xml
    1207.test_peptides_to_blast_8mers_swissprot_taxid_9606.fsa
  ./logs
    bind-pred_1207.test_12-08-2019.log
```
