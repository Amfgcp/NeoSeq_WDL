# NeoSeq_WDL
## Binding prediction installation and test
1- (Optional) Create & activate conda environment with python 3.6 and docopt (newer python versions may work, but were not tested):
    
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
