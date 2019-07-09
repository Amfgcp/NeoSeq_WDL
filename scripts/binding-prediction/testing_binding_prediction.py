from mhctools import NetMHCpan
from mhctools import NetMHC

# 1207net = HLA-A0201,HLA-A2402,HLA-B2705
# 1207pan = HLA-A02:01,HLA-A24:02,HLA-B27:05,HLA-B44:05,HLA-C02:02

def run_netMHCpan():
    # Run NetMHCpan for alleles HLA-A*01:01 and HLA-A*02:01
    predictor = NetMHCpan(alleles=["A*02:01", "hla-a0101"])
    predictor = NetMHC(alleles=["A*02:01", "hla-a0101"])

    # scan the short proteins 1L2Y and 1L3Y for epitopes
    protein_sequences = {
      "chr3_52533601_ACGCCGGGGCAGTGGG/A-L6-WT25-MUT25-M1": "NLYIQSSGRPPPSWLKDGGP",
      "chr8_127738959_G/A-1": "NLYIQSSGRPPPSWLKDGGP",
      "1L2Y": "NLYIQWLKDGGPSSGRPPPS",
      "1L3Y": "ECDTINCERYNGQVCGGPGRGLCFCGKCRCHPGFEGSACQA"
    }

    binding_predictions = predictor.predict_subsequences(protein_sequences, peptide_lengths=[9])

    # flatten binding predictions into a Pandas DataFrame
    df = binding_predictions.to_dataframe()

    # epitope collection is sorted by percentile rank
    # of binding predictions
    for binding_prediction in binding_predictions:
        print("Binder: %s" % (binding_prediction,))
        # if binding_prediction.affinity < 10000:

if __name__ == '__main__':
    run_netMHCpan()
