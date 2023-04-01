# tempus_tech_challenge

Bioinformatics technical challenge

## Instructions

Each variant in [`test_vcf_data.txt`](test_vcf_data.txt) must be annotated with the following pieces of information:

1. Depth of sequence coverage at the site of variation.
2. Number of reads supporting the variant.
3. Percentage of reads supporting the variant versus those supporting reference reads.
4. Using the [VEP hgvs API](https://rest.ensembl.org/#VEP), get the gene of the variant, type of variation (substitution, insertion, CNV, etc.) and their effect (missense, silent, intergenic, etc.).
5. The minor allele frequency of the variant if available.
6. Any additional annotations that you feel might be relevant.

## Solution

* Exploratory prototyping in [`prototype.ipynb`](prototype.ipynb)
* Utility functions in [`utils.py`](utils.py)
* I chose to employ the [PyVCF](https://pyvcf.readthedocs.io/en/latest/) package for easy parsing of VCF records
