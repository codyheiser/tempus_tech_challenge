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

* Utility functions in [`utils.py`](utils.py)
* I chose to employ the [PyVCF](https://pyvcf.readthedocs.io/en/latest/) package for easy parsing of VCF records
* Exploratory prototyping and interactive annotation in [`prototype.ipynb`](prototype.ipynb)
* `cli.py` wraps the functions in a command-line interface:

```
usage: cli.py [-h] [-j N_JOBS] [-o [OUTFILE]] file

positional arguments:
  file                  VCF file to annotate

optional arguments:
  -h, --help            show this help message and exit
  -j N_JOBS, --n-jobs N_JOBS
                        Number of jobs for parallel. If 1, run serially.
  -o [OUTFILE], --outfile [OUTFILE]
                        Output file for annotated VCF. Default 'test_annotations.csv'.
```

### Output (columns in `test_annotations.csv`)

Result of `python cli.py test_vcf_data.txt -j 5 -o test_annotations.csv`:

1. `"NR"` - Depth of sequence coverage at the site of variation.
2. `"NV"` - Number of reads supporting the variant.
3. `"VAF"` - Percentage of reads supporting the variant versus `"RAF"` - those supporting reference reads.
4. Using the [VEP hgvs API](https://rest.ensembl.org/#VEP):
 * `"gene_id"` & `"gene_symbol"` - the gene of the variant
 * `"ALT_type"` - type of variation (substitution, insertion, CNV, etc.)
 * `"most_severe_consequence"` - variant effect (missense, silent, intergenic, etc.).
5. `"minor_allele_freq"` - The minor allele frequency of the variant if available.
6. Any additional annotations that you feel might be relevant:
 * `"snp_id"` and `"minor_allele"` for any SNPs (see **5**)
 * `"gene_biotype"` - type of gene the variant is found on
 * `"gene_description"` - description of gene the variant is found on
