import argparse

from joblib import Parallel, delayed
from utils import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file",
        type=str,
        help="VCF file to annotate",
    )
    parser.add_argument(
        "-j",
        "--n-jobs",
        type=int,
        help="Number of jobs for parallel. If 1, run serially.",
        default=-1,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        help="Output file for annotated VCF. Default 'test_annotations.csv'.",
        nargs="?",
        default="test_annotations.csv",
    )
    args = parser.parse_args()

    vcf_read = vcf.Reader(open(args.file, "r"))  # open file
    if args.n_jobs == 1:  # serial processing
        out = []
        for i, record in enumerate(vcf_read):
            out.append(process_record(record, i))
    else:  # parallel with joblib
        out = Parallel(n_jobs=args.n_jobs)(
            delayed(process_record)(record, i) for i, record in enumerate(vcf_read)
        )

    out_df = pd.concat(out)  # concatenate outputs into master df

    calc_frequencies(out_df)  # calculate allele frequencies

    out_df.to_csv(args.outfile)  # save to outfile

    print("\nDone! Total records: {}".format(i))
