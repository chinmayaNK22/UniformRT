## Fetch column indexes for the columns of interest from Skyline output file
from itertools import islice


def get_dia_header_idx(infile):
    with open(infile) as file:
        for i in islice(file, 0, 1):
            split_i = i.rstrip().split('\t')
            pep_idx = split_i.index("Peptide")
            try:
                mod_pep_idx = split_i.index("Modified Sequence")
            except:
                mod_pep_idx = split_i.index("Peptide")
            pro_idx = split_i.index("Protein Name")
            qvalue_idx = split_i.index("Detection Q Value")
            z = split_i.index("Precursor Charge")
            raw_file = split_i.index("Replicate")
            miss_cleave_idx = split_i.index("Missed Cleavages")
            mz = split_i.index("Precursor Mz")
            rt = split_i.index("Peptide Retention Time")
            return pep_idx, pro_idx, qvalue_idx, z, miss_cleave_idx, mz, rt, raw_file, mod_pep_idx


if __name__== "__main__":
    get_dia_header_idx(args.infile[0])
