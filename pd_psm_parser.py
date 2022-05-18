## Fetch column ids for the columns of interest to be used for

def get_dda_header_idx(infile):
    for i in open(infile):
        split_i = i.rstrip().split('\t')
        if '"' in split_i[0]:
            try:                
                seq = split_i.index('"Annotated Sequence"')
                mod = split_i.index('"Modifications"')
                charge = split_i.index('"Charge"')
                rt = split_i.index('"RT [min]"')
                mz = split_i.index('"m/z [Da]"')
                pro = split_i.index('"Master Protein Accessions"')
                scan = split_i.index('"First Scan"')
                file = split_i.index('"Spectrum File"')
                #searchid = split_i.index('"Search ID"')
                return seq, mod, pro, charge, rt, mz, scan, file
                break
            except:
                seq = split_i.index('"Sequence"')
                mod = split_i.index('"Modifications"')
                charge = split_i.index('"Charge"')
                rt = split_i.index('"RT [min]"')
                pro = split_i.index('"Master Protein Accessions"')
                mz = split_i.index('"m/z [Da]"')
                scan = split_i.index('"First Scan"')
                file = split_i.index('"Spectrum File"')
                #searchid = split_i.index('"Search ID"')
                return seq, mod, pro, charge, rt, mz, scan, file
                break
        else:
            try:                
                seq = split_i.index('Annotated Sequence')
                mod = split_i.index('Modifications')
                charge = split_i.index('Charge')
                rt = split_i.index('RT [min]')
                mz = split_i.index('m/z [Da]')
                pro = split_i.index('Master Protein Accessions')
                scan = split_i.index('First Scan')
                file = split_i.index('Spectrum File')
                #searchid = split_i.index('Search ID')
                return seq, mod, pro, charge, rt, mz, scan, file
                break
            except:
                seq = split_i.index('Sequence')
                mod = split_i.index('Modifications')
                charge = split_i.index('Charge')
                rt = split_i.index('RT [min]')
                mz = split_i.index('m/z [Da]')
                pro = split_i.index('Master Protein Accessions')
                scan = split_i.index('First Scan')
                file = split_i.index('Spectrum File')
                #searchid = split_i.index('Search ID')
                return seq, mod, pro, charge, rt, mz, scan, file
                break

if __name__== "__main__":
    get_header_idx(args.infile[0])
