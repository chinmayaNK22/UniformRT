from itertools import islice
from statistics import variance, stdev, mean
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from pd_psm_parser import get_dda_header_idx
from skyline_parser import get_dia_header_idx

parser = argparse.ArgumentParser(description='''Variation from run-to-run retention time (RT) deviation in LC-MS derived DDA or DIA results''')

parser.add_argument('experiment', metavar='-ex', type=str, nargs='+', help='Mention the type of acquisition method used (DDA/DIA)')

parser.add_argument('infile', metavar='-ip', type=str, nargs='+', help='This tool can take Peptide Spectrum Match (PSM) table from Proteome Discoverer or DIA-Spectral Library search output from Skyline in tab delimitted format (txt/tsv)')

parser.add_argument('min_rts', metavar='-rt', type=str, nargs='+', help='Minimum number of RTs to be considered for the analysis')

args = parser.parse_args()


def DIA_RT_consistency(infile):
    dicts = {}
    a = get_dia_header_idx(infile)
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            pep = split_i[a[0]]
            mod_pep = split_i[a[-1]]
            pro = split_i[a[1]]
            charge = split_i[a[3]]
            condition = split_i[a[-2]]
            rt = split_i[a[-3]]
            fdr = split_i[a[2]]
            if fdr < '0.01':
                if condition not in dicts:
                    dicts[condition] = [pep + '@' + mod_pep + '@' + charge + '@' + rt]
                    
                else:
                    dicts[condition].append(pep + '@' + mod_pep + '@' + charge + '@' + rt)
                    
            else:
                pass
                #print (pep, pro, fdr, condition, rt)
    return dicts

def DDA_RT_consistency(infile):
    dicts = {}
    a = get_dda_header_idx(infile)
    with open(infile) as file:
        for i in islice(file, 1, None):
            split_i = i.rstrip().split('\t')
            pep = ""
            pro = ""
            rt = ""
            file = ""
            condition = []
            if '"' in split_i[4]:
                mod = split_i[a[1]].strip('"')
                pep = split_i[a[0]].strip('"').split('.')[1].upper()
                z = split_i[a[3]].strip('"')
                rt = split_i[a[4]].strip('"')
                condition = split_i[a[-1]].strip('"')
            else:
                mod = split_i[a[1]]
                pep = split_i[a[0]].split('.')[1].upper()
                rt = split_i[a[4]]
                z = split_i[a[3]]
                condition = split_i[a[-1]]
            if condition not in dicts:
                dicts[condition] = [pep + '@' + mod + '@' + z + '@' + rt]
            else:
                dicts[condition].append(pep + '@' + mod + '@' + z + '@' + rt)

    return dicts

def rt_variation_check(experiment, infile, min_rts):
    if experiment == 'DIA':
        file_pep_rt = DIA_RT_consistency(infile)
    elif experiment == 'DDA':
        file_pep_rt = DDA_RT_consistency(infile)

    output = []
    final = {}
    for k, v in file_pep_rt.items():
        precursors = {}
        for idx, i in enumerate(v):
            if i.split('@')[0] + '_' + i.split('@')[1] + '_' + i.split('@')[2] not in precursors:
                precursors[i.split('@')[0] + '_' + i.split('@')[1] + '_' + i.split('@')[2]] = [i.split('@')[-1]]
            else:
                precursors[i.split('@')[0] + '_' + i.split('@')[1] + '_' + i.split('@')[2]].append(i.split('@')[-1])

        for l, m in precursors.items():
            if '_'.join(k.split('_')[0:2]) not in final:
                final['_'.join(k.split('_')[0:2])] = [l + '_' + m[0]]
            else:
                final['_'.join(k.split('_')[0:2])].append(l + '_' + m[0])

    cv_count = 0
    cnt = 0
    dev_count = 0
    for k, v in final.items():
        precursors = {}
        for idx, i in enumerate(v):
            if '_'.join(i.split('_')[0:3]) not in precursors:
                precursors['_'.join(i.split('_')[0:3])] = [float(i.split('_')[-1])]
            else:
                precursors['_'.join(i.split('_')[0:3])].append(float(i.split('_')[-1]))

        for l, m in precursors.items():
            if len(m) >= int(min_rts):
                var = variance(m)
                std_dev = stdev(m)
                output.append([k, l] + [str(rts) for rts in m] + [str(mean(m)), str(std_dev), str(var)])
                cnt += 1
                if var < 2.0:
                    cv_count += 1
                if std_dev < 0.5:
                    dev_count += 1

    print ("There are " + str(cv_count) + " peptide precursors with Retention Time (RT) %CV < 2% out of " + str(cnt) + ".")
    print ("There are " + str(dev_count) + " peptide precursors with Retention Time (RT) deviation < 0.5 min (30 sec) out of " + str(cnt) + ".")
    
    runs = []
    for j in range(0, int(min_rts)):
        runs.append("RUN "+ str(j+1))

    header = ['Condition', 'Peptide Precursor'] + runs + ['Average RT (min)', 'RT Std. Dev. (min)', 'RT %CV']
    outfile = ""
    if infile.split('.')[-1] == 'txt':
        outfile = "{0}_RT_reproducibility_output.txt".format(infile.rstrip('txt').rstrip('.'))
    elif infile.split('.')[-1] == 'tsv':
        outfile = "{0}_RT_reproducibility_output.txt".format(infile.rstrip('tsv').rstrip('.'))

    with open(outfile, 'w') as outf:
        outf.write('\t'.join(header) + '\n')
        outf.writelines('\t'.join(i) + '\n' for i in output)

if __name__== "__main__":
    if int(args.min_rts[0]) >= 2:
        rt_variation_check(args.experiment[0], args.infile[0], args.min_rts[0])
    else:
        raise ValueError("The statistics cant be applied with " + args.min_rts[0] +  ", minimum of 2 replicate RT values should be used.")
