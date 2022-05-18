# RetentionTime-Consistency-Calculator
LC-MS/MS derived peptide retention time deviation calculator across replicates for DDA and DIA derived result files.

## How to use RetentionTime_Consistency_Calculator in Windows/Linux
```
usage: RetentionTime_Consistency_Calculator.py [-h]
                                                -ex [-ex ...] -ip [-ip ...]
                                                -rt [-rt ...]

Variation from run-to-run retention time (RT) deviation in LC-MS derived DDA
or DIA results

positional arguments:
  -ex         Mention the type of acquisition method used (DDA/DIA)
  -ip         This tool can take Peptide Spectrum Match (PSM) table from
              Proteome Discoverer or DIA-Spectral Library search output from
              Skyline in tab delimitted format (txt/tsv)
  -rt         Minimum number of RTs to be considered for the analysis

optional arguments:
  -h, --help  show this help message and exit
```
