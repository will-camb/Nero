# This script need adjusting for file names, locations, and parameters. It is intended as a guide for converting
# sparsepainter results into input files to run sourcefind to infer genome-wide ancestry proportions.
import pandas as pd


# Assemble target chunklength file (summed across chromosomes and both haps)
df_temp = pd.read_csv(f"output_files/SparsePainter.22_chunklength.txt.gz", sep=" ", index_col='indnames')
df_temp.loc[:] = 0
for chrom in range(1, 23):
    try:
        temp = pd.read_csv(f"output_files/SparsePainter.{chrom}_chunklength.txt.gz", sep=" ", index_col='indnames')
        df_temp = df_temp + temp
    except FileNotFoundError:
        print(f"No file found for chromosome {chrom}")
        continue
df_temp.reset_index(inplace=True)
df_temp['Recipient'] = df_temp['indnames'].str.rsplit('_', 1).str[0]
df_temp = df_temp.groupby("Recipient").sum()
df_temp.to_csv(f"output_files/SparsePainter_ALL_chunklength.txt", sep=" ")

# Assemble ref chunklength file (summed across chromosomes and both haps)
df_temp = pd.read_csv(f"output_files/SparsePainter.22_refvsref_chunklength.txt.gz", sep=" ", index_col='indnames')
df_temp.loc[:] = 0
for chrom in range(1, 23):
    try:
        temp = pd.read_csv(f'output_files/SparsePainter.{chrom}_refvsref_chunklength.txt.gz', sep=" ",
                           index_col='indnames')
        df_temp = df_temp + temp
    except FileNotFoundError:
        print(f"No file found for chromosome {chrom}")
        continue
df_temp.reset_index(inplace=True)
df_temp['Recipient'] = df_temp['indnames'].str.rsplit('_', 1).str[0]
df_temp = df_temp.groupby("Recipient").sum()
df_temp.to_csv(f"output_files/SparsePainter_refvsref_ALL_chunklength.txt", sep=" ")

# Make Sourcefind copyvectors file
reference = pd.read_csv(f"output_files/SparsePainter_refvsref_ALL_chunklength.txt", sep=" ")
target = pd.read_csv(f"output_files/SparsePainter_ALL_chunklength.txt", sep=" ")
copyvectors = reference.append(target)
copyvectors.to_csv(f"input.file.copyvectors.txt", index=False, sep=" ")

# ref_samples_popfile is the -popfile input to sparsepainter ("Population file of reference individuals that contains
# two columns without headers. The first column is the names of all the reference samples (must be in the same order as
# reffile). The second column is the population labels of the reference samples, which can be either strings or
# numbers"):
ref_samples_popfile = pd.read_csv(f"ref_samples_popfile.txt", header=None, sep=" ")
ref_samples_popfile[2] = 1
ref_samples_popfile.columns = ['Recipient', 'Population', 'Include']
target['Population'] = target['Recipient']
target['Include'] = 1
ref_samples_popfile.append(target)[['Recipient', 'Population', 'Include']].to_csv(f"input.file.ids.txt", header=False,
                                                                                  index=False, sep=" ")
with open('SourcefindParamfile.txt', 'w') as f:
    print("self.copy.ind: 0", file=f)
    print("num.surrogates: 5", file=f)
    print("exp.num.surrogates: 3", file=f)
    print("input.file.ids: input.file.ids.txt", file=f)
    print("input.file.copyvectors: input.file.copyvectors.txt", file=f)
    print("save.file: output", file=f)
    print("copyvector.popnames: " + ' '.join(ref_samples_popfile['Population'].unique()), file=f)
    print("surrogate.popnames: " + ' '.join(ref_samples_popfile['Population'].unique()), file=f)
    print('target.popnames: ' + ' '.join(target.Recipient), file=f)
    #     print('target.popnames: ' + ' '.join(ref_samples_popfile.Recipient), file=f) # Alternative to
    #     infer ref samples
    print('num.slots: 100', file=f)
    print('num.iterations: 2000', file=f)  # Quick:2000, Accurate:20000
    print('num.burnin: 500', file=f)  # Quick:500, Accurate:5000
    print('num.thin: 50', file=f)  # Quick:50, Accurate:500

with open('run_model.sh', 'w') as f:
    print('# !/bin/bash', file=f)
    print('copyvectors = input.file.copyvectors.txt', file=f)
    print('idfile = input.file.ids.txt', file=f)
    print('paramfile = IA.SourcefindParamfile.txt', file=f)
    print('SOURCEFIND = $PATHTO/sourcefindv2.R', file=f)
    print('Rscript ${SOURCEFIND} --chunklengths ${copyvectors} --parameters ${paramfile} --idfile ${idfile} --output '
          'output', file=f)
# Now you can run sourcefind with a command like: bash run_model.sh
