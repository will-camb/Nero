import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file",
                    help="Should be named in form anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-phasefile",
                    help="Should be named in form chr#.merged.phase.gz",
                    required=True)
parser.add_argument("-idfile",
                    help="Directory of idfile used in chromopainter ('individual pop 1')",
                    required=True)

# Read copyprobs file, using dtype dictionary 
col_names = pd.read_csv(str(args.copyprobs_file), sep=" ", nrows=0).columns
types_dict = {'0': str}
types_dict.update({col: 'int8' for col in col_names if col not in types_dict})
anc_copyprobs = pd.read_csv(str(args.copyprobs_file), sep=" ", dtype=types_dict)
anc_copyprobs.set_index("0", inplace=True)

# Extract ancestry of interest
anc_copyprobs = anc_copyprobs.iloc[0::10, :] # Extract CHG only
# CHG = anc_copyprobs.iloc[0::10, :]
# EHG = anc_copyprobs.iloc[1::10, :]
# FarmerAnatolian = anc_copyprobs.iloc[2::10, :]
# FarmerEarly = anc_copyprobs.iloc[3::10, :]
# FarmerLate = anc_copyprobs.iloc[4::10, :]
# FarmerMiddle = anc_copyprobs.iloc[5::10, :]
# EastAsian = anc_copyprobs.iloc[6::10, :]
# African = anc_copyprobs.iloc[7::10, :]
# WHG = anc_copyprobs.iloc[8::10, :]
# Yamnaya = anc_copyprobs.iloc[9::10, :]

copyprobs_haps = list()
for h in range(int(anc_copyprobs.shape[0] / 2)):
    copyprobs_haps.extend([1, 2])
# Read idfile, with individuals in same order as phasefile
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
# Read phasefile, and add columns as positions and index as sample names from idfile
phase = pd.read_csv(str(args.phasefile), header=None, sep=" ", dtype='int8')
phase.columns = anc_copyprobs.columns.tolist()[::-1]  # Phasefile cols should be ascending
phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
phase_haps = list()
for h in range(int(phase.shape[0] / 2)):
    phase_haps.extend([1, 2])
# Option to extract a subset of SNPs; as default we take all
list_of_SNPs = phase.columns.tolist()
phase_temp = phase[list_of_SNPs]
phase_temp = phase_temp.reset_index().rename(columns={'index': 'ID'})
phase_temp['haps'] = phase_haps
anc_copyprobs_temp = anc_copyprobs[list_of_SNPs]
anc_copyprobs_temp = anc_copyprobs_temp.reset_index().rename(columns={'0': 'ID'})
anc_copyprobs_temp['haps'] = copyprobs_haps
# Merge copyprobs and phasefile as 3D dataframe
merged_phase_copyprobs = pd.merge(phase_temp, anc_copyprobs_temp, on=['ID', 'haps'])
# Re-order columns and rename etc
reorder = [[str(x) + "_x", str(x) + "_y"] for x in list_of_SNPs]
reorder = [item for sublist in reorder for item in sublist]
merged_phase_copyprobs = merged_phase_copyprobs[['ID', 'haps'] + reorder]
iterables = [["ID"] + list_of_SNPs, ["phase", "copyprobs"]]
merged_phase_copyprobs.columns = pd.MultiIndex.from_product(iterables, names=['first', 'second'])
merged_phase_copyprobs.set_index('ID', inplace=True)
# Now we have a 3D dataframe with phase and copyprobs values for each painted position, with the index as (individual#, haplotype#)
# We can extract the sum of the phasefile values at each position and place in a dictionary:
phase_sum = merged_phase_copyprobs.loc[:, (slice(None), 'phase')].sum().tolist()
for n, i in enumerate(list_of_SNPs):
    phase_sum_dict[i] = phase_sum[n]

# If we have a mapping of whether we want the 0 or 1 values from the phasefile, we can access the corresponeding values in the copyprobs column for a given SNP:
# Where i is the SNP position we want (as a string)
# phase_val is a 0 or 1 depending on which we want to select
merged_phase_copyprobs.loc[merged_phase_copyprobs[i]['phase'] == phase_val].loc[:,(i, 'copyprobs')]
