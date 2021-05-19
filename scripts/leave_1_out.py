import pandas as pd
import argparse
import seaborn as sns
import os
import sys

sns.set(color_codes=True)

parser = argparse.ArgumentParser()
parser.add_argument("-chr",
                    help="The chromosome that the sites of interest are on",
                    required=True)
parser.add_argument("-sites_file",
                    help="The positions of the sites for which average minor allele ancestry is to be calculated, "
                         "one per row",
                    required=True)
args = parser.parse_args()

sites_list = pd.read_csv(args.sites_file, header=None)[0].tolist()
sitesDf = pd.read_csv(args.sites_file, header=None)
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
phase = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/phasefiles_transformed/" + str(args.chr) +
                    "merged.phase.gz", header=None, sep=" ", dtype='int8')
cols = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/WHG." +
                   str(args.chr) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", nrows=0).columns.tolist()
del cols[0]
cols = [int(x) for x in cols]
cols_reversed = cols[::-1]

phase.columns = cols_reversed
phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]




wrong_right_map = pd.DataFrame(cols)
wrong_right_map.loc[:, 1] = cols_reversed
wrong_right_map.columns = ['Wrong', 'Right']
mapped_positions = pd.merge(sitesDf, wrong_right_map, left_on=0, right_on='Right')  # Check this
cols_mapped = mapped_positions['Wrong'].tolist()
types_dict = {'0': str}
types_dict.update({col: 'int8' for col in cols_mapped if col not in types_dict})

for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
    copyprobs = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/" + str(anc) + "." +
                            str(args.chr) + ".master_all_copyprobsperlocus.txt.gz",
                            sep=" ", dtype=types_dict, usecols=types_dict.keys())
    hapslist = list()
    for h in range(int(copyprobs.shape[0] / 2)):
        hapslist.extend([1, 2])
    copyprobs["haps"] = hapslist





    haps = pd.read_csv("/willerslev/ukbiobank/migraine_analysis/output_files/" + rsID + ".0.haps", header=None)
    merged = pd.merge(haps, copyprobs, left_on=[0, 1], right_on=['0', 'haps'])
    merged.set_index(0, inplace=True)
    merged.drop([1, '0', 'haps'], inplace=True, axis=1)
    df = pd.DataFrame(merged.mean().tolist())
    merged_cols = pd.DataFrame([int(col) for col in merged.columns.tolist()])
    pd.merge(wrong_right_map, merged_cols, left_on='Wrong', right_on=0)
    df[1] = pd.merge(wrong_right_map, merged_cols, left_on='Wrong', right_on=0)['Right']  # Saved in correct order
    df.to_csv("/willerslev/ukbiobank/migraine_analysis/ancestry_plots/" + rsID + "." + str(anc) + ".0.snp.subset.csv",
              header=False,
              index=False)
