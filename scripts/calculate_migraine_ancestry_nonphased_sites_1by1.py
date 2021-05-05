import pandas as pd
import argparse
import seaborn as sns
import os
import sys

sns.set(color_codes=True)

parser = argparse.ArgumentParser()
parser.add_argument("-chr",
                    help="The chromosome",
                    required=True)
parser.add_argument("-rsID",
                    help="The chromosome",
                    required=True)
args = parser.parse_args()


if not os.path.isfile("output_files/" + str(args.rsID) + ".hom.1.samples"):
    print("No valid samples file in directory output_files for " + str(args.rsID))
    sys.exit(0)

for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
    print("Processing " + str(anc) + " for " + str(args.rsID))
    if os.path.isfile("ancestry_plots/" + args.rsID + "." + str(anc) + ".1.snp.subset.csv"):
        print("Already have ancestry plots for " + str(anc) + " for " + str(args.rsID))
        continue
    cols = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/" + str(anc) + "." +
                       str(args.chr) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", nrows=0).columns.tolist()
    del cols[0]
    cols = [int(x) for x in cols]
    cols_reversed = cols[::-1]
    wrong_right_map = pd.DataFrame(cols)
    wrong_right_map.loc[:, 1] = cols_reversed
    wrong_right_map.columns = ['Wrong', 'Right']
    position_mapping = pd.read_csv("rsID-position_mapping", sep="\t", header=None)
    position = position_mapping.loc[position_mapping[0] == args.rsID][1].item()
    closest_position = min(cols, key=lambda x: abs(x - position))  # Get closest position in cols, because not
    # all SNPs were painted
    mapped_position = cols[cols_reversed.index(closest_position)]
    max_position = int(mapped_position) + 6000000
    min_position = int(mapped_position) - 6000000
    cols_subset = [str(col) for col in cols if col in range(min_position, max_position)]
    print("Calculating averages for " + str(len(cols_subset)) + " SNPs")
    types_dict = {'0': str}
    types_dict.update({col: 'int8' for col in cols_subset if col not in types_dict})
    copyprobs = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/" + str(anc) + "." +
                            str(args.chr) + ".master_all_copyprobsperlocus.txt.gz",
                            sep=" ", dtype=types_dict, usecols=types_dict.keys())
    samples = pd.read_csv("output_files/" + str(args.rsID) + ".hom.1.samples", header=None)
    merged = pd.merge(samples, copyprobs, left_on=0, right_on='0')
    merged.set_index(0, inplace=True)
    merged.drop('0', inplace=True, axis=1)
    df = pd.DataFrame(merged.mean().tolist())
    merged_cols = pd.DataFrame([int(col) for col in merged.columns.tolist()])
    pd.merge(wrong_right_map, merged_cols, left_on='Wrong', right_on=0)
    df[1] = pd.merge(wrong_right_map, merged_cols, left_on='Wrong', right_on=0)['Right']  # Saved in correct order
    df.to_csv("ancestry_plots/" + str(args.rsID) + "." + str(anc) + ".1.snp.subset.csv",
              header=False,
              index=False)
