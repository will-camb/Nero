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
args = parser.parse_args()

mapping = pd.read_csv("/willerslev/ukbiobank/clinvar_analysis/rsID-chr_mapping", header=None, sep="\t", dtype=str)
rsID_list = mapping[mapping[2] == str(args.chr)][0].tolist()
rsID_files = list()
for rsID in rsID_list:
    rsID_files.append("/willerslev/ukbiobank/clinvar_analysis/output_files/" + str(rsID) + ".haps")
rsID_exists = [f for f in rsID_files if os.path.isfile(f)]
if not rsID_exists:
    print("No rsIDs with valid input for chr " + str(args.chr))
    sys.exit(0)
for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
    cols = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/" + str(anc) + "." +
                       str(args.chr) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", nrows=0).columns.tolist()
    del cols[0]
    cols = [int(x) for x in cols]
    cols_reversed = cols[::-1]
    wrong_right_map = pd.DataFrame(cols)
    wrong_right_map.loc[:, 1] = cols_reversed
    wrong_right_map.columns = ['Wrong', 'Right']
    for rsID in rsID_list:
        if not os.path.exists("/willerslev/ukbiobank/clinvar_analysis/output_files/" + rsID + ".haps"):
            continue
        # if os.path.exists("/willerslev/ukbiobank/clinvar_analysis/ancestry_plots/" + rsID + "." + str(anc) + ".snp"
        #                                                                                                      ".subset"
        #                                                                                                      ".csv"):
        #     continue
        print("Processing " + str(rsID))
        position = pd.read_csv("/willerslev/ukbiobank/clinvar_analysis/output_files/" + rsID + ".position",
                               skiprows=1, nrows=1, sep="\t").loc[:, 'position'][0]
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
        hapslist = list()
        for h in range(int(copyprobs.shape[0] / 2)):
            hapslist.extend([1, 2])
        copyprobs["haps"] = hapslist
        haps = pd.read_csv("/willerslev/ukbiobank/clinvar_analysis/output_files/" + rsID + ".haps", header=None)
        merged = pd.merge(haps, copyprobs, left_on=[0, 1], right_on=['0', 'haps'])
        merged.set_index(0, inplace=True)
        merged.drop([1, '0', 'haps'], inplace=True, axis=1)
        df = pd.DataFrame(merged.mean().tolist())
        merged_cols = pd.DataFrame([int(col) for col in merged.columns.tolist()])
        pd.merge(wrong_right_map, merged_cols, left_on='Wrong', right_on=0)
        df[1] = pd.merge(wrong_right_map, merged_cols, left_on='Wrong', right_on=0)['Right']  # Saved in correct order
        df.to_csv("/willerslev/ukbiobank/clinvar_analysis/ancestry_plots/" + rsID + "." + str(anc) + ".snp.subset.csv",
                  header=False,
                  index=False)
