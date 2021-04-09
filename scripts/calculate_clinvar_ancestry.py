import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
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
    col_names = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/" + str(anc) + "." +
                            str(args.chr) + ".master_all_copyprobsperlocus.txt.gz",
                            sep=" ",
                            nrows=0).columns
    types_dict = {'0': str}
    types_dict.update({col: 'int8' for col in col_names if col not in types_dict})
    copyprobs = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/" + str(anc) + "." +
                            str(args.chr) + ".master_all_copyprobsperlocus.txt.gz",
                            sep=" ",
                            dtype=types_dict)
    hapslist = list()
    for h in range(int(copyprobs.shape[0] / 2)):
        hapslist.extend([1, 2])
    copyprobs["haps"] = hapslist
    for rsID in rsID_list:
        if not os.path.exists("/willerslev/ukbiobank/clinvar_analysis/output_files/" + rsID + ".haps"):
            continue
        if os.path.exists("/willerslev/ukbiobank/clinvar_analysis/ancestry_plots/" + rsID + "." + str(anc) + ".csv"):
            continue
        haps = pd.read_csv("/willerslev/ukbiobank/clinvar_analysis/output_files/" + rsID + ".haps", header=None)
        merged = pd.merge(haps, copyprobs, left_on=[0, 1], right_on=['0', 'haps'])
        merged.set_index("0", inplace=True)
        merged.drop([0, 1, 'haps'], inplace=True, axis=1)
        df = pd.DataFrame(merged.mean().tolist())
        df[1] = merged.columns.tolist()
        df.to_csv("/willerslev/ukbiobank/clinvar_analysis/ancestry_plots/" + rsID + "." + str(anc) + ".csv",
                  header=False,
                  index=False)
        # fig, ax = plt.subplots(figsize=[50, 30])
        # ax.set_title(rsID + " " + str(anc), {'fontsize': 40}, pad=10)
        # sns_plot = sns.scatterplot(x=df[1], y=df[0], s=30, ax=ax)
        # ax.set(ylim=(0, 8))
        # ax.set_ylabel('Average Painting Probability', size=40)
        # ax.set_xlabel('Position', size=40)
        # plt.axhline(merged.mean().mean(), color='red', ls='--')
        # sns_plot.figure.savefig("/willerslev/ukbiobank/clinvar_analysis/ancestry_plots/" + rsID + "." + str(anc) + ".png")
