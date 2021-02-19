import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(color_codes=True)

parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file",
                    help="The copyprobs file to calculate ",
                    required=True)

args = parser.parse_args()


col_names = pd.read_csv(args.copyprobs_file, sep=" ", nrows=0).columns
types_dict = {'0': str}
types_dict.update({col: 'int8' for col in col_names if col not in types_dict})
copyprobs = pd.read_csv(args.copyprobs_file, sep=" ", dtype=types_dict)
copyprobs.set_index("0", inplace=True)
copyprobs.columns = copyprobs.columns.astype(int)
df = pd.DataFrame(copyprobs.mean().tolist())
df[1] = copyprobs.columns.tolist()
df.to_csv("average_ancestry_per_snp/average_ancestry_per_snp_" + args.copyprobs_file, header=False, index=False)
fig, ax = plt.subplots(figsize=[50, 30])
ax.set_title(str(args.copyprobs_file), {'fontsize': 40}, pad=10)
sns_plot = sns.scatterplot(x=df[1], y=df[0], s=30, ax=ax)
ax.set(ylim=(0, 8))
ax.set_ylabel('Average Painting Probability', size=40)
ax.set_xlabel('Position', size=40)
plt.axhline(copyprobs.mean().mean(), color='red', ls='--')
sns_plot.figure.savefig("average_ancestry_per_snp/average_ancestry_per_snp_" + args.copyprobs_file + ".png")
