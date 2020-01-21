import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input", help="Path and filename of cp copyprobsperlocus.out.csv file", required=True)
parser.add_argument("-out", help="Path to save cp best guess csv to", required=True)
args = parser.parse_args()

df = pd.read_csv(args.input, delim_whitespace=True)

rowstoignore = list(df[df['pos'].str.contains("HAP") == True].index)
introws = list()
for i in list(df.index):
    if i not in rowstoignore:
        introws.append(i)

labellist = list()
for i in range(2 * (len(df.columns) - 1)):
    for j in range(239546):
        labellist.append("HAP_" + str(i))

df['hap_id'] = labellist
df_reduced = df.iloc[introws]
df_reduced.iloc[:, 1:13] = df_reduced.iloc[:, 1:13].astype(float)
df_reduced['ancestry'] = df_reduced.iloc[:, 1:13].idxmax(axis=1)

df_final_ancestry = df_reduced[['hap_id', 'pos', 'ancestry']]
df_final_ancestry.to_csv(args.out + '/df_final_ancestry.csv')
