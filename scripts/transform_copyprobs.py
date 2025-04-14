import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-chr",
                    help="the chromosome, NB must be named in form $chr.master_all_copyprobsperlocus.txt",
                    required=True)
parser.add_argument("-phasefiles_location",
                    help="the location of phasefiles, named in form filtered.$chr.phase",
                    required=True)
args = parser.parse_args()

positions = pd.read_csv(f"{args.phasefiles_location}/{args.chr}.merged.phase.gz", skiprows=2, nrows=1, sep=" ",
                        header=None).T.drop(0)
anc_copyprobs = pd.read_csv(f"{args.chr}.master_all_copyprobsperlocus.txt.gz", header=None, index_col=0)
print("Done reading files")
anc_copyprobs = anc_copyprobs[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
anc_copyprobs.columns = positions[0].tolist()[::-1]

# counts = pd.DataFrame(anc_copyprobs.index.value_counts())
# to_drop = counts.loc[counts[0] != 20].index.tolist()
# anc_copyprobs.drop(to_drop, inplace=True)

anc_copyprobs.to_csv(f"transformed_copyprobs/{args.chr}.master_all_copyprobsperlocus.txt.gz", compression='gzip')
print(
    f"Saving transformed file (space separated, positions added as header in correct (reverse) order) as "
    f"transformed_copyprobs/{args.chr}.master_all_copyprobsperlocus.txt.gz")