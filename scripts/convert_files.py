import pandas as pd
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-directory",
                    required=True)
parser.add_argument("-chr1",
                    required=True)
parser.add_argument("-chr2",
                    required=True)
args = parser.parse_args()


for i in range(int(args.chr1), int(args.chr2)):
    my_file = Path(str(args.directory) + "/temp.Farmer."+str(i)+".master_all_copyprobsperlocus.txt.gz")
    positions = pd.read_csv(str(args.directory) + "/phasefiles/" + str(i) + ".merged.phase",
                            skiprows=2,
                            nrows=1,
                            sep=" ",
                            header=None).T.drop(0)
    if not my_file.is_file():
        print("Reading Farmer chr" + str(i) + " in " + str(args.directory))
        copyprobsFarmer = pd.read_csv(str(args.directory) + "/temp.FarmerAnatolian."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              pd.read_csv(str(args.directory) + "/temp.FarmerEarly."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              pd.read_csv(str(args.directory) + "/temp.FarmerLate."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                              pd.read_csv(str(args.directory) + "/temp.FarmerMiddle."+str(i)+".master_all_copyprobsperlocus.txt.gz",
                              header=None,
                              index_col=0)[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
        copyprobsFarmer.columns = positions[0].tolist()
        copyprobsFarmer.to_csv(str(args.directory) + "/temp.Farmer."+str(i)+".master_all_copyprobsperlocus.txt", sep=" ")
        print("Farmer done, now converting other ancestries to correct format")
    for anc in ["CHG", "EHG", "African", "EastAsian", "WHG", "Yamnaya"]:
        print("Processing " + str(anc) + " chromosome " + str(i) + " in " + str(args.directory))
        anc_file = Path(str(args.directory) + "/temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt")
        if not anc_file.is_file():
            anc_copyprobs = pd.read_csv(str(args.directory) + "/temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt.gz", header=None, index_col=0)
            anc_copyprobs = anc_copyprobs[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
            anc_copyprobs.columns = positions[0].tolist()
            anc_copyprobs.to_csv(str(args.directory) + "/temp." + str(anc) + "." + str(i) + ".master_all_copyprobsperlocus.txt", sep=" ")
            print(str(anc) + " done!")
