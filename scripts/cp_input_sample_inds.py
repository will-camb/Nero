import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-pop_ids",
                    help="The pop_id file that you wan to sample from",
                    required=True)
parser.add_argument("-subsample_idfile",
                    help="A file containing the individuals to keep, one per line",
                    required=True)
parser.add_argument("-path",
                    help="Path to save new idfile to",
                    required=True)
args = parser.parse_args()
path = args.path
if not os.path.exists(path):
    os.makedirs(path)

#Load id files
subsample_idfile = pd.read_csv(args.subsample_idfile, sep=" ", header=None)
subsample_idfile[1] = "keep"
pop_idfile = pd.read_csv(args.pop_ids, sep=" ", header=None)

df = pop_idfile.reset_index().merge(subsample_idfile, how="left").set_index('index')
df2=df[df[1]!='keep']

#Make new idfile
with open(os.path.join(path, "pop_ids_subsampled"), "w") as file:
    for index, row in subsample_idfile.iterrows():
        file.write(str(row[0]) + " " + str(row[0]) + " 1\n")

# list of index value (in pop_ids) of individuals to remove
lines_to_delete = df2.index.to_list()

# for phase file
string = ''.join(str(2*(i+1)+2) + "d;" + str(2*(i+1)+3) + "d;" for i in lines_to_delete)[:-1]
print("To keep only individuals specified in subsample_idfile in phasefile: sed '" + string + "' phasefile > subsampled_phasefile")
print("NB you must alter the first line of subsampled_phasefile to reflect the new number of haplotypes (i.e. 2 x no individuals)")
