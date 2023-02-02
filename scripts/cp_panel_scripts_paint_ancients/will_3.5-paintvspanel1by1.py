import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-idfile",
                    help="file containing all IDs to process; only those in pop 'allAncients' will be processed in current implementation",
                    required=True)
parser.add_argument("-cp_panel_scripts",
                    help="the location of the panel scripts which will be copied to each temp directory",
                    required=True)
args = parser.parse_args()
idfile = pd.read_csv(args.idfile, sep=" ", header=None)
idfile = idfile[idfile[1]=='allAncients']
idfile.reset_index(inplace=True)

myfile = open('paintvspanel1by1_commands.txt', 'w')
for index, row in idfile.iterrows():
    myfile.write(f"bash paintsample1by1.sh {row[0]} {row['index']+1} {(row['index']+2)*2} {args.cp_panel_scripts} \n")
myfile.close()

