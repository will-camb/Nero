import pandas as pd 
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("-phasefile",
                    help="phasefile to extract sites from",
                    required=True)
parser.add_argument("-o",
                    help="location to save all_copyprobsperlocus.out to",
                    required=True)
args = parser.parse_args()

phaseDF = pd.read_csv(args.phasefile, delim_whitespace=True, nrows=3)

list_of_sites=phaseDF.iloc[2,:].tolist()
pd.DataFrame()