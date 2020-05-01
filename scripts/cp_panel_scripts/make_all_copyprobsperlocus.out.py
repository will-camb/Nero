import pandas as pd 
import argparse 
import os

parser = argparse.ArgumentParser()
parser.add_argument("-phasefile",
                    help="phasefile to extract sites from",
                    required=True)
parser.add_argument("-chr",
                    help="which chromosome",
                    required=True)
parser.add_argument("-o",
                    help="location to save all_copyprobsperlocus.out to",
                    required=True)
args = parser.parse_args()
phaseDF = pd.read_csv(args.phasefile, delim_whitespace=True, nrows=3)
all_copyprobsperlocusDF = pd.DataFrame(phaseDF.iloc[2, :].tolist()).set_index(0, drop=True)
all_copyprobsperlocusDF.to_csv(os.path.join(args.o,'/', args.chr, '_all_copyprobsperlocus.txt'), sep=' ', index_label='pos')
