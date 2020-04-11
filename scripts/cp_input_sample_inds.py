import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-subsample_idfile",
                    help="A file containing the individuals to keep, one per line",
                    required=True
parser.add_argument("-phasefile",
                    help="The phasefile that you want to sample from",
                    required=True)
parser.add_argument("-recomb_file",
                    help="The recombfile that you want to sample from",
                    required=True)
parser.add_argument("-pop_ids",
                    help="The pop_id file that you wan to sample from",
                    required=True)
parser.add_argument("-cp_out",
                    help="Path to save sub-sampled output files to e.g. Documents/msprime2cptest",
                    required=True)
parser.add_argument("-seed",
                    help="If this is None, automatically generated",
                    required=False)
args = parser.parse_args()

subsample_idfile = pd.read_csv(args.subsample_idfile, sep=" ", header=None)
subsample_idfile[1] = "keep"
pop_idfile = pd.read_csv(args.pop_ids, sep=" ", header=None)
df = pd.merge(subsample_idfile, pop_idfile, how='outer')