import pandas as pd
import argparse
from pathlib import Path
import sys

pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file",
                    help="Should be named in form anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-bp",
                    required=True)
parser.add_argument("-chr",
                    help="Chromosome number",
                    required=True)
parser.add_argument("-anc",
                    help="Ancestry being analysed",
                    required=True)
parser.add_argument("-reverse_cols",
                    help="Whether to reverse column names; can take two values: True or False",
                    required=True)
args = parser.parse_args()
if args.reverse_cols == 'True':
    reverse_cols = True
else:
    reverse_cols = False

if reverse_cols:  # Can't run on incorrectly ordered cols
    sys.exit()
else:  # For correctly labelled cols in copyprobs:
    def impute(impute_snp_list_):
        global painted_snp_list
        anc_copyprobs = pd.read_csv(args.copyprobs_file, sep=" ", usecols=['ID']).set_index('ID', drop=True)
        df = pd.DataFrame(index=anc_copyprobs.index)
        for snp in impute_snp_list_:
            if snp in painted_snp_list:
                temp_anc_copyprobs = pd.read_csv(args.copyprobs_file, sep=" ",
                                                 usecols=['ID', str(snp)],
                                                 dtype={'ID': str, str(snp): 'int8'}).set_index("ID", drop=True)
                df.loc[:, snp] = temp_anc_copyprobs[str(snp)]
            else:
                painted_snp_list.append(snp)
                painted_snp_list.sort(reverse=True)
                snp1 = painted_snp_list[painted_snp_list.index(snp)-1]
                snp2 = painted_snp_list[painted_snp_list.index(snp)+1]
                painted_snp_list.remove(snp)
                temp_anc_copyprobs = pd.read_csv(args.copyprobs_file, sep=" ",
                                                 usecols=['ID', str(snp1), str(snp2)],
                                                 dtype={'ID': str, str(snp1): 'int8', str(snp2): 'int8'})\
                    .set_index("ID", drop=True)
                D1 = abs(snp - snp1)
                D2 = abs(snp - snp2)
                temp_anc_copyprobs.loc[:, snp] = (1 - (D1 / (D1 + D2))) * temp_anc_copyprobs.iloc[:, 0] + (
                        1 - (D2 / (D1 + D2))) * temp_anc_copyprobs.iloc[:, 1]
                df.loc[:, snp] = temp_anc_copyprobs[snp]
        return df


    print("You have indicated that the columns in copyprobs are correctly labelled")
    Path("imputed/").mkdir(parents=True, exist_ok=True)
    if Path("imputed/temp." + str(args.anc) + "." + str(args.chr) + ".master_all_copyprobsperlocus.txt.gz").is_file():
        print("temp." + str(args.anc) + "." + str(args.chr) + ".master_all_copyprobsperlocus.txt.gz already imputed")
        sys.exit()
    painted_snp_list = pd.read_csv(str(args.copyprobs_file), sep=" ", nrows=0).columns.tolist()  # descending=correct
    del painted_snp_list[0]
    painted_snp_list = [int(x) for x in painted_snp_list]
    impute_snp_list = [args.bp]
    impute_snp_list = [int(x) for x in impute_snp_list]
    if len(impute_snp_list) == 0:
        print("No SNPs to be imputed for " + str(args.chr) + ", exiting")
        sys.exit()
    impute_snp_list.sort(reverse=True)  # descending
    imputed_snps = impute(impute_snp_list)  # NB saves in correct order!
    imputed_snps.to_csv("imputed/temp." + str(args.bp) + "." + str(args.anc) + ".master_all_copyprobsperlocus.txt.gz",
                        compression='gzip')
