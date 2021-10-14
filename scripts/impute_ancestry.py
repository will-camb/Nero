import pandas as pd
import argparse
from pathlib import Path
pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file",
                    help="Should be named in form anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-snp_file",
                    help="File containing details of SNPs to be imputed. Should be two columns: chr and position",
                    required=True)
parser.add_argument("-chr",
                    help="Chromosome number",
                    required=True)
parser.add_argument("-anc",
                    help="Ancestry being analysed",
                    required=True)
args = parser.parse_args()


def impute(impute_snp_list_):
    global painted_snp_list
    df = pd.DataFrame(index=anc_copyprobs.index)
    for snp in impute_snp_list_:
        if snp in painted_snp_list:
            df.loc[:, snp] = anc_copyprobs[str(snp)]
        else:
            painted_snp_list.append(snp)
            painted_snp_list.sort(reverse=True)
            temp_anc_copyprobs = anc_copyprobs.iloc[:, [painted_snp_list.index(snp) - 1, painted_snp_list.index(snp)]]
            D1 = abs(snp - temp_anc_copyprobs.columns.astype(int).tolist()[0])
            D2 = abs(snp - temp_anc_copyprobs.columns.astype(int).tolist()[1])
            temp_anc_copyprobs.loc[:, snp] = (1 - (D1 / (D1 + D2))) * temp_anc_copyprobs.iloc[:, 0] + (
                    1 - (D2 / (D1 + D2))) * temp_anc_copyprobs.iloc[:, 1]
            df.loc[:, snp] = temp_anc_copyprobs[snp]
    return df


Path("imputed/").mkdir(parents=True, exist_ok=True)
col_names = pd.read_csv(str(args.copyprobs_file), sep=" ", nrows=0).columns
types_dict = {'0': str}
types_dict.update({col: 'int8' for col in col_names if col not in types_dict})
anc_copyprobs = pd.read_csv(str(args.copyprobs_file), sep=" ", dtype=types_dict)
anc_copyprobs.set_index("0", inplace=True)
anc_copyprobs.columns = anc_copyprobs.columns.tolist()[::-1]  # Reverse column names because they are mis-labelled;
#  so when correct, column names are descending
painted_snp_list = anc_copyprobs.columns.astype(int).tolist()
impute_snps = pd.read_csv(args.snp_file)
impute_snps = impute_snps.loc[impute_snps.iloc[:, 0] == int(args.chr)]
impute_snp_list = impute_snps.iloc[:, 1].tolist()
impute_snp_list.sort(reverse=True)
imputed_snps = impute(impute_snp_list)
imputed_snps.columns = imputed_snps.columns.tolist()[::-1]  # Reverse column names to incorrect order for other scripts
imputed_snps.to_csv("imputed/temp." + str(args.anc) + "." + str(args.chr) + ".master_all_copyprobsperlocus.txt.gz",
                    compression='gzip')
