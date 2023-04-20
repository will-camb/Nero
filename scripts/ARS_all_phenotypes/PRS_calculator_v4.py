#  This script is for already pruned SNPs
#  Only runs one phenotype at a time
#  Looks at homozygous individuals
#  Requirements: imputed painting scores for pruned SNPs, homozygous individuals for each variant

import pandas as pd
import argparse
import os
import numpy as np
import sys
pd.options.mode.chained_assignment = None


def analyse_anc(anc_copyprobs_temp_, anc, chrom, pval, iteration, list_of_SNPs_, GWAS_variants_pval_):
    skipped_snps = 0
    alt_mapping = {}
    output = pd.DataFrame(columns=['pos', 'alt', 'ref', 'alt_average', 'ref_average', 'OR', 'rsID'])
    for n, i in enumerate(list_of_SNPs_):
        i = int(i)
        try:
            GWAS_variants_pval_.reset_index(drop=True, inplace=True)
            EUR_maf = GWAS_variants_pval_.loc[GWAS_variants_pval_['BP'] == i].effect_allele_frequency.item()
            OR = GWAS_variants_pval_.loc[GWAS_variants_pval_['BP'] == i].OR.item()
            rsid = GWAS_variants_pval_.loc[GWAS_variants_pval_['BP'] == i].SNP.item()
            samples_0 = pd.read_csv(args.output_files + str(rsid) + ".hom.0.samples", header=None)
            samples_1 = pd.read_csv(args.output_files + str(rsid) + ".hom.1.samples", header=None)
        except Exception as e:
            skipped_snps += 1
            print("Couldn't find rsID/EUR_maf/output file for position:" + str(i) + ", so skipping: " + str(e))
            continue
        try:
            anc_copyprobs_temp_i = anc_copyprobs_temp_[['ID', str(i)]]
        except Exception as e:
            print("Couldn't find SNP with pos " + str(i) + " in anc_copyprobs: " + str(e))
            continue
        if EUR_maf == -1:
            if samples_0.shape[0] < (len(copyprobs_haps) / 8):
                alt_mapping[i] = 0
            else:
                alt_mapping[i] = 1
        elif EUR_maf <= 0.5:
            #  minor is effect
            if samples_0.shape[0] < (len(copyprobs_haps) / 8):
                alt_mapping[i] = 0
            else:
                alt_mapping[i] = 1
        elif EUR_maf > 0.5:
            if samples_0.shape[0] < (len(copyprobs_haps) / 8):
                alt_mapping[i] = 1
            else:
                alt_mapping[i] = 0
        anc_copyprobs_temp_i[str(i)] = anc_copyprobs_temp_i[str(i)].astype('float64')
        if alt_mapping[i] == 0:
            alt_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_0[0].tolist())][str(i)].sum()
            ref_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_1[0].tolist())][str(i)].sum()
            alt_average = alt_sum / len(samples_0[0].tolist())
            ref_average = ref_sum / len(samples_1[0].tolist())
        elif alt_mapping[i] == 1:
            alt_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_1[0].tolist())][str(i)].sum()
            ref_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_0[0].tolist())][str(i)].sum()
            alt_average = alt_sum / len(samples_1[0].tolist())
            ref_average = ref_sum / len(samples_0[0].tolist())
        output.loc[n] = [i, alt_sum, ref_sum, alt_average, ref_average, OR, rsid]
    output['maf'] = output['alt'] / (output['alt'] + output['ref'])
    output['maf'] = output['maf'].fillna(0)
    output['maf_average'] = output['alt_average'] / (output['alt_average'] + output['ref_average'])
    output['maf_average'] = output['maf_average'].fillna(0)
    output['beta'] = np.log(output['OR'].astype(float))  # Transform OR to beta
    output['maf_x_beta'] = output['maf'] * output['beta']
    PRS = output['maf_x_beta'].sum()
    maf_GWAS = output
    number_of_SNPs = len(list_of_SNPs_)
    results_list.append(
        [args.copyprobs_file_imputed, "phenotype", anc, chrom, PRS, number_of_SNPs, skipped_snps, pval, iteration,
         maf_GWAS['maf'].tolist(), maf_GWAS['maf_average'].tolist(), maf_GWAS['beta'].tolist(),
         maf_GWAS['pos'].tolist(), maf_GWAS['rsID'].tolist()])


parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file_imputed",
                    help="Should be named in form temp.anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-output_files",
                    help="Location of output files from get_nonphased_samples.sh",
                    required=True)
parser.add_argument("-variants",
                    help="GWAS variants file containing SNP details - "
                         "SNP CHR BP other_allele effect_allele effect_allele_frequency P OR",
                    required=True)
parser.add_argument("-chr",
                    help="Chromosome number",
                    required=True)
parser.add_argument("-anc",
                    help="Ancestry being analysed",
                    required=True)
parser.add_argument("-bootstrap",
                    help="Whether to bootstrap or not; can take two values: True or False",
                    required=True)
parser.add_argument("-reverse_cols",
                    help="Whether to reverse column names; can take two values: True or False. NB this refers to "
                         "copyprobs_file_original AND copyprobs_file_imputed",
                    required=True)
args = parser.parse_args()
if args.bootstrap == 'True':
    bootstrap = True
else:
    bootstrap = False
if args.reverse_cols == 'True':
    reverse_cols = True
else:
    reverse_cols = False

# Read in copyprobs
if not os.path.exists(str(args.copyprobs_file_imputed)):
    print("***Can't find copyprobs_file_imputed for " + args.anc + " chr" + str(args.chr) + "***")
    sys.exit()
if reverse_cols:  # For copyprobs files with incorrectly ordered columns
    print("You have indicated cols in copyprobs files are incorrectly ordered")
    col_names = pd.read_csv(str(args.copyprobs_file_imputed), sep=" ", nrows=0).columns
    types_dict = {'0': str}
    types_dict.update({col: 'float16' for col in col_names if col not in types_dict})
    anc_copyprobs = pd.read_csv(str(args.copyprobs_file_imputed), sep=" ", dtype=types_dict)
    anc_copyprobs.set_index("0", inplace=True)
    anc_copyprobs.columns = anc_copyprobs.columns.tolist()[::-1]  # Reverse column names because they are mis-labelled;
#  so when correct, column names are descending
else:  # For correctly labelled cols in copyprobs
    print("You have indicated cols in copyprobs files are correctly ordered")
    col_names = pd.read_csv(str(args.copyprobs_file_imputed), nrows=0).columns
    types_dict = {'ID': str}
    types_dict.update({col: 'float16' for col in col_names if col not in types_dict})
    anc_copyprobs = pd.read_csv(str(args.copyprobs_file_imputed), dtype=types_dict)
    anc_copyprobs.set_index("ID", inplace=True)
copyprobs_haps = list()
for h in range(int(anc_copyprobs.shape[0] / 2)):
    copyprobs_haps.extend([1, 2])
print("*** Successfully loaded copyprobs file for " + args.anc + " chr" + str(args.chr) + "***")
variants = pd.read_csv(args.variants,
                       usecols=['SNP', 'CHR', 'BP', 'other_allele', 'effect_allele', 'effect_allele_frequency', 'P', 'OR'],
                       dtype={'SNP': 'string', 'CHR': 'int8', 'BP': 'int32', 'other_allele': 'string',
                              'effect_allele': 'string', 'effect_allele_frequency': 'float', 'P': 'float', 'OR': 'float'},
                       sep=" ")
variants = variants.loc[variants['CHR'] == int(args.chr)]
results_list = list()
# for p in [0.0001, 0.001, 0.01]:
for p in [5e-8]:
    GWAS_variants_pval = variants.loc[variants['P'] < p]
    if GWAS_variants_pval.empty:
        print("No SNPs pass p-val threshold for " + str(args.chr))
        continue
    list_of_SNPs = GWAS_variants_pval['BP'].astype('str').drop_duplicates().tolist()
    anc_copyprobs_temp = anc_copyprobs[list_of_SNPs]
    if reverse_cols:
        anc_copyprobs_temp = anc_copyprobs_temp.reset_index().rename(columns={'0': 'ID'})
    else:
        anc_copyprobs_temp = anc_copyprobs_temp.reset_index()
    anc_copyprobs_temp['haps'] = copyprobs_haps
    if bootstrap:
        for bs in range(1000):
            anc_copyprobs_temp_bs = anc_copyprobs_temp.sample(n=anc_copyprobs_temp.shape[0], replace=True)
            analyse_anc(anc_copyprobs_temp_bs, str(args.anc), args.chr, p, bs, list_of_SNPs, GWAS_variants_pval)
    else:
        analyse_anc(anc_copyprobs_temp, str(args.anc), args.chr, p, bootstrap, list_of_SNPs, GWAS_variants_pval)

print("***Success! Now writing results to output file***")
if not os.path.exists("PRS_calculations_v4"):
    open("PRS_calculations_v4", 'a').close()
PRS_calculations = pd.DataFrame.from_records(results_list)
PRS_calculations.to_csv("PRS_calculations_v4", mode='a', header=False, index=False,
                        sep=" ")
