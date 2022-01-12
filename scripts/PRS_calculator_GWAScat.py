import pandas as pd
import argparse
import os
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'


def analyse_anc(merged_phase_copyprobs_temp, anc, chrom, pval, iteration, phenotype_file, list_of_SNPs_):
    global list_of_SNPs_phase
    global GWAS_variants_pval
    global anc_copyprobs_temp
    skipped_snps = 0
    if len(list_of_SNPs_phase) != 0:
        phase_sum = merged_phase_copyprobs_temp.loc[:, (slice(None), 'phase')].sum().tolist()
        phase_sum_dict = {}
        for n, i in enumerate(list_of_SNPs_phase):
            phase_sum_dict[i] = phase_sum[n]
    alt_mapping = {}
    output = pd.DataFrame(columns=['pos', 'alt', 'ref', 'alt_average', 'ref_average', 'OR', 'rsID'])
    for n, i in enumerate(list_of_SNPs_):  # i is string!
        if str(i) in list_of_SNPs_phase:  # for phased snps
            print(i + " is phased and painted (in phasefile)")
            try:
                GWAS_variants_pval.reset_index(drop=True, inplace=True)
                EUR_maf = GWAS_variants_pval.loc[(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()]['RISK ALLELE FREQUENCY'].item()
                OR = GWAS_variants_pval.loc[(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()]['OR or BETA'].item()
                rsid = GWAS_variants_pval.loc[(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()]['SNPS']
                GWAS_variants_pval.drop([(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()], inplace=True)
            except Exception as e:
                skipped_snps += 1
                print("Couldn't find EUR_maf for position:" + str(i) + ", so skipping: " + str(e))
                continue
            #  Check if alt/effect allele is major or minor
            if EUR_maf == -1:  # For all of these, alt/effect=minor
                if phase_sum_dict[i] > (len(copyprobs_haps) / 2):
                    alt_mapping[i] = 0
                else:
                    alt_mapping[i] = 1
            elif EUR_maf <= 0.5:
                #  minor is alt/effect
                #  Find if minor is 1 or 0 in phase
                if phase_sum_dict[i] > (len(copyprobs_haps) / 2):
                    alt_mapping[i] = 0
                else:
                    alt_mapping[i] = 1
            elif EUR_maf > 0.5:
                if phase_sum_dict[i] > (len(copyprobs_haps) / 2):
                    alt_mapping[i] = 1
                else:
                    alt_mapping[i] = 0
            alt_sum = merged_phase_copyprobs_temp.loc[merged_phase_copyprobs_temp[i]['phase'] == alt_mapping[i]].loc[:,
                      (i, 'copyprobs')].sum()
            ref_sum = merged_phase_copyprobs_temp.loc[merged_phase_copyprobs_temp[i]['phase'] != alt_mapping[i]].loc[:,
                      (i, 'copyprobs')].sum()
            alt_average = alt_sum / merged_phase_copyprobs_temp.loc[merged_phase_copyprobs_temp[i]['phase'] == alt_mapping[i]].loc[:,
                      (i, 'copyprobs')].shape[0]
            ref_average = ref_sum / merged_phase_copyprobs_temp.loc[merged_phase_copyprobs_temp[i]['phase'] != alt_mapping[i]].loc[:,
                      (i, 'copyprobs')].shape[0]
            output.loc[n] = [i, alt_sum, ref_sum, alt_average, ref_average, OR, rsid]
        else:  # for non-phased snps
            print(str(i) + " is not phased, so looking at homozygous individuals")
            try:
                GWAS_variants_pval.reset_index(drop=True, inplace=True)
                EUR_maf = GWAS_variants_pval.loc[(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()]['RISK ALLELE FREQUENCY'].item()
                OR = GWAS_variants_pval.loc[(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()]['OR or BETA'].item()
                rsid = GWAS_variants_pval.loc[(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()]['SNPS']
                GWAS_variants_pval.drop([(GWAS_variants_pval['CHR_POS'].astype(int) == int(i)).idxmax()], inplace=True)
                # Get first instance of rsID; when there is more than 1 variant at the same position,
                # take first then delete row so we take second the second time
                samples_0 = pd.read_csv(args.output_files + str(rsid) + ".hom.0.samples", header=None)
                samples_1 = pd.read_csv(args.output_files + str(rsid) + ".hom.1.samples", header=None)
            except Exception as e:
                skipped_snps += 1
                print("Couldn't find rsID/EUR_maf/output file for position:" + str(i) + ", so skipping: " + str(e))
                continue
            try:
                anc_copyprobs_temp_i = anc_copyprobs_temp[['ID', str(i)]]
            except Exception as e:
                print("Couldn't find SNP with pos " + str(i) + " in anc_copyprobs: " + str(e))
            if EUR_maf == -1:  # For all of these, alt/effect=minor
                if samples_0.shape[0] < (len(copyprobs_haps) / 8):
                    alt_mapping[i] = 0
                else:
                    alt_mapping[i] = 1
            elif EUR_maf <= 0.5:
                #  minor is alt/effect
                if samples_0.shape[0] < (len(copyprobs_haps) / 8):
                    alt_mapping[i] = 0
                else:
                    alt_mapping[i] = 1
            elif EUR_maf > 0.5:
                if samples_0.shape[0] < (len(copyprobs_haps) / 8):
                    alt_mapping[i] = 1
                else:
                    alt_mapping[i] = 0
            if alt_mapping[i] == 0:
                alt_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_0[0].tolist())][i].sum()
                ref_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_1[0].tolist())][i].sum()
                alt_average = alt_sum / len(samples_0[0].tolist())
                ref_average = ref_sum / len(samples_1[0].tolist())
            elif alt_mapping[i] == 1:
                alt_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_1[0].tolist())][i].sum()
                ref_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_0[0].tolist())][i].sum()
                alt_average = alt_sum / len(samples_1[0].tolist())
                ref_average = ref_sum / len(samples_0[0].tolist())
            output.loc[n] = [i, alt_sum, ref_sum, alt_average, ref_average, OR, rsid]
    output['maf'] = output['alt'] / (output['alt'] + output['ref'])
    output['maf'] = output['maf'].fillna(0)
    output['maf_average'] = output['alt_average'] / (output['alt_average'] + output['ref_average'])
    output['maf_average'] = output['maf_average'].fillna(0)
    output['beta'] = np.log(output['OR'].astype(float))  # Need to transform OR to beta
    output['maf_x_beta'] = output['maf'] * output['beta']
    PRS = output['maf_x_beta'].sum()
    maf_GWAS = output
    number_of_SNPs = len(list_of_SNPs_)
    results_list.append(
        [args.copyprobs_file_imputed, phenotype_file, anc, chrom, PRS, number_of_SNPs, skipped_snps, pval, iteration,
         maf_GWAS['maf'].tolist(), maf_GWAS['maf_average'].tolist(), maf_GWAS['beta'].tolist(), maf_GWAS['pos'].tolist(), maf_GWAS['rsID'].tolist()])


parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file_original",
                    help="Should be named in form anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-copyprobs_file_imputed",
                    help="Should be named in form anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-phasefile",
                    help="Should be named in form chr#.merged.phase.gz",
                    required=True)
parser.add_argument("-idfile",
                    help="Directory of idfile used in chromopainter ('individual pop 1') for phasefile",
                    required=True)
parser.add_argument("-variants",
                    help="GWAS variants file containing SNP details - chr, pos, rsID, OR/BETA",
                    required=True)
parser.add_argument("-phenotypes",
                    help="File with list of phenotypes being looked at; NB extension must be .gz and not.bgz!",
                    required=True)
parser.add_argument("-output_files",
                    help="Location of output files from get_nonphased_samples.sh",
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
parser.add_argument("-LD_prune",
                    help="Whether to prune SNPs for LD - lowest p-value per Pickrell LD block; "
                         "can take two values: True or False",
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
if args.LD_prune == 'True':
    LD_prune = True
else:
    LD_prune = False

# Read in copyprobs
if reverse_cols:  # For copyprobs files with incorrectly ordered columns
    print("You have indicated cols in copyprobs files are incorrectly ordered")
    col_names = pd.read_csv(str(args.copyprobs_file_imputed), sep=" ", nrows=0).columns
    types_dict = {'0': str}
    types_dict.update({col: 'float' for col in col_names if col not in types_dict})
    anc_copyprobs = pd.read_csv(str(args.copyprobs_file_imputed), sep=" ", dtype=types_dict)
    anc_copyprobs.set_index("0", inplace=True)
    anc_copyprobs.columns = anc_copyprobs.columns.tolist()[::-1]  # Reverse column names because they are mis-labelled;
#  so when correct, column names are descending
else:  # For correctly labelled cols in copyprobs
    print("You have indicated cols in copyprobs files are correctly ordered")
    col_names = pd.read_csv(str(args.copyprobs_file_imputed), nrows=0).columns
    types_dict = {'ID': str}
    types_dict.update({col: 'float' for col in col_names if col not in types_dict})
    anc_copyprobs = pd.read_csv(str(args.copyprobs_file_imputed), dtype=types_dict)
    anc_copyprobs.set_index("ID", inplace=True)
copyprobs_haps = list()
for h in range(int(anc_copyprobs.shape[0] / 2)):
    copyprobs_haps.extend([1, 2])
print("*** Successfully loaded copyprobs file for " + args.anc + " chr" + str(args.chr) + "***")
# Read in idfile
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
# Read in phasefile
phase = pd.read_csv(str(args.phasefile), header=None, sep=" ", dtype='int8')
print("*** Successfully loaded phasefile file for " + args.anc + " chr" + str(args.chr) + "***")
if reverse_cols:
    phase.columns = pd.read_csv(str(args.copyprobs_file_original), sep=" ", nrows=0).columns.tolist()[1:]
else:
    phase.columns = pd.read_csv(str(args.copyprobs_file_original), sep=" ", nrows=0).columns.tolist()[1:][::-1]
# Phasefile cols should be ascending
phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
phase_haps = list()
for h in range(int(phase.shape[0] / 2)):
    phase_haps.extend([1, 2])
# Load variants file to find which is alt allele (not always minor!)
variants = pd.read_csv(args.variants,
                       dtype={'CHR_ID': 'string', 'CHR_POS': 'string', 'OR or BETA': 'float', 'P-VALUE': 'float',
                              'RISK ALLELE FREQUENCY': 'float'})
variants = variants.loc[variants['CHR_ID'] == str(args.chr)]
ldetect = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/PRS_calculation/ldetect-data/EUR/fourier_ls-chr"
                      + str(args.chr) + ".bed", sep='\t')
ldetect.rename(columns=lambda x: x.strip(), inplace=True)
results_list = list()
# LOOP through GWAS files
phenotypes = pd.read_csv(args.phenotypes, header=None, dtype=str)[0].tolist()
phenotypes = [s.strip() for s in phenotypes]
for file in phenotypes:
    print("***Looking at " + str(file) + "***")
    GWAS_variants = variants
    if GWAS_variants.empty:
        print("GWAS_variants is empty!")
        continue
    for p in [5e-8]:  # NB Should this be changed??
        GWAS_variants_pval = GWAS_variants.loc[GWAS_variants['P-VALUE'] < p]
        if GWAS_variants_pval.empty:
            print("No SNPs pass p-val threshold for " + str(args.chr))
            continue
        if LD_prune:  # Choose lowest p-value per LD block
            print("Finding lowest p-value SNP per independent LD block")
            best_per_block = pd.DataFrame(columns=GWAS_variants_pval.columns)
            for index, row in ldetect.iterrows():
                block = GWAS_variants_pval.loc[(GWAS_variants_pval['CHR_POS'].astype(int) >= row['start']) & (
                        GWAS_variants_pval['CHR_POS'].astype(int) < row['stop'])].reset_index()
                try:
                    block = block.reset_index()
                    # best_per_block.loc[index] = block.iloc[block['beta'].abs().idxmax()]  # Take highest beta
                    best_per_block.loc[index] = block.iloc[block['P-VALUE'].idxmin()]  # Take lowest p-val
                except ValueError:  # For when there are no SNPs that pass p-val threshold in the block
                    continue
            list_of_SNPs = best_per_block['CHR_POS'].drop_duplicates().tolist()
            list_of_SNPs_unique = list_of_SNPs
        else:
            print("Including all SNPs - not running LD pruning")
            list_of_SNPs_unique = GWAS_variants_pval['CHR_POS'].drop_duplicates().tolist()
            list_of_SNPs = GWAS_variants_pval['CHR_POS'].tolist()  # There are some repeats of positions
        phase_temp = phase.loc[:, phase.columns.intersection([str(x) for x in list_of_SNPs_unique])]
        list_of_SNPs_phase = phase_temp.columns.tolist()
        phase_temp = phase_temp.reset_index().rename(columns={'index': 'ID'})
        phase_temp['haps'] = phase_haps
        anc_copyprobs_temp = anc_copyprobs[list_of_SNPs_unique]
        if reverse_cols:
            anc_copyprobs_temp = anc_copyprobs_temp.reset_index().rename(columns={'0': 'ID'})
        else:
            anc_copyprobs_temp = anc_copyprobs_temp.reset_index()
        anc_copyprobs_temp['haps'] = copyprobs_haps
        merged_phase_copyprobs = pd.merge(phase_temp, anc_copyprobs_temp, on=['ID', 'haps'])
        reorder = [[str(x) + "_x", str(x) + "_y"] for x in list_of_SNPs_phase]
        reorder = [item for sublist in reorder for item in sublist]
        merged_phase_copyprobs = merged_phase_copyprobs[['ID', 'haps'] + reorder]
        iterables = [["ID"] + list_of_SNPs_phase, ["phase", "copyprobs"]]
        merged_phase_copyprobs.columns = pd.MultiIndex.from_product(iterables, names=['first', 'second'])
        merged_phase_copyprobs.set_index('ID', inplace=True)
        if bootstrap:
            for bs in range(1000):
                temp = merged_phase_copyprobs.sample(n=merged_phase_copyprobs.shape[0], replace=True)
                analyse_anc(temp, str(args.anc), args.chr, p, bs, file, list_of_SNPs)
        else:
            analyse_anc(merged_phase_copyprobs, str(args.anc), args.chr, p, bootstrap, file, list_of_SNPs)

print("***Success! Now writing results to output file***")
if not os.path.exists("PRS_calculations_" + str(args.phenotypes)):
    open("PRS_calculations_" + str(args.phenotypes), 'a').close()
PRS_calculations = pd.DataFrame.from_records(results_list)
PRS_calculations.to_csv("PRS_calculations_" + str(args.phenotypes), mode='a', header=False, index=False,
                        sep=" ")
