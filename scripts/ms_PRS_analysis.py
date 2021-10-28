import pandas as pd
import argparse
import os
import numpy as np


def analyse_anc_phased(merged_phase_copyprobs_temp, anc, chrom, pval, iteration, phenotype_file):
    global list_of_SNPs
    global list_of_SNPs_phase
    global rsID_mapping
    global anc_copyprobs_temp
    skipped_snps = 0
    if len(list_of_SNPs_phase) != 0:
        phase_sum = merged_phase_copyprobs_temp.loc[:, (slice(None), 'phase')].sum().tolist()
        phase_sum_dict = {}
        for n, i in enumerate(list_of_SNPs_phase):
            phase_sum_dict[i] = phase_sum[n]
    alt_mapping = {}
    output = pd.DataFrame(index=list_of_SNPs, columns=['alt', 'ref'])
    for i in list_of_SNPs:
        if str(i) in list_of_SNPs_phase:  # for phased snps
            print(str(i) + " is phased and painted (in phasefile)")
            try:
                EUR_maf = GWAS_variants.loc[GWAS_variants['POS (hg19)'] == i].EUR.item()
            except ValueError:
                skipped_snps += 1
                list_of_SNPs.remove(i)
                print("Couldn't find EUR_maf for position:" + str(i) + ", so skipping")
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
            output.at[i, 'alt'] = alt_sum
            output.at[i, 'ref'] = ref_sum
        else:  # for non-phased snps
            print(str(i) + " is not phased, so looking at homozygous individuals")
            try:
                EUR_maf = GWAS_variants.loc[GWAS_variants['POS (hg19)'] == i].EUR.item()
                rsid = rsID_mapping.loc[rsID_mapping['Position'] == int(i)]["rsID"].item()
                samples_0 = pd.read_csv(args.output_files + str(rsid) + ".hom.0.samples", header=None)
                samples_1 = pd.read_csv(args.output_files + str(rsid) + ".hom.1.samples", header=None)
            except (ValueError, IOError):
                skipped_snps += 1
                list_of_SNPs.remove(i)
                print("Couldn't find rsID/EUR_maf/output file for position:" + str(i) + ", so skipping")
                continue
            anc_copyprobs_temp_i = anc_copyprobs_temp[['ID', str(i)]]
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
            elif alt_mapping[i] == 1:
                alt_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_1[0].tolist())][i].sum()
                ref_sum = anc_copyprobs_temp_i.loc[anc_copyprobs_temp_i['ID'].isin(samples_0[0].tolist())][i].sum()
            output.at[i, 'alt'] = alt_sum
            output.at[i, 'ref'] = ref_sum
    output['maf'] = output['alt'] / (output['alt'] + output['ref'])
    output['maf'] = output['maf'].fillna(0)
    output = output.reset_index().rename(columns={'index': 'pos'})
    maf_GWAS = pd.merge(output, GWAS_variants, left_on='pos', right_on='POS (hg19)')
    maf_GWAS['beta'] = np.log(maf_GWAS['OR'])  # Need to transform OR to beta
    maf_GWAS['maf_x_beta'] = maf_GWAS['maf'] * maf_GWAS['beta']
    PRS = maf_GWAS['maf_x_beta'].sum()
    number_of_SNPs = len(list_of_SNPs)
    results_list.append(
        [args.copyprobs_file_imputed, phenotype_file, anc, chrom, PRS, number_of_SNPs, skipped_snps, pval, iteration,
         maf_GWAS['maf'].tolist(), maf_GWAS['beta'].tolist(), list_of_SNPs])


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
                    help="Directory of idfile used in chromopainter ('individual pop 1')",
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
args = parser.parse_args()
if args.bootstrap == 'True':
    bootstrap = True
else:
    bootstrap = False

# Read in copyprobs
col_names = pd.read_csv(str(args.copyprobs_file_imputed), sep=" ", nrows=0).columns
types_dict = {'0': str}
types_dict.update({col: 'float' for col in col_names if col not in types_dict})
anc_copyprobs = pd.read_csv(str(args.copyprobs_file_imputed), sep=" ", dtype=types_dict)
anc_copyprobs.set_index("0", inplace=True)
anc_copyprobs.columns = anc_copyprobs.columns.tolist()[::-1]  # Reverse column names because they are mis-labelled;
#  so when correct, column names are descending
copyprobs_haps = list()
for h in range(int(anc_copyprobs.shape[0] / 2)):
    copyprobs_haps.extend([1, 2])
print("*** Successfully loaded copyprobs file for " + args.anc + " chr" + str(args.chr) + "***")
# Read in idfile
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
# Read in phasefile
phase = pd.read_csv(str(args.phasefile), header=None, sep=" ", dtype='int8')
phase.columns = pd.read_csv(str(args.copyprobs_file_original), sep=" ", nrows=0).columns.tolist()[1:]
# Phasefile cols should be ascending
phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
phase_haps = list()
for h in range(int(phase.shape[0] / 2)):
    phase_haps.extend([1, 2])
# Load variants file to find which is alt allele (not always minor!)
variants = pd.read_csv("ms_all_snps_annotated.csv",
                       dtype={'chr': 'string', 'POS (hg19)': 'string', 'A1': 'string', 'A2': 'string', 'OR': 'float',
                              'EUR': 'float', 'P (joined)': 'float'})
variants = variants.loc[variants['chr'] == str(args.chr)]
ldetect = pd.read_csv("/willerslev/ukbiobank/painting_results_aggregate/PRS_calculation/ldetect-data/EUR/fourier_ls-chr"
                      + str(args.chr) + ".bed", sep='\t')
ldetect.rename(columns=lambda x: x.strip(), inplace=True)
rsID_mapping = pd.read_csv("ms_all_snps_rsids_annotated.csv")
rsID_mapping = rsID_mapping.loc[rsID_mapping['chr'].astype(int) == int(args.chr)]
print(rsID_mapping.shape)
results_list = list()
# LOOP through GWAS files
phenotypes = pd.read_csv(args.phenotypes, header=None, dtype=str)[0].tolist()
phenotypes = [s.strip() for s in phenotypes]
for file in phenotypes:
    print("***Looking at " + str(file) + "***")
    GWAS_variants = variants
    print(GWAS_variants.shape)
    # GWAS_variants = variants.loc[variants['POS (hg19)'].isin(phase.columns.tolist())]
    if GWAS_variants.empty:
        print("GWAS_variants is empty!")
        continue
    for p in [5e-8]:
        GWAS_variants_pval = GWAS_variants.loc[GWAS_variants['P (joined)'] < p]
        if GWAS_variants_pval.empty:
            print("No SNPs pass p-val threshold for " + str(args.chr))
            continue
        #  Drop SNPs from GWAS_variants that are neither painted or imputed
        for index, row in GWAS_variants.iterrows():
            try:
                rsID = rsID_mapping.loc[rsID_mapping['Position']==row['POS (hg19)']]['rsID'].item()
            except ValueError:
                GWAS_variants.drop(index, inplace=True)  # For SNPs with no rsID found
                continue
            if not row['POS (hg19)'].isin(anc_copyprobs.columns.tolist()):  # Check dtypes!
                if not os.path.exists("/willerslev/ukbiobank/painting_results_aggregate/PRS_calculation/"
                                      "multiple_sclerosis/non_phased_snps/output_files/" + rsID + ".hom.1.samples"):
                    GWAS_variants.drop(index, inplace=True)
        GWAS_variants.reset_index(inplace=True)
        best_per_block = pd.DataFrame(columns=GWAS_variants.columns)
        for index, row in ldetect.iterrows():
            block = GWAS_variants_pval.loc[(GWAS_variants_pval['POS (hg19)'].astype(int) >= row['start']) & (
                    GWAS_variants_pval['POS (hg19)'].astype(int) < row['stop'])].reset_index()
            try:
                block = block.reset_index()
                # best_per_block.loc[index] = block.iloc[block['beta'].abs().idxmax()]  # Take highest beta
                best_per_block.loc[index] = block.iloc[block['P (joined)'].idxmin()]  # Take lowest p-val
            except ValueError:  # For when there are no SNPs that pass p-val threshold in the block
                continue
        list_of_SNPs = best_per_block['POS (hg19)'].drop_duplicates().tolist()
        print(len(list_of_SNPs))
        phase_temp = phase.loc[:, phase.columns.intersection([str(x) for x in list_of_SNPs])]
        list_of_SNPs_phase = phase_temp.columns.tolist()
        print(len(list_of_SNPs_phase))
        phase_temp = phase_temp.reset_index().rename(columns={'index': 'ID'})
        phase_temp['haps'] = phase_haps
        anc_copyprobs_temp = anc_copyprobs[list_of_SNPs]
        anc_copyprobs_temp = anc_copyprobs_temp.reset_index().rename(columns={'0': 'ID'})
        print(anc_copyprobs_temp.shape)
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
                analyse_anc_phased(temp, str(args.anc), args.chr, p, bs, file)
        else:
            analyse_anc_phased(merged_phase_copyprobs, str(args.anc), args.chr, p, bootstrap, file)

print("***Success! Now writing results to output file***")
if not os.path.exists("PRS_calculations_bootstrapped_ms_" + str(args.phenotypes)):
    open("PRS_calculations_bootstrapped_ms_" + str(args.phenotypes), 'a').close()
PRS_calculations = pd.DataFrame.from_records(results_list)
PRS_calculations.to_csv("PRS_calculations_bootstrapped_ms_" + str(args.phenotypes), mode='a', header=False, index=False,
                        sep=" ")
