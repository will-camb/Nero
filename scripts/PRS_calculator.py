import pandas as pd
import argparse
import os
from collections import defaultdict
import random


def analyse_anc(anc_copyprobs, anc, chrom):
    global PRS
    anc_copyprobs = anc_copyprobs[[str(a) for a in list_of_SNPs]]
    if anc_copyprobs.empty:
        print("Error 1!")
    PRS = 0
    for i in list_of_SNPs:
        anc_copyprobs_i = pd.DataFrame(zip(anc_copyprobs.index, haps, anc_copyprobs[str(i)]))
        phase_i = pd.DataFrame(zip(phase.index, haps, phase[i]))
        df = pd.merge(anc_copyprobs_i, phase_i, left_on=[0, 1], right_on=[0, 1]).set_index(0)
        df.columns = ['haps', str(anc), 'phase']
        if df.empty:
            print("Error 2!")
        minor = df.loc[df['phase'] == '0'][anc].sum() #Sum painting prob for minor allele
        major = df.loc[df['phase'] == '1'][anc].sum() #Sum painting prob for major allele
        anc_maf = minor / (minor + major)
        Beta = GWAS_n.loc[GWAS_n['pos'] == i].beta.item()
        PRS += anc_maf * Beta
    number_of_SNPs = len(list_of_SNPs)
    anc_dict[anc].append([chrom, PRS, number_of_SNPs])


parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_directory",
                    help="Directory of copyprobs files; should be named in form anc.chr.master_all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-phasefile_directory",
                    help="Directory of phasefiles; should be named in form chr#.merged.phase",
                    required=True)
parser.add_argument("-idfile",
                    help="Directory of idfile used in chromopainter ('individual pop 1')",
                    required=True)
parser.add_argument("-file_name",
                    help="Phenotype being looked at; this should be the file name in the Neale Lab google sheet; NB "
                         "extension must be .gz and not.bgz!",
                    required=True)
parser.add_argument("-pruned_dir",
                    help="Directory containing pruned.in files",
                    required=True)
parser.add_argument("-pruned_in",
                    help="A file containing to SNPs to keep after LD pruning with a specific threshold eg 05.pruned.in",
                    required=True)
parser.add_argument("-p_value",
                    help="The p-value cutoff for including effect sizes from Neale Lab GWAS",
                    required=True)
args = parser.parse_args()
GWAS = pd.read_csv(args.file_name, sep='\t', usecols=['beta', 'pval', 'variant'])
GWAS[['chr', 'pos', '1', '2']] = GWAS['variant'].str.split(':', expand=True)
GWAS = GWAS[GWAS['pval'] <= float(args.p_value)]
if not os.path.exists("PRS_calculations"):
    open("PRS_calculations", 'a').close()
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
idfile = idfile[~idfile.index.isin(idfile.tail(318).index.tolist())]
anc_dict = defaultdict(list)
PRS = 0
print("Done loading of GWAS results for " + str(args.file_name))
for n in random.sample(range(1,23), 22):
# for n in range(21, 23):
    print("Processing chromosome " + str(n))
    positions = pd.read_csv("/willerslev/ukbiobank/phasefiles/" + str(n) + ".merged.phase",
                            skiprows=2,
                            nrows=1,
                            sep=" ",
                            header=None).T.drop(0)
    GWAS_n = GWAS[GWAS['chr'] == str(n)]
    GWAS_n = GWAS_n.loc[GWAS_n['pos'].astype(str).isin(positions[0].astype(str).tolist())]
    GWAS_n.pos = pd.to_numeric(GWAS_n.pos)
    pruned_in = pd.read_csv(str(args.pruned_dir) + "/" + str(n) + "." + str(args.pruned_in), header=None)
    GWAS_n = GWAS_n.loc[GWAS_n['pos'].isin(pruned_in[0].astype(str).tolist())]
    if GWAS_n.empty:
        print("No SNPs for chr" + str(n) + " that pass LD and p-value thresholds, \
        moving onto next chromosome")
        continue
    list_of_SNPs = GWAS_n['pos'].tolist()
    # phase = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase", skiprows=3, header=None)
    phase = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase.gz",
                        header=None,
                        sep=" ",
                        dtype='int8')
    # phase = phase[~phase.index.isin(phase.tail(636).index.tolist())]
    # phase = phase[0].dropna().apply(lambda x: pd.Series(list(x)))# Memory error here
    phase.columns = positions[0]
    phase = phase[list_of_SNPs]
    phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
    haps = list()
    for h in range(int(phase.shape[0] / 2)):
        haps.extend([1, 2])
    for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
        print("Calculating MAF/Beta for " + str(anc) + ", chr" + str(n))
        col_names = pd.read_csv(str(args.copyprobs_directory) + "/" + str(anc) + "." + str(n) + ".master_all_copyprobsperlocus.txt.gz",
                                sep=" ",
                                nrows=0).columns
        types_dict = {'0': str}
        types_dict.update({col: 'int8' for col in col_names if col not in types_dict})
        anc_copyprobs = pd.read_csv(str(args.copyprobs_directory) + "/" + str(anc) + "." + str(n) + ".master_all_copyprobsperlocus.txt.gz",
                                    sep=" ",
                                    dtype=types_dict)
        anc_copyprobs.set_index("0", inplace=True)
        analyse_anc(anc_copyprobs, anc, n)

print("*** Success! Now calculating PRS per anc and writing to file ***")
PRS_dict = dict()
for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
    PRS_per_anc = 0
    total_number_of_SNPs = 0
    for chrom, PRS, number_of_SNPs in anc_dict[anc]:
        PRS_per_anc += PRS
        total_number_of_SNPs += number_of_SNPs
    PRS_dict[anc] = PRS_per_anc
    PRS_dict[anc + "_perSNP"] = PRS_per_anc/total_number_of_SNPs
PRS_dict['Phenotype'] = args.file_name
PRS_dict['LD_pruned'] = args.pruned_in
PRS_dict['p_value'] = args.p_value
PRS_calculations_temp = pd.DataFrame([PRS_dict])
PRS_calculations_temp = PRS_calculations_temp[["Phenotype", "LD_pruned", "p_value", "CHG", "EHG", "Farmer",
                                               "African", "EastAsian", "WHG", "Yamnaya", "CHG_perSNP", "EHG_perSNP",
                                               "Farmer_perSNP", "African_perSNP", "EastAsian_perSNP", "WHG_perSNP",
                                               "Yamnaya_perSNP"]]
PRS_calculations_temp.to_csv('PRS_calculations', mode='a', index=False, header=True, sep=" ")
print("*** Success! See PRS_calculations for results ***")
