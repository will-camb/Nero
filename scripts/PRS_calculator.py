import pandas as pd
import argparse
import os
from collections import defaultdict


def analyse_anc(anc_copyprobs, anc, chrom):
    global PRS
    anc_copyprobs = anc_copyprobs[[str(a) for a in list_of_SNPs]]
    PRS = 0
    for i in list_of_SNPs:
        anc_copyprobs_i = pd.DataFrame(zip(anc_copyprobs.index, haps, anc_copyprobs[str(i)]))
        phase_i = pd.DataFrame(zip(phase.index, haps, phase[i]))
        df = pd.merge(anc_copyprobs_i, phase_i, left_on=[0, 1], right_on=[0, 1]).set_index(0)
        df.columns = ['haps', str(anc), 'phase']
        minor = df.loc[df['phase'] == '0'][anc].sum() #Sum painting prob for minor allele
        major = df.loc[df['phase'] == '1'][anc].sum() #Sum painting prob for major allele
        anc_maf = minor / (minor + major + 1)
        Beta = GWAS_n.loc[GWAS_n['pos'] == i].beta.item()
        vector = vectors[vectors['start.pos'].astype(int) == i]['vector'].item()
        PRS += anc_maf * Beta * vector
    anc_dict[anc].append([chrom, PRS])


parser = argparse.ArgumentParser()
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
args = parser.parse_args()
# Load Neale Lab effect size for this phenotype
GWAS = pd.read_csv(args.file_name, sep='\t')
GWAS[['chr', 'pos', '1', '2']] = GWAS['variant'].str.split(':', expand=True)
GWAS = GWAS[GWAS['pval'] <= 0.05]
if not os.path.exists("PRS_calculations"):
    open("PRS_calculations", 'a').close()
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
idfile = idfile[~idfile.index.isin(idfile.tail(318).index.tolist())]
anc_dict = defaultdict(list)
PRS = 0
print("Done loading of GWAS results for " + str(args.file_name))
for n in range(1, 23):
    print("Processing chromosome " + str(n))
    positions = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase",
                            skiprows=2,
                            nrows=1,
                            sep=" ",
                            header=None).T.drop(0)
    vectors = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase",
                          skiprows=2,
                          nrows=1,
                          sep=" ",
                          header=None).T.drop(0)
    vectors.columns = ['start.pos']
    vectors['previous'] = vectors['start.pos'].shift(1)
    vectors['next'] = vectors['start.pos'].shift(-1)
    vectors = vectors.fillna(0)
    vectors['vector'] = (vectors.next - vectors.previous) / 2
    vectors.loc[vectors.index[-1], 'vector'] = 0
    vectors.loc[vectors.index[-1], 'vector'] = vectors.vector.mean()
    vectors = vectors[['start.pos', 'vector']]
    vectors['start.pos'] = vectors['start.pos'].astype(str)
    GWAS_n = GWAS[GWAS['chr'] == str(n)]
    GWAS_n = GWAS_n.loc[GWAS_n['pos'].astype(str).isin(positions[0].astype(str).tolist())]
    GWAS_n.pos = pd.to_numeric(GWAS_n.pos)
    list_of_SNPs = GWAS_n['pos'].tolist()
    phase = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase", skiprows=3, header=None)
    phase = phase[~phase.index.isin(phase.tail(636).index.tolist())]
    phase = phase[0].dropna().apply(lambda x: pd.Series(list(x)))
    phase.columns = positions[0]
    phase = phase[list_of_SNPs]
    phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
    haps = list()
    for h in range(int(phase.shape[0] / 2)):
        haps.extend([1, 2])
    for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
        if anc == 'CHG':
            print("Calculating MAF/Beta/vectors for CHG, chr" + str(n))
            anc_copyprobs = pd.read_csv("temp.CHG." + str(n) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", index_col=0)
            analyse_anc(anc_copyprobs, anc, n)
        elif anc == 'EHG':
            print("Calculating MAF/Beta/vectors for EHG, chr" + str(n))
            anc_copyprobs = pd.read_csv("temp.EHG." + str(n) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", index_col=0)
            analyse_anc(anc_copyprobs, anc, n)
        elif anc == 'Farmer':
            print("Calculating MAF/Beta/vectors for Farmer, chr" + str(n))
            anc_copyprobs = pd.read_csv("temp.Farmer." + str(n) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", index_col=0)
            analyse_anc(anc_copyprobs, anc, n)
        elif anc == 'African':
            print("Calculating MAF/Beta/vectors for African, chr" + str(n))
            anc_copyprobs = pd.read_csv("temp.African." + str(n) + ".master_all_copyprobsperlocus.txt.gz", sep=" ",
                                        index_col=0)
            analyse_anc(anc_copyprobs, anc, n)
        elif anc == 'EastAsian':
            print("Calculating MAF/Beta/vectors for EastAsian, chr" + str(n))
            anc_copyprobs = pd.read_csv("temp.EastAsian." + str(n) + ".master_all_copyprobsperlocus.txt.gz", sep=" ",
                                        index_col=0)
            analyse_anc(anc_copyprobs, anc, n)
        elif anc == 'WHG':
            print("Calculating MAF/Beta/vectors for WHG, chr" + str(n))
            anc_copyprobs = pd.read_csv("temp.WHG." + str(n) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", index_col=0)
            analyse_anc(anc_copyprobs, anc, n)
        elif anc == 'Yamnaya':
            print("Calculating MAF/Beta/vectors for Yamnaya, chr" + str(n))
            anc_copyprobs = pd.read_csv("temp.Yamnaya." + str(n) + ".master_all_copyprobsperlocus.txt.gz", sep=" ", index_col=0)
            analyse_anc(anc_copyprobs, anc, n)
print("*** Success! Now calculating PRS per anc and writing to file ***")
# PRS_dict = dict()
# for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
#     PRS = 0
#     for j, k, l in anc_dict[anc]:
#         PRS += j * k * l
#     PRS_dict[anc] = PRS
PRS_dict = dict()
for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
    PRS_per_anc = 0
    for chrom, PRS in anc_dict[anc]:
        PRS_per_anc += PRS
    PRS_dict[anc] = PRS_per_anc

PRS_dict['Phenotype'] = args.file_name
PRS_calculations_temp = pd.DataFrame([PRS_dict])
PRS_calculations_temp = PRS_calculations_temp[["Phenotype", "CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]]
PRS_calculations_temp.to_csv('PRS_calculations', mode='a', index=False, header=False, sep=" ")
print("*** Success! See PRS_calculations for results ***")
