import pandas as pd
import argparse
import os
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_directory",
                    help="Directory of all_copyprobsperlocus files; should be named in form n.all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-phasefile_directory",
                    help="Directory of phasefiles; should be named in form chr#.merged.phase",
                    required=True)
parser.add_argument("-file_name",
                    help="Phenotype being looked at; this should be the file name in the Neale Lab google sheet; NB extension must be .gz and not.bgz!",
                    required=True)
args = parser.parse_args()
#Load Neale Lab effect size for this phenotype
GWAS = pd.read_csv(args.file_name, sep='\t')
GWAS[['chr','pos','1','2']] = GWAS['variant'].str.split(':',expand=True)
#If it doesn't already exist, make master PRS_calculations file
if not os.path.exists("PRS_calculations"):
    open("PRS_calculations", 'a').close()
anc_dict = defaultdict(list)
for n in range(1,23):
    copyprobs = pd.read_csv(str(args.copyprobs_directory) + "/" + str(n) + ".all_copyprobsperlocus.txt.gz", header=None, skiprows=1, index_col=0)
    # for final copyprobs file that has been cleaned:
    # copyprobs = pd.read_csv(str(args.copyprobs_directory) + "/" + str(n) + ".all_copyprobsperlocus.txt.gz", header=None, index_col=0)
    copyprobs = copyprobs[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
    positions = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase", skiprows=2, nrows=1, sep=" ", header=None).T.drop(0)
    copyprobs.columns = positions[0]
    #Select every 10th row starting at row 9, i.e. WHG ancestry
    CHG = copyprobs.iloc[0::10]
    EHG = copyprobs.iloc[1::10]
    Farmer = copyprobs.iloc[5::10] + copyprobs.iloc[4::10] + copyprobs.iloc[2::10] + copyprobs.iloc[3::10]
    WHG = copyprobs.iloc[8::10]
    Yamnaya = copyprobs.iloc[9::10]
    # Load phase file to get genotype for each hap
    phase = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase", skiprows=3, header=None)
    phase = phase[0].dropna().apply(lambda x: pd.Series(list(x)))
    #Take only UKBB individuals in phase
    phase = phase[~phase.index.isin(phase.tail(636).index.tolist())]
    phase.columns = positions[0]
    # Select chr, and SNPs that have been painted
    GWAS_n = GWAS[GWAS['chr'] == str(n)]
    GWAS_n = GWAS_n.loc[GWAS_n['pos'].astype(str).isin(positions[0].astype(str).tolist())]
    # Select SNPs with a p-value less than 0.05
    GWAS_n = GWAS_n[GWAS_n['pval'] <= 0.05]
    GWAS_n.pos = pd.to_numeric(GWAS_n.pos)
    list_of_SNPs = GWAS_n['pos'].tolist()
    for i in list_of_SNPs:
        df = pd.DataFrame(list(zip(CHG[i].tolist(), EHG[i].tolist(), Farmer[i].tolist(), WHG[i].tolist(), Yamnaya[i].tolist(), phase[i].tolist())), columns = ['CHG', 'EHG', 'Farmer', 'WHG', 'Yamnaya', 'phase'])
        for anc in ["CHG", "EHG", "Farmer", "WHG", "Yamnaya"]:
            df1 = df.loc[df[anc] >= 6]
            if df1.shape[0] > 1:
                MAF = df1.loc[df1['phase'] == '0'].shape[0] / df1.shape[0]
                Beta = GWAS_n.loc[GWAS_n['pos'] == i].beta.item()
                anc_dict[anc].append([MAF, Beta])
PRS_dict = dict()
for anc in ["CHG", "EHG", "Farmer", "WHG", "Yamnaya"]:
    PRS = 0
    for j,k in anc_dict[anc]:
        PRS += j*k
    PRS_dict[anc] = PRS
PRS_dict['Phenotype'] = args.file_name
PRS_calculations_temp = pd.DataFrame([PRS_dict])
PRS_calculations_temp=PRS_calculations_temp[["Phenotype", "CHG", "EHG", "Farmer", "WHG", "Yamnaya"]]
PRS_calculations_temp.to_csv('PRS_calculations', mode='a', index=False, header=False, sep=" ")
#NB calling all tmp files the same will prevent parallelisation

