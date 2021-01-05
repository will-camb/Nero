import pandas as pd
import argparse
import os
from collections import defaultdict


def analyse_anc(anc_copyprobs, anc):
    if anc != 'Farmer':
        anc_copyprobs = anc_copyprobs[2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
    anc_copyprobs.columns = positions[0].tolist()
    for i in list_of_SNPs:
        anc_copyprobs_i = pd.DataFrame(zip(anc_copyprobs.index, haps, anc_copyprobs[i]))
        phase_i = pd.DataFrame(zip(phase.index, haps, phase[i]))
        df = pd.merge(anc_copyprobs_i, phase_i, left_on=[0,1], right_on=[0,1]).set_index(0)
        df.columns=['haps', str(anc), 'phase']
        df1 = df.loc[df[anc]>=6]
        if df1.shape[0] > 1:
            MAF = df1.loc[df1['phase']=='0'].shape[0]/df1.shape[0]
            Beta = GWAS_n.loc[GWAS_n['pos']==i].beta.item()
            anc_dict[anc].append([MAF,Beta])


parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_directory",
                    help="Directory of all_copyprobsperlocus files; should be named in form n.all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-phasefile_directory",
                    help="Directory of phasefiles; should be named in form chr#.merged.phase",
                    required=True)
parser.add_argument("-idfile",
                    help="Directory of idfile used in chromopainter ('individual pop 1')",
                    required=True)
parser.add_argument("-file_name",
                    help="Phenotype being looked at; this should be the file name in the Neale Lab google sheet; NB extension must be .gz and not.bgz!",
                    required=True)
args = parser.parse_args()
# Load Neale Lab effect size for this phenotype
GWAS = pd.read_csv(args.file_name, sep='\t')
GWAS[['chr', 'pos', '1', '2']] = GWAS['variant'].str.split(':', expand=True)
# If it doesn't already exist, make master PRS_calculations file
if not os.path.exists("PRS_calculations"):
    open("PRS_calculations", 'a').close()
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
idfile = idfile[~idfile.index.isin(idfile.tail(318).index.tolist())]
anc_dict = defaultdict(list)
for n in range(1, 23):
    print("Processing chromosome" + str(n))
    copyprobs = pd.read_csv(str(args.copyprobs_directory) + "/" + str(n) + ".all_copyprobsperlocus.txt.gz", header=None,
                            skiprows=1, index_col=0)
    positions = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase", skiprows=2, nrows=1,
                            sep=" ", header=None).T.drop(0)
    # Select chr, and SNPs that have been painted from GWAS file
    GWAS_n = GWAS[GWAS['chr'] == str(n)]
    GWAS_n = GWAS_n.loc[GWAS_n['pos'].astype(str).isin(positions[0].astype(str).tolist())]
    # Select SNPs with a p-value less than 0.05
    GWAS_n = GWAS_n[GWAS_n['pval'] <= 0.05]
    GWAS_n.pos = pd.to_numeric(GWAS_n.pos)
    list_of_SNPs = GWAS_n['pos'].tolist()
    # Load phase file to get genotype for each hap
    # NB Need to change this to ensure ordering of haps is the same in copyprobs and phasefile
    phase = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase", skiprows=3, header=None)
    phase = phase[0].dropna().apply(lambda x: pd.Series(list(x)))
    # Take only UKBB individuals in phase
    phase = phase[~phase.index.isin(phase.tail(636).index.tolist())]
    phase.columns = positions[0]
    phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
    haps = list()
    for h in range(int(phase.shape[0] / 2)):
        haps.extend([1, 2])
    for anc in ["CHG", "EHG", "Farmer", "WHG", "Yamnaya"]:
        if anc == 'CHG':
            analyse_anc(copyprobs.iloc[0::10], anc)
        elif anc == 'EHG':
            analyse_anc(copyprobs.iloc[1::10], anc)
        elif anc == 'Farmer':
            anc_copyprobs = copyprobs.iloc[5::10][2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                            copyprobs.iloc[4::10][2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                            copyprobs.iloc[2::10][2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric) + \
                            copyprobs.iloc[3::10][2].dropna().apply(lambda x: pd.Series(list(x))).apply(pd.to_numeric)
            analyse_anc(anc_copyprobs, anc)
        elif anc == 'WHG':
            analyse_anc(copyprobs.iloc[8::10], anc)
        elif anc == 'Yamnaya':
            analyse_anc(copyprobs.iloc[9::10], anc)
PRS_dict = dict()
for anc in ["CHG", "EHG", "Farmer", "WHG", "Yamnaya"]:
    PRS = 0
    for j, k in anc_dict[anc]:
        PRS += j * k
    PRS_dict[anc] = PRS
PRS_dict['Phenotype'] = args.file_name
PRS_calculations_temp = pd.DataFrame([PRS_dict])
PRS_calculations_temp = PRS_calculations_temp[["Phenotype", "CHG", "EHG", "Farmer", "WHG", "Yamnaya"]]
PRS_calculations_temp.to_csv('PRS_calculations', mode='a', index=False, header=False, sep=" ")
# NB calling all tmp files the same will prevent parallelisation
