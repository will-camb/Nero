import pandas as pd
import wget
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_directory",
                    help="directory for all_copyprobsperlocus files; should be named in form n.all_copyprobsperlocus.txt.gz",
                    required=True)
parser.add_argument("-phasefile_directory",
                    help="directory for phasefiles; should be named in form n.merged.phase",
                    required=True)
parser.add_argument("-o",
                    help="location to save all_copyprobsperlocus.out to",
                    required=True)
args = parser.parse_args()

#Download Neale Lab effect size estimates from Dropbox
url = 'https://www.dropbox.com/s/x02h89mis8fe6md/age.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0'
output = os.getcwd() + 'temp.tsv.gz'
filename = wget.download(url, out=output)
#Load GWAS file, split to create new pos column
GWAS = pd.read_csv(filename, sep='\t')
GWAS[['chr','pos','1','2']] = GWAS['variant'].str.split(':',expand=True)

output_file_name = args.chr + '.all_copyprobsperlocus.txt'
if not os.path.isfile(os.path.join(args.o, file_name)):
    all_copyprobsperlocusDF = pd.DataFrame(columns=pd.read_csv(args.phasefile, nrows=2, header=None).iloc[1])
    all_copyprobsperlocusDF.to_csv(os.path.join(args.o, file_name), sep=' ', index=False)

for n in range(1,23):
    copyprobs = pd.read_csv(str(args.copyprobs_directory) + "/" + str(n) + ".all_copyprobsperlocus.txt.gz", header=None, skiprows=1, index_col=0)
    copyprobs = copyprobs[2].dropna().apply(lambda x: pd.Series(list(x)))
    positions = pd.read_csv(str(args.phasefile_directory) + "/" + str(n) + ".merged.phase", skiprows=2, nrows=1, sep=" ", header=None).T.drop(0)
    copyprobs.columns = positions[0]
    #Select every 10th row starting at row 9, i.e. WHG ancestry
    WHG = copyprobs.iloc[8::10]
    # CHG =
    # Yamnaya =
    # EHG =
    # Farmer =

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
    list_MAF_Beta = list()

    for i in list_of_SNPs:
        df = pd.DataFrame(list(zip(WHG[i].tolist(), phase[i].tolist())), columns=['WHG', 'phase'])
        df['WHG'] = pd.to_numeric(df['WHG'])
        df['phase'] = pd.to_numeric(df['phase'])
        df = df.loc[df['WHG'] >= 6]
        if df.shape[0] > 5:
            MAF = df.loc[df['phase'] == 0].shape[0] / df.shape[0]
            Beta = GWAS17.loc[GWAS17['pos'] == i].beta.item()
            list_MAF_Beta.append((MAF, Beta))

    PRS = 0
    for j, k in list_MAF_Beta:
        PRS += j * k



os.remove("temp.tsv.gz")