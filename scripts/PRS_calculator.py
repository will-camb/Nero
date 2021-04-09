import pandas as pd
import argparse
import os

# Could print warning when maf is > 0.5?

def analyse_anc(anc_copyprobs_temp, anc, chrom, pval, LD_value, phenotype_file):
    global list_of_SNPs
    PRS = 0
    skipped_snps = 0
    for i in list_of_SNPs:
        anc_copyprobs_temp_i = pd.DataFrame(zip(anc_copyprobs_temp.index, haps, anc_copyprobs_temp[str(i)]))
        phase_i = pd.DataFrame(zip(phase_temp.index, haps, phase_temp[str(i)]))
        df = pd.merge(anc_copyprobs_temp_i, phase_i, left_on=[0, 1], right_on=[0, 1]).set_index(0)
        df.columns = ['haps', anc, 'phase']
        minor_df = df.loc[df['phase'] == 0]
        if minor_df.empty:
            minor = 0
        else:
            minor = minor_df[anc].sum()
        if isinstance(minor, str):
            print("Minor.dtype is string for " + str(i) + str(anc) + str(chrom))
            skipped_snps += 1
            continue
        major_df = df.loc[df['phase'] == 1]
        if major_df.empty:
            major = 0
        else:
            major = major_df[anc].sum()  # Sum painting prob for major allele
        if isinstance(major, str):
            print("Major.dtype is string for " + str(i) + str(anc) + str(chrom))
            skipped_snps += 1
            continue
        if minor + major != 0:
            anc_maf = minor/(minor+major)
        else:
            anc_maf = 0
        if GWAS_pruned.loc[GWAS_pruned['pos'] == i].empty:
            print("Can't find SNP " + str(i) + " for " + str(anc) + str(chrom))
            skipped_snps += 1
            continue
        try:
            beta = GWAS_pruned.loc[GWAS_pruned['pos'] == i].beta.item()
        except ValueError:  # This skips SNPs with more then one entry in the GWAS file
            print("ValueError for beta for SNP at position " + str(i) + " for " + str(anc) + " chr" + str(chrom) + "!")
            skipped_snps += 1
            continue
        PRS += anc_maf * beta
    number_of_SNPs = len(list_of_SNPs)
    results_list.append(
        [args.copyprobs_file, phenotype_file, anc, chrom, PRS, number_of_SNPs, skipped_snps, pval, LD_value])


parser = argparse.ArgumentParser()
parser.add_argument("-copyprobs_file",
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
parser.add_argument("-pruned_dir",
                    help="Directory containing pruned.in files",
                    required=True)
parser.add_argument("-chr",
                    help="Chromosome number",
                    required=True)
parser.add_argument("-anc",
                    help="Ancestry being analysed",
                    required=True)
args = parser.parse_args()

# Read in copyprobs
col_names = pd.read_csv(str(args.copyprobs_file), sep=" ", nrows=0).columns
types_dict = {'0': str}
types_dict.update({col: 'int8' for col in col_names if col not in types_dict})
anc_copyprobs = pd.read_csv(str(args.copyprobs_file), sep=" ", dtype=types_dict)
anc_copyprobs.set_index("0", inplace=True)
anc_copyprobs.columns = anc_copyprobs.columns.tolist()[::-1]  # Reverse column names because they are mis-labelled
print("*** Successfully loaded copyprobs file for " + args.anc + " chr" + str(args.chr) + "***")

# Read in idfile
idfile = pd.read_csv(args.idfile, header=None, sep=" ", index_col=0)
# idfile = idfile[~idfile.index.isin(idfile.tail(318).index.tolist())]

# Read in phasefile
positions = anc_copyprobs.columns.tolist()[::-1]
phase = pd.read_csv(str(args.phasefile), header=None, sep=" ", dtype='int8')
phase.columns = positions
phase.index = [val for val in idfile.index.unique().tolist() for _ in (0, 1)]
haps = list()
for h in range(int(phase.shape[0] / 2)):
    haps.extend([1, 2])

results_list = list()
# anc_dict = defaultdict(list)
if not os.path.exists("PRS_calculations"):
    open("PRS_calculations", 'a').close()

# LOOP through GWAS files
phenotypes = pd.read_csv(args.phenotypes, header=None, dtype=str)[0].tolist()
phenotypes = [s.strip() for s in phenotypes]
for file in phenotypes:
    GWAS = pd.read_csv(file, sep='\t', usecols=['beta', 'pval', 'variant'])
    print("***Looking at " + str(file) + "***")
    GWAS[['chr', 'pos', '1', '2']] = GWAS['variant'].str.split(':', expand=True)
    GWAS = GWAS[GWAS['chr'] == str(args.chr)]
    GWAS = GWAS.loc[GWAS['pos'].isin(positions)]  # NB about 80% of painted SNPs are in GWAS file
    for p in [0.01, 0.05, 0.1]:
        GWAS_pval = GWAS[GWAS['pval'] <= p]
        for LD in ["05.merged.filtered.prune.in", "5.merged.filtered.prune.in", "9.merged.filtered.prune.in"]:
            pruned_in = pd.read_csv(str(args.pruned_dir) + "/" + str(args.chr) + "." + str(LD), header=None)
            GWAS_pruned = GWAS_pval.loc[GWAS_pval['pos'].isin(pruned_in[0].astype(str).tolist())]
            list_of_SNPs = GWAS_pruned['pos'].tolist()
            phase_temp = phase[list_of_SNPs]
            anc_copyprobs_temp1 = anc_copyprobs[list_of_SNPs]
            analyse_anc(anc_copyprobs_temp1, str(args.anc), args.chr, p, LD, file)

print("***Success! Now writing results to output file***")
PRS_calculations = pd.DataFrame.from_records(results_list)
# PRS_calculations = pd.DataFrame([anc_dict])
PRS_calculations.to_csv('PRS_calculations', mode='a', header=False, index=False, sep=" ")

# PRS_dict = dict()
# for anc in ["CHG", "EHG", "Farmer", "African", "EastAsian", "WHG", "Yamnaya"]:
#     PRS_per_anc = 0
#     total_number_of_SNPs = 0
#     for chrom, PRS, number_of_SNPs in anc_dict[anc]:
#         PRS_per_anc += PRS
#         total_number_of_SNPs += number_of_SNPs
#     PRS_dict[anc] = PRS_per_anc
#     PRS_dict[anc + "_perSNP"] = PRS_per_anc/total_number_of_SNPs
# PRS_dict['Phenotype'] = args.file_name
# PRS_dict['LD_pruned'] = args.pruned_in
# PRS_dict['p_value'] = args.p_value
# PRS_calculations_temp = pd.DataFrame([PRS_dict])
# PRS_calculations_temp = PRS_calculations_temp[["Phenotype", "LD_pruned", "p_value", "CHG", "EHG", "Farmer",
#                                                "African", "EastAsian", "WHG", "Yamnaya", "CHG_perSNP", "EHG_perSNP",
#                                                "Farmer_perSNP", "African_perSNP", "EastAsian_perSNP", "WHG_perSNP",
#                                                "Yamnaya_perSNP"]]
# PRS_calculations_temp.to_csv('PRS_calculations', mode='a', index=False, header=True, sep=" ")
# print("*** Success! See PRS_calculations for results ***")
