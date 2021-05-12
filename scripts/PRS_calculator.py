import pandas as pd
import argparse
import os


def analyse_anc(anc_copyprobs_temp, anc, chrom, pval, LD_value, phenotype_file):
    global list_of_SNPs
    global alt
    global ref
    PRS = 0
    skipped_snps = 0
    for i in list_of_SNPs:
        anc_copyprobs_temp_i = pd.DataFrame(zip(anc_copyprobs_temp.index, haps, anc_copyprobs_temp[str(i)]))
        phase_i = pd.DataFrame(zip(phase_temp.index, haps, phase_temp[str(i)]))
        df = pd.merge(anc_copyprobs_temp_i, phase_i, left_on=[0, 1], right_on=[0, 1]).set_index(0)
        df.columns = ['haps', anc, 'phase']

        #  Check if the minor allele = alt allele

        try:
            minor_allele = GWAS_pruned.loc[GWAS_pruned['pos'] == i].minor_allele.item()
        except ValueError:  # This skips SNPs with more then one entry in the GWAS file
            print("ValueError for minor_allele for SNP at position " + str(i) + " for " + str(anc) + " chr" + str(chrom) + "!")
            skipped_snps += 1
            continue
        try:
            alt_allele = GWAS_pruned.loc[GWAS_pruned['pos'] == i].alt.item()
        except ValueError:  # This skips SNPs with more then one entry in the GWAS file
            print("ValueError for alt_allele for SNP at position " + str(i) + " for " + str(anc) + " chr" + str(chrom) + "!")
            skipped_snps += 1
            continue

        # if GWAS_pruned.loc[GWAS_pruned['pos'] == i].minor_allele.empty or GWAS_pruned.loc[GWAS_pruned['pos'] == i].alt.empty:
        #     print("Can't find minor/alt info for SNP " + str(i) + " for " + str(anc) + str(chrom))
        #     skipped_snps += 1
        #     continue
        if minor_allele == alt_allele:
            #  minor is alt
            #  Find if minor is 1 or 0 in phase
            if df['phase'].sum() > (df.shape[0] / 2):
                alt = 0
                ref = 1
        else:
            #  major is alt
            if df['phase'].sum() > (df.shape[0] / 2):
                alt = 1
                ref = 0
        alt_df = df.loc[df['phase'] == alt]
        if alt_df.empty:
            alt_sum = 0
        else:
            alt_sum = alt_df[anc].sum()
        if isinstance(alt_sum, str):
            print("Minor.dtype is string for " + str(i) + str(anc) + str(chrom))
            skipped_snps += 1
            continue
        ref_df = df.loc[df['phase'] == ref]
        if ref_df.empty:
            ref_sum = 0
        else:
            ref_sum = ref_df[anc].sum()  # Sum painting prob for major allele
        if isinstance(ref_sum, str):
            print("Major.dtype is string for " + str(i) + str(anc) + str(chrom))
            skipped_snps += 1
            continue
        if alt_sum + ref_sum != 0:
            anc_maf = alt_sum / (alt_sum + ref_sum)
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
# positions = anc_copyprobs.columns.tolist()[::-1]  # Is this re-mislabelling? Pretty sure yes
positions = anc_copyprobs.columns.tolist()
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

# Load variants file to find which is alt allele (not always minor!)
variants = pd.read_csv("/willerslev/datasets/UKBiobank/NealeV2/variants.tsv.gz",
                       sep='\t',
                       usecols=['variant', 'chr', 'alt'],
                       dtype={'variant': 'string', 'chr': 'string', 'alt': 'string'})
variants = variants.loc[variants['chr'] == str(args.chr)]
alt = 1
ref = 0
# LOOP through GWAS files
phenotypes = pd.read_csv(args.phenotypes, header=None, dtype=str)[0].tolist()
phenotypes = [s.strip() for s in phenotypes]
for file in phenotypes:
    GWAS = pd.read_csv(file,
                       sep='\t',
                       usecols=['variant', 'beta', 'pval', 'minor_allele'],
                       dtype={'variant': 'string', 'beta': 'float', 'pval': 'float', 'minor_allele': 'string'})
    print("***Looking at " + str(file) + "***")
    GWAS[['chr', 'pos', '1', '2']] = GWAS['variant'].str.split(':', expand=True)
    GWAS = GWAS[GWAS['chr'] == str(args.chr)]
    GWAS = GWAS.loc[GWAS['pos'].isin(positions)]  # NB about 80% of painted SNPs are in GWAS file

    GWAS_variants = pd.merge(variants[['variant', 'alt']], GWAS, on='variant')

    for p in [0.01, 0.05, 0.1]:
        GWAS_pval = GWAS_variants[GWAS_variants['pval'] <= p]
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
