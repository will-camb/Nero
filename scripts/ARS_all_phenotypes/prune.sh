#!/bin/bash
set -o errexit

if [ "$#" -ne "2" ]; then
  echo "Usage: master.sh <phenotype> <sum_stats>"
  echo "<phenotype>: The phenotype being calculated"
  echo "<sum_stats>: Name of sum stats file for the specified phenotype, in folder called 'summary_stats'; NB if Neale Lab should already be combined with file variants.tsv.gz. Cols: SNP CHR BP other_allele effect_allele effect_allele_frequency P OR"
  exit 0
fi

# Define and check user-inputted variables
phenotype="$1"
sum_stats="$2"
CHR_LIST=$(seq 1 22)
VCF_PATH="/maps/projects/lundbeck/data/1000genomes_2015_nature" # Named as "$chr.1000g.freeze9.umich.GRCh37.snps.biallelic.pass.vcf.gz"
LD_IDS="/projects/lundbeck/people/jdn321/GBR-FIN-TSI_ids.txt"
path_sites="/projects/lundbeck/people/jdn321/ancestral_paths_new_sites.tsv.gz"
imputed_calls="/datasets/ukb-AUDIT/painting_results_aggregate/PRS_calculation/all_phenotypes/ukb_mfi_chrALL_v3.txt.gz"
CLUMP_P1=5e-8
CLUMP_R2=0.05
CLUMP_KB=250

mkdir -p intermediate_files || true
cd intermediate_files
mkdir "$phenotype" || true
cp ../summary_stats/"$sum_stats" "$phenotype"
cd "$phenotype"

function calculate_freq() {
  for chr in $CHR_LIST; do plink --silent --vcf "$VCF_PATH"/"$chr".1000g.freeze9.umich.GRCh37.snps.biallelic.pass.vcf.gz --keep $LD_IDS --freq --out "$chr".freq & done
  wait
}

function load_sum_stats() {
  PYTHON_ARG1="$1" PYTHON_ARG2="$2" PYTHON_ARG3="$3" python3 - <<END
import pandas as pd
import os
import sys

# Function to get effect allele frequency
def get_effect_allele_freq(row):
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if row['A1'] == row['effect_allele']:
        return row['MAF']
    elif row['A2'] == row['effect_allele']:
        return 1 - row['MAF']
    elif mapping[row['A1']] == row['effect_allele']:
        return row['MAF']
    else:
        return None

variants = pd.read_csv(os.environ['PYTHON_ARG1'], dtype={'SNP': 'string', 'CHR': 'int8', 'BP': 'int32', 'other_allele': 'string',
                              'effect_allele': 'string', 'effect_allele_frequency': 'float', 'P': 'float', 'OR': 'float'},
                              sep=" ")

df_freq = pd.DataFrame()
for i in range(1,23):
     try:
         df_temp = pd.read_csv(str(i)+".freq.frq", delim_whitespace=True)
         df_freq = pd.concat([df_freq, df_temp])
     except FileNotFoundError:
             continue
df_freq = df_freq.loc[df_freq['MAF']>0] # Just keep SNPs with MAF>0
if df_freq.empty:
  sys.exit()

# Merge the allele frequencies data with the variants data
variants = pd.merge(variants, df_freq, on=['CHR', 'SNP'])

# Calculate the effect allele frequency
variants['effect_allele_frequency'] = variants.apply(get_effect_allele_freq, axis=1)
n_dropped = variants.shape[0] - variants.dropna(subset=["effect_allele_frequency"]).shape[0]
print(f"Dropping {str(n_dropped)} SNPs from summary stats because couldn't find allele frequency information")
variants = variants.dropna(subset=["effect_allele_frequency"]) # Drop SNPs without effect allele frequency information
variants.drop(['A1', 'A2', 'MAF', 'NCHROBS'], axis=1, inplace=True)

# Filter for path sites only
paths = pd.read_csv(os.environ['PYTHON_ARG2'], sep="\t")
paths['temp'] = paths['[2]POS'].astype(str) + "chr" + paths['# [1]CHROM'].astype(str)
variants['temp'] = variants['BP'].astype(str) + "chr" + variants['CHR'].astype(str)
variants = pd.merge(paths['temp'], variants)
variants.drop(['temp'], axis=1, inplace=True)

# Filter for MAF>0/INFO>0.8 in UKB data
# imputed_calls = pd.read_csv(os.environ['PYTHON_ARG3'], sep="\t", header=None, usecols=[1,5])
# imputed_calls = imputed_calls.loc[imputed_calls[7] > 0.8] # INFO
# imputed_calls = imputed_calls.loc[imputed_calls[5] > 0] # MAF
# variants = variants.loc[variants['SNP'].isin(imputed_calls[1].tolist())]

variants.to_csv(os.environ['PYTHON_ARG1'], index=False, sep=" ")
if variants.shape[0] == 0:
  print("No variants left after filtering or restricting to path sites, exiting")
  sys.exit()
for i in range(1,23):
  temp = variants.loc[variants['CHR']==i]
  temp['SNP'].to_csv(str(i)+".rsids", header=False, index=False)
END
wait
}

function convert_and_clump() {
  for chr in $CHR_LIST; do plink --silent --vcf "$VCF_PATH"/"$chr".1000g.freeze9.umich.GRCh37.snps.biallelic.pass.vcf.gz --keep $LD_IDS --extract "$chr".rsids -make-bed --out "$chr".plink & done
  wait
  for chr in $CHR_LIST; do cut -f 2 "$chr".plink.bim | sort | uniq -d >"$chr".dups & done
  wait
  for chr in $CHR_LIST; do plink --silent --bfile "$chr".plink --exclude "$chr".dups --clump "$sum_stats" --clump-p1 $CLUMP_P1 --clump-r2 $CLUMP_R2 --clump-kb $CLUMP_KB --out "$chr".clumped & done
  wait
  rm *.plink.bed &
  rm *.plink.bim &
  rm *.plink.fam &
  rm *.plink.nosex &
  rm *.dups &
  rm *.rsids &
  rm *.frq &
}

function aggregate_results() {
  PYTHON_ARG="$1" python3 - <<END
import pandas as pd
import sys
import os

df_all=pd.DataFrame()
for i in range(1,23):
     try:
             df = pd.read_csv(str(i)+".clumped.clumped", delim_whitespace=True)
             df_all = pd.concat([df_all,df])
     except FileNotFoundError:
             continue
if df_all.empty:
  sys.exit()
df_all.to_csv("all_clumped.csv", index=False, sep=" ")
variants = pd.read_csv(os.environ['PYTHON_ARG'], dtype={'SNP': 'string', 'CHR': 'int8', 'BP': 'int32', 'other_allele': 'string',
                              'effect_allele': 'string', 'effect_allele_frequency': 'float', 'P': 'float', 'OR': 'float'},
                              sep=" ")

# Merge the clumped data with the variants data
merged = pd.merge(variants, df_all, on=['CHR', 'BP', 'SNP'])

try:
  clumped = pd.read_csv("all_clumped.csv", sep=" ")
except pd.errors.EmptyDataError:
  sys.exit()

# Save the data to a CSV file
merged.rename({'P_x':'P'}, axis=1).drop('P_y', axis=1).to_csv("all_clumped_annotated.csv", index=False, sep=" ")
END
  rm *.clumped &
  rm *.clumped.nosex &
  wait
}

# a.	Load sum stats
# b.	Filter for painted sites
echo Calculating allele frequencies
calculate_freq
echo Loading and filtering sum stats for phenotype "$phenotype"
load_sum_stats "$sum_stats" "$path_sites" "$imputed_calls"
echo "Done loading sum stats for phenotype $phenotype"

echo 'Now converting files to plink format, and clumping'
convert_and_clump

echo 'Done clumping, now aggregating'
#   f.	Aggregate results
aggregate_results "$sum_stats"
# if no results to aggregate, exit script
if [[ ! -f all_clumped_annotated.csv ]]; then
  echo 'File all_clumped_annotated.csv does not exist, exiting. Likely due to no SNPs passing p-value threshold. If this is unexpected, check format of your input summary stats file'
  exit 1
fi

cp all_clumped_annotated.csv ../../pruned_files/all_clumped_annotated_"$phenotype".csv
