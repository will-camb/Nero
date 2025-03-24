#!/bin/bash
set -o errexit

# Run separately for each phenotype
# Steps:
# 1.	Filter and LD-prune markers (optional)
# 2.	Impute ancestry for missing SNPs
# 3.	Get homozygous individuals for effect allele
# 4.	Run PRS calculator
# 5.	Run accelerated CI

# Requirements: in the master folder you need the following
# GBR-FIN-TSI_ids.txt containing a list of 1KG individuals to use for LD pruning
# impute_ancestry_bp.py
# run_PRS_calculator_v4.sh
# PRS_calculator_v4.py
# get_nonphased_samples.sh
# sum_Steppe.py
# name2id_UKBB
# UKBB_samples
# ukb_mfi_chrALL_v3.txt.gz - file containing INFO scores for all UKB imputed genotypes (https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1965) concatenated into one file

if [ "$#" -ne "7" ]; then
  echo "Usage: master.sh <phenotype> <sum_stats> <copyprobs_location> <prune> <restrict_to_painted_sites>"
  echo "<phenotype>: The phenotype being calculated"
  echo "<sum_stats>: Name of sum stats file for the specified phenotype, in folder called 'summary_stats'; NB if Neale Lab should already be combined with file variants.tsv.gz. Cols: SNP CHR BP other_allele effect_allele effect_allele_frequency P OR"
  echo "<copyprobs_location>: e.g. /willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/reversed_cols/"
  echo "<prune>: Whether to LD-prune the summary stats using plink; must take values either 'true' or 'false'"
  echo "<restrict_to_painted_sites>: Whether to restrict to painted sites; this option speeds up compute time significantly as no need to impute painting scores, but reduces the number of SNPs available for computation. Must take values 'true' or 'false'. NB Only runs if also LD-pruning!"
  echo "<restrict_to_path_sites>: Whether to restrict to sites included in the paths analysis (Pearson et al 2023). Must take values 'true' or 'false'. NB Only runs if also LD-pruning!"
  echo "<filter_maf_info>: Whether to filter for MAF (>0.01)and INFO (>0.8) in UKB imputed data; this is a good idea if pruning, but may elimitae fine-mapped SNPs. Must take values 'true' or 'false'"
  exit 0
fi

# Define and check user-inputted variables
phenotype="$1"
sum_stats="$2"
copyprobs_location="$3"
prune="$4"
restrict_to_painted_sites="$5"
restrict_to_path_sites="$6"
filter_maf_info="$7"
if [ "$prune" != true ] && [ "$prune" != false ]; then
  echo variable prune must take value true or false! Terminating
  exit 1
fi
if [ "$restrict_to_painted_sites" != true ] && [ "$restrict_to_painted_sites" != false ]; then
  echo variable restrict_to_painted_sites must take value true or false! Terminating
  exit 1
fi
if [ "$restrict_to_path_sites" != true ] && [ "$restrict_to_path_sites" != false ]; then
  echo variable restrict_to_path_sites must take value true or false! Terminating
  exit 1
fi
if [ "$filter_maf_info" != true ] && [ "$filter_maf_info" != false ]; then
  echo variable filter_maf_info must take value true or false! Terminating
  exit 1
fi

# Define other constants
CHR_LIST=$(seq 1 22)
ANC_LIST="CHG WHG EHG Farmer Yamnaya African EastAsian"
VCF_PATH="/willerslev/datasets/1000_Genomes_phase3_v5a/individual_chromosomes" # Named as chr"$chr".1kg.phase3.v5a.vcf.gz
LD_IDS="/datasets/ukb-AUDIT/painting_results_aggregate/PRS_calculation/all_phenotypes/GBR-FIN-TSI_ids.txt"
imputed_calls="/datasets/ukb-AUDIT/painting_results_aggregate/PRS_calculation/all_phenotypes/ukb_mfi_chrALL_v3.txt.gz"
#/maps/datasets/ukb-AUDIT/imputation_bgen/ukb_mfi_chrALL_v3.txt
path_sites="/datasets/ukb-AUDIT/painting_results_aggregate/PRS_calculation/all_phenotypes/ancestral_paths_merged_filtered.sites.gz"
CLUMP_P1=5e-8
CLUMP_R2=0.05
CLUMP_KB=250

mkdir -p intermediate_files || true
cd intermediate_files
mkdir "$phenotype" || true
cp ../summary_stats/"$sum_stats" "$phenotype"
cd "$phenotype"

function load_sum_stats() {
  PYTHON_ARG1="$1" PYTHON_ARG2="$2" PYTHON_ARG3="$3" PYTHON_ARG4="$4" PYTHON_ARG5="$5" PYTHON_ARG6="$6" PYTHON_ARG7="$7" python3 - <<END
import pandas as pd
import os
import sys

variants = pd.read_csv(os.environ['PYTHON_ARG1'], dtype={'SNP': 'string', 'CHR': 'int8', 'BP': 'int32', 'other_allele': 'string',
                              'effect_allele': 'string', 'effect_allele_frequency': 'float', 'P': 'float', 'OR': 'float'},
                              sep=" ")
# OPTIONAL filter for MAF and INFO in UKB data
if os.environ['PYTHON_ARG4'] == 'true':
  imputed_calls = pd.read_csv(os.environ['PYTHON_ARG5'], sep="\t", header=None)
  imputed_calls = imputed_calls.loc[imputed_calls[7] > 0.8] # INFO
  imputed_calls = imputed_calls.loc[imputed_calls[5] > 0.01] # MAF
  variants = variants.loc[variants['SNP'].isin(imputed_calls[1].tolist())]

# OPTIONAL filter for painted sites only
if os.environ['PYTHON_ARG3'] == 'true':
  painted_sites = []
  for i in range(1,23):
    col_names = pd.read_csv(f"{os.environ['PYTHON_ARG2']}/Yamnaya.{i}.master_all_copyprobsperlocus.txt.gz", sep=" ", nrows=0).columns.tolist()
    col_names = [x + f"chr{i}" for x in col_names]
    col_names.pop(0)
    painted_sites = painted_sites + col_names
  variants['temp'] = variants['BP'].astype(str) + "chr" + variants['CHR'].astype(str)
  variants = variants.loc[variants['temp'].isin(painted_sites)]
  variants.drop(['temp'], axis=1, inplace=True)

# OPTIONAL filter for path sites only
if os.environ['PYTHON_ARG6'] == 'true':
  paths = pd.read_csv(os.environ['PYTHON_ARG7'], sep="\t")
  paths['temp'] = paths['[2]POS'].astype(str) + "chr" + paths['# [1]CHROM'].astype(str)
  variants['temp'] = variants['BP'].astype(str) + "chr" + variants['CHR'].astype(str)
  variants = pd.merge(paths['temp'], variants)
  variants.drop(['temp'], axis=1, inplace=True)

variants.to_csv(os.environ['PYTHON_ARG1'], index=False, sep=" ")
if variants.shape[0] == 0:
  print("No variants left after filtering or restricting to painted/path sites, exiting")
  sys.exit()
for i in range(1,23):
  temp = variants.loc[variants['CHR']==i]
  temp['SNP'].to_csv(str(i)+".rsids", header=False, index=False)
END
}

function run_ld_pruning() {
  for chr in $CHR_LIST; do plink --vcf "$VCF_PATH"/chr"$chr".1kg.phase3.v5a.vcf.gz --extract "$chr".rsids --keep $LD_IDS -make-bed --out "$chr".plink & done
  wait
  for chr in $CHR_LIST; do cut -f 2 "$chr".plink.bim | sort | uniq -d >"$chr".dups & done
  wait
  for chr in $CHR_LIST; do plink --bfile "$chr".plink --exclude "$chr".dups --clump "$sum_stats" --clump-p1 $CLUMP_P1 --clump-r2 $CLUMP_R2 --clump-kb $CLUMP_KB --out "$chr".clumped & done
  wait
  for chr in $CHR_LIST; do plink --bfile "$chr".plink --freq --out "$chr".freq & done
  wait
  rm *.plink.* &
  rm *.dups &
  rm *.rsids &
}

function aggregate_results() {
  PYTHON_ARG="$1" python3 - <<END
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

df_all=pd.DataFrame()
for i in range(1,23):
     try:
             df = pd.read_csv(str(i)+".clumped.clumped", delim_whitespace=True)
             df_freq = pd.read_csv(str(i)+".freq.frq", delim_whitespace=True)
             df = df.merge(df_freq, how='left', left_on=['SNP','CHR'], right_on=['SNP','CHR'])
             df_all = pd.concat([df_all,df])
     except FileNotFoundError:
             continue
if df_all.empty:
  sys.exit()
df_all.to_csv("all_clumped.csv", index=False, sep=" ")
variants = pd.read_csv(os.environ['PYTHON_ARG'], dtype={'SNP': 'string', 'CHR': 'int8', 'BP': 'int32', 'other_allele': 'string',
                              'effect_allele': 'string', 'effect_allele_frequency': 'float', 'P': 'float', 'OR': 'float'},
                              sep=" ")

# Merge the allele frequencies data with the variants data
merged = pd.merge(variants, df_all, on=['CHR', 'BP', 'SNP'])

# Calculate the effect allele frequency
merged['effect_allele_frequency'] = merged.apply(get_effect_allele_freq, axis=1)
n_dropped = merged.shape[0] - merged.dropna(subset=["effect_allele_frequency"]).shape[0]
print(f"Dropping {str(n_dropped)} SNPs because couldn't find allele frequency information")
merged = merged.dropna(subset=["effect_allele_frequency"]) # Drop SNPs without effect allele frequency information

try:
  clumped = pd.read_csv("all_clumped.csv", sep=" ")
except pd.errors.EmptyDataError:
  sys.exit()

# Save the data to a CSV file
merged.rename({'P_x':'P'}, axis=1).drop('P_y', axis=1).to_csv("all_clumped_annotated.csv", index=False, sep=" ")
END
rm *.frq
}

function make_impute_sites_commands() {
  PYTHON_ARG="$1" python3 <<END
import pandas as pd
import os

df = pd.read_csv("all_clumped_annotated.csv", sep=" ")
with open("impute_commands", "w") as myfile:
    for index,row in df.iterrows():
        for anc in ['CHG', 'WHG','EHG','Farmer','Yamnaya','African','EastAsian']:
          myfile.write(f"python3 ../../impute_ancestry_bp.py -copyprobs_file {os.environ['PYTHON_ARG']}/{anc}.{str(row['CHR'])}.master_all_copyprobsperlocus.txt.gz -bp {str(row['BP'])} -chr {str(row['CHR'])} -anc {anc} -reverse_cols False\n")
END
}

function concat_impute_results() {
  python3 <<END
import pandas as pd

variants = pd.read_csv("all_clumped_annotated.csv", sep=" ")
random_variant = variants.iloc[0]['BP'].item()
for anc in ['CHG', 'WHG','EHG','Farmer','Yamnaya','African','EastAsian']:
    for i in range(1,23):
        new_df = pd.DataFrame(index=pd.read_csv(f"imputed/temp.{random_variant}.Yamnaya.master_all_copyprobsperlocus.txt.gz", index_col="ID").index)
        variants_i = variants.loc[variants['CHR']==i]
        sites = variants_i['BP'].tolist()
        for site in sites:
            temp = pd.read_csv(f"imputed/temp.{site}.{anc}.master_all_copyprobsperlocus.txt.gz", index_col="ID")
            new_df = pd.concat([new_df, temp], axis=1)
        new_df.to_csv(f"imputed/temp.{anc}.{i}.master_all_copyprobsperlocus.txt.gz")
END
}

function make_get_nonphased_samples_commands() {
  python3 <<END
import pandas as pd

variants = pd.read_csv("../all_clumped_annotated.csv", sep=" ")
mylist = []
for index,row in variants.iterrows():
  mylist.append(f"bash get_nonphased_samples.sh {row['SNP']} {str(row['CHR'])}")
#variants['temp'] = variants[f"bash get_nonphased_samples.sh {str(variants['SNP'])} {str(variants['CHR'])}"]
variants['temp'] = mylist
variants['temp'].to_csv("get_nonphased_samples_commands.sh", index=False, header=False)
END
}

function make_unimputed_copyprobs() {
 PYTHON_ARG1="$1" PYTHON_ARG2="$2" PYTHON_ARG3="$3" python3 <<END
import pandas as pd
import os
import sys

variants = pd.read_csv("all_clumped_annotated.csv", sep=" ")
variants_i = variants.loc[variants['CHR']==int(os.environ['PYTHON_ARG2'])]
if variants_i.empty:
  sys.exit()
sites = ["ID"] + variants_i['BP'].astype(str).tolist()
types_dict = {'ID': 'str'}
types_dict.update({snp: 'int8' for snp in sites if snp not in types_dict})
copyprobs = pd.read_csv(f"{os.environ['PYTHON_ARG1']}/{os.environ['PYTHON_ARG3']}.{os.environ['PYTHON_ARG2']}.master_all_copyprobsperlocus.txt.gz", sep=" ", usecols=sites, dtype=types_dict)
copyprobs.to_csv(f"imputed/temp.{os.environ['PYTHON_ARG3']}.{os.environ['PYTHON_ARG2']}.master_all_copyprobsperlocus.txt.gz", index=False)
END
}

# Step 1: LD-prune markers
# a.	Load sum stats
# b.	Filter for MAF, INFO, OPTIONAL painted sites
# c.	Filter 1KG pops for these
# d.	Remove duplicates
# e.	Run clumping method
# f.	Aggregate results

# a.	Load sum stats
# b.	Filter for MAF, INFO, OPTIONAL painted sites
echo Loading sum stats for phenotype "$phenotype"
if [ "$filter_maf_info" = true ]; then echo Filtering for MAF and INFO in UKB data ; fi
if [ "$filter_maf_info" = false ]; then echo You have specified no filtering for MAF and INFO in UKB data ; fi
load_sum_stats "$sum_stats" "$copyprobs_location" "$restrict_to_painted_sites" "$filter_maf_info" "$imputed_calls" "$restrict_to_path_sites" "$path_sites"
echo **Done loading sum stats for phenotype "$phenotype"**

if [ "$prune" = true ]; then
  # c.	Filter 1KG pops for these
  # d.	Remove duplicates
  # e.	Run clumping method
  run_ld_pruning
  echo '**Done clumping, now aggregating**'
#   f.	Aggregate results
  aggregate_results "$sum_stats"
  # if no results to aggregate, exit script
  if [[ ! -f all_clumped_annotated.csv ]] ; then
  echo 'File all_clumped_annotated.csv does not exist, exiting. Likely due to no SNPs passing p-value threshold. If this is unexpected, check format of your input summary stats file'
  exit 1
  fi
fi

if [ "$prune" = false ]; then
  cp "$sum_stats" all_clumped.csv
  cp "$sum_stats" all_clumped_annotated.csv
  echo Not running LD pruning. Saving original sumstats file "$sum_stats" as all_clumped.csv and all_clumped_annotated.csv
fi
echo '**Done aggregating, now imputing ancestry if required (i.e. if restrict_to_painted_sites is false)**'

# Exit if number of SNPs is below threshold:
num_snps=$(tail -n +2 all_clumped_annotated.csv | wc -l)
echo "Total number of SNPs used: $num_snps"
if [ "$num_snps" -le 5 ]; then
  echo 'Not enough SNPs in all_clumped_annotated.csv (<=5), probably because not enough pass genome-wide significance. Exiting script.'
  cd ../
  rm -r "$phenotype"
  exit 1
fi

# 2.  Impute ancestry for missing SNPs
# a.  Impute sites separately
# b.  Concatenate results to master copyprobs file
if [ "$restrict_to_painted_sites" = false ]; then
  # a.  Impute sites separately
  touch impute_commands
  make_impute_sites_commands "$copyprobs_location"
  cat impute_commands | parallel -j 30
  # b.  Concatenate results to master copyprobs file
  concat_impute_results
fi

if [ "$restrict_to_painted_sites" = true ]; then
  mkdir -p imputed || true
  export -f make_unimputed_copyprobs
  for chr in $CHR_LIST; do for anc in $ANC_LIST; do
  echo make_unimputed_copyprobs "$copyprobs_location" "$chr" "$anc" >> make_unimputed_copyprobs_commands; done; done
  cat make_unimputed_copyprobs_commands | parallel  -j 30
  rm make_unimputed_copyprobs_commands
fi
echo '**Done imputation if specified, now calculating homozygous individuals**'

# 3.	Get homozygous individuals
mkdir non_phased_snps || true
cp ../../name2id_UKBB non_phased_snps
cp ../../UKBB_samples non_phased_snps
cp ../../get_nonphased_samples.sh non_phased_snps
cd non_phased_snps
make_get_nonphased_samples_commands
cat get_nonphased_samples_commands.sh | parallel -j 30
cd ../
echo '**Done! Now running PRS calculator**'

# 4.	Run PRS calculator
cp ../../PRS_calculator_v4.py PRS_calculator_v4.py
bash ../../run_PRS_calculator_v4.sh imputed non_phased_snps/output_files/ all_clumped_annotated.csv False False &>stdout.file
mv PRS_calculations_v4 PRS_calculations_v4_"$phenotype"
mkdir ../../results_files || true
cp PRS_calculations_v4_"$phenotype" ../../results_files/PRS_calculations_v4_"$phenotype"
echo '*** Job complete! Results are in PRS_calculations_v4_"$phenotype" ***'
