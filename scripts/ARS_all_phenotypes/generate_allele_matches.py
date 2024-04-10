import pandas as pd

imputed_calls = "/maps/datasets/ukb-AUDIT/imputation_bgen/ukb_mfi_chrALL_v3.txt"
copyprobs_location = "/maps/datasets/ukb-AUDIT/painting_results_aggregate/copyprobs_per_anc/reversed_cols"
pruned_SNPs = "/maps/datasets/ukb-AUDIT/painting_results_aggregate/PRS_calculation/all_phenotypes/intermediate_files/parkinsons/all_clumped_annotated.csv"

# Gte painted sites
painted_sites = []
painted_chr = []
for i in range(1, 23):
    col_names = pd.read_csv(f"{copyprobs_location}/Yamnaya.{i}.master_all_copyprobsperlocus.txt.gz", sep=" ",
                            nrows=0).columns.tolist()
    col_names.pop(0)
    chrs = [i for x in col_names]
    painted_sites = painted_sites + col_names
    painted_chr = painted_chr + chrs
painted = pd.DataFrame(painted_sites, painted_chr, columns=['BP']).reset_index(names='CHR')
painted['BP'] = painted['BP'].astype(int)
painted['CHR'] = painted['CHR'].astype(int)

# Get imputed sites with allele frequency info
imputed_calls = pd.read_csv(imputed_calls, sep="\t", header=None)
imputed_calls = imputed_calls[imputed_calls[3].str.len() == 1]  # Biallelic
imputed_calls = imputed_calls[imputed_calls[4].str.len() == 1]
imputed_calls['CHR'] = imputed_calls[0].str.split(":", expand=True)[0]
imputed_calls = imputed_calls.loc[imputed_calls['CHR'].isin([str(x) for x in range(1, 23)])]
imputed_calls['CHR'] = imputed_calls['CHR'].astype(int)
imputed_calls.rename({2: 'BP'}, axis=1, inplace=True)
imputed_calls['BP'] = imputed_calls['BP'].astype(int)

# Merge to get painted sites with frequencies
imputed_calls_painted = pd.merge(imputed_calls, painted, on=['CHR', 'BP'])

painted_sites = []
for i in range(1, 23):
    col_names = pd.read_csv(f"{copyprobs_location}/Yamnaya.{i}.master_all_copyprobsperlocus.txt.gz", sep=" ",
                            nrows=0).columns.tolist()
    col_names = [x + f"chr{i}" for x in col_names]
    col_names.pop(0)
    painted_sites = painted_sites + col_names
imputed_calls['temp'] = imputed_calls['BP_imputed'].astype(str) + "chr" + imputed_calls['CHR_imputed'].astype(str)
imputed_calls_painted = imputed_calls.loc[imputed_calls['temp'].isin(painted_sites)]


imputed_calls.drop(['temp'], axis=1, inplace=True)





imputed_calls.rename({'BP': 'BP_imputed'}, axis=1, inplace=True)
imputed_calls.rename({'CHR': 'CHR_imputed'}, axis=1, inplace=True)
imputed_calls = imputed_calls.drop_duplicates(subset=['BP_imputed', 'CHR_imputed'])

imputed_calls_painted = pd.merge(imputed_calls, painted, left_on=['CHR_imputed', 'BP_imputed'], right_on=['CHR', 'BP'])

pruned_SNPs = pd.read_csv(pruned_SNPs, sep=" ")




variants['temp'] = variants['BP'].astype(str) + "chr" + variants['CHR'].astype(str)
variants = variants.loc[variants['temp'].isin(painted_sites)]
variants.drop(['temp'], axis=1, inplace=True)
