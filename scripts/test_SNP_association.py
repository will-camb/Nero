import argparse
import os
from statistics import stdev, mean
import sys

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import t
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

# disable pandas chained assignment warning
pd.options.mode.chained_assignment = None


# Define a function to calculate odds ratios and confidence intervals
def calc_odds_ratio(beta, ci):
    odds_ratio = np.exp(beta)
    lower_ci = np.exp(beta - 1.96 * ci)
    upper_ci = np.exp(beta + 1.96 * ci)
    return odds_ratio, lower_ci, upper_ci


# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-chr", help="Chromosome number", required=True)
parser.add_argument("-anc", help="Ancestry being analyzed", required=True)
parser.add_argument("-position", help="Position in hg19", required=True)
args = parser.parse_args()

# Load copyprobs data
copyprobs_file = f"/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/reversed_cols/parquet_files/{args.anc}.{args.chr}.master_all_copyprobsperlocus.test.parquet"
try:
    copyprobs = pd.read_parquet(copyprobs_file, columns=['ID', args.position])
except (Exception,):
    print(f"Couldn't find copyprobs column corresponding to chr{args.chr}:{args.position}, so skipping")
    sys.exit()
# anc_copyprobs = pd.read_csv(
#     f"/willerslev/ukbiobank/painting_results_aggregate/copyprobs_per_anc/reversed_cols/{args.anc}.{args.chr}.master_all_copyprobsperlocus.txt.gz",
#     sep=" ", usecols=['ID', args.position], dtype={'ID': str, args.position: 'int8'})

# Create haps column
copyprobs['haps'] = np.tile([1, 2], len(copyprobs) // 2)

# Create ID-hap column
copyprobs['ID-hap'] = copyprobs['ID'] + '_' + copyprobs['haps'].astype(str)

# Load phase data
phase_file = f"/willerslev/ukbiobank/painting_results_aggregate/phasefiles_transformed/{args.chr}.merged.positions.phase.gz"
phase = pd.read_csv(phase_file, sep=" ", dtype='int8', usecols=[args.position]).squeeze("columns")
# Load ID file
id_file = "/willerslev/ukbiobank/painting_results_aggregate/ordered_all_pop_ids_mapped_modern"
id_df = pd.read_csv(id_file, sep=" ", header=None, index_col=0)

# Create haps column
phase_haps = np.tile([1,2], len(phase) //2)

# Create ID-hap column
phase = pd.DataFrame({
    'ID': [val for val in id_df.index.unique().tolist() for _ in (0, 1)],
    args.position: phase.tolist(),
    'haps': phase_haps
})
phase['ID-hap'] = phase['ID'] + '_' + phase['haps'].astype(str)

# Filter out rows with low copyprobs
copyprobs = copyprobs[copyprobs[args.position] >= 5]

# Filter phase data to only keep rows with high copyprobs
phase = phase[phase['ID-hap'].isin(copyprobs['ID-hap'].tolist())]

# Check if there is variation at selected sites
if len(phase[args.position].unique()) == 1:  # No variation so can't test association
    table = pd.DataFrame({
        "Feature": args.position,
        "Coefficient": "NA",
        "Beta Value": "NA",
        "Standard Error": "NA",
        "t-value": "NA",
        "p-value": "NA",
        "Odds Ratio": "NA",
        "95% CI Lower": "NA",
        "95% CI Upper": "NA",
    }, index=[0])
    table['Ancestry'] = args.anc
    table['Chromosome'] = args.chr
    output_path = 'SNP_associations_results_monomorphic.csv'
    table.head(1).to_csv(output_path, mode='a', header=not os.path.exists(output_path), index=False)
    print("Site is monomorphic at haplotypes passing the painting probability for this ancestry. Exiting as cannot "
          "calculate association")
    sys.exit()

height = pd.read_csv("/willerslev/ukbiobank/ukb42425.subsets/ukb42425.standingheight.csv")
height['eid'] = height['eid'].astype(str)
pcs = pd.read_csv("/willerslev/ukbiobank/ukb42425.subsets/ukb42425.pcs.csv.gz")
pcs['eid'] = pcs['eid'].astype(str)
age_sex = pd.read_csv("/willerslev/ukbiobank/ukb42425.subsets/ukb42425.age-sex.csv.gz")
age_sex['eid'] = age_sex['eid'].astype(str)
name2id = pd.read_csv("/willerslev/ukbiobank/name2id_UKBB", sep=" ", dtype={0: str}, header=None)

merged = pd.merge(phase, name2id, left_on='ID', right_on=2)
merged = pd.merge(merged, height, left_on=0, right_on='eid')
merged = pd.merge(merged, pcs, left_on=0, right_on='eid')
merged = pd.merge(merged, age_sex, left_on=0, right_on='eid')
merged = merged[
    ['50-0.0', args.position, '21003-0.0', '31-0.0', '22009-0.1', '22009-0.2', '22009-0.3', '22009-0.4', '22009-0.5',
     '22009-0.6', '22009-0.7', '22009-0.8', '22009-0.9', '22009-0.10', '22009-0.11', '22009-0.12', '22009-0.13',
     '22009-0.14', '22009-0.15', '22009-0.16', '22009-0.17', '22009-0.18', '22009-0.19', '22009-0.20']]
merged = merged.dropna()


# SKLEARN LINEAR REGRESSION
# Define the predictor variables (genetic variant, age, sex, and PCs 1-20)
X = merged[
    [args.position, '21003-0.0', '31-0.0', '22009-0.1', '22009-0.2', '22009-0.3', '22009-0.4', '22009-0.5', '22009-0.6',
     '22009-0.7', '22009-0.8', '22009-0.9', '22009-0.10', '22009-0.11', '22009-0.12', '22009-0.13', '22009-0.14',
     '22009-0.15', '22009-0.16', '22009-0.17', '22009-0.18', '22009-0.19', '22009-0.20']]

for col in X.columns.tolist():
    X[col] = X[col].astype('float')

scaler = StandardScaler()
X = scaler.fit_transform(X)

# Define the outcome variable (height)
y = merged['50-0.0']

# Fit the linear regression model
model = LinearRegression().fit(X, y)

# Get the coefficients and standard errors for each feature
coefficients = model.coef_
se = np.sqrt(np.diag(np.linalg.inv(np.dot(X.T, X))))

# Calculate the t-values and p-values for each feature
n_features = X.shape[1]
t_values = coefficients / se
p_values = 2 * t.sf(np.abs(t_values), len(y) - n_features - 1)

# Calculate the 95% confidence intervals for each coefficient
ci = 1.96 * se

# Calculate the beta values
beta_values = coefficients * np.std(X, axis=0) / np.std(y)

# Convert the coefficients to odds ratios and calculate confidence intervals
odds_ratios, lower_ci, upper_ci = calc_odds_ratio(beta_values, ci)

# Create a table of coefficients, beta values, p-values, odds ratios, and confidence intervals
table = pd.DataFrame({
    "Feature": merged[
        [args.position, '21003-0.0', '31-0.0', '22009-0.1', '22009-0.2', '22009-0.3', '22009-0.4', '22009-0.5',
         '22009-0.6',
         '22009-0.7', '22009-0.8', '22009-0.9', '22009-0.10', '22009-0.11', '22009-0.12', '22009-0.13', '22009-0.14',
         '22009-0.15', '22009-0.16', '22009-0.17', '22009-0.18', '22009-0.19', '22009-0.20']].columns,
    "Coefficient": coefficients,
    "Beta Value": beta_values,
    "Standard Error": se,
    "t-value": t_values,
    "p-value": p_values,
    "Odds Ratio": odds_ratios,
    "95% CI Lower": lower_ci,
    "95% CI Upper": upper_ci
})
table['Ancestry'] = args.anc
table['Chromosome'] = args.chr
# Print the table
output_path = 'SNP_associations_results_sklearn.csv'
table.head(1).to_csv(output_path, mode='a', header=not os.path.exists(output_path), index=False)

# STATSMODEL LINEAR REGRESSION
model = sm.OLS(y, sm.add_constant(X)).fit()

table = pd.DataFrame({
    'Feature': args.position,
    'beta': model.params,
    'se': model.bse,
    't': model.tvalues,
    'pval': model.pvalues,
    'CI low': model.conf_int(alpha=0.05, cols=None)[0],
    'CI high': model.conf_int(alpha=0.05, cols=None)[1],
    'Ancestry': args.anc,
    'Chromosome': args.chr
    })

# Print the table
output_path = 'SNP_associations_results_statsmodel.csv'
table.iloc[[1]].to_csv(output_path, mode='a', header=not os.path.exists(output_path), index=False)