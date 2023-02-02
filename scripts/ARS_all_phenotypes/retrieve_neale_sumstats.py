import pandas as pd
import numpy as np
import argparse
from pathlib import Path


def get_effect_frequency(minor, alt, minor_AF):
    if minor == alt:
        return minor_AF
    else:
        return 1 - minor_AF


parser = argparse.ArgumentParser()
parser.add_argument("-phenotype",
                    help="The phenotype being retrieved",
                    required=True)
parser.add_argument("-neale_lab_sumstats",
                    help="Path to sumstats file, sep='\t'",
                    required=True)
parser.add_argument("-neale_lab_variants",
                    help="Path to variants file, sep='\t'",
                    required=True)
args = parser.parse_args()

Path("summary_stats/").mkdir(parents=True, exist_ok=True)
print(f"Retrieving GWAS file for {args.phenotype}")
GWAS = pd.read_csv(args.neale_lab_sumstats, sep='\t')
GWAS[['CHR', 'BP', '1', '2']] = GWAS['variant'].str.split(':', expand=True)
print(f"Retrieving variants file")
variants = pd.read_csv(args.neale_lab_variants,
                       sep='\t',
                       usecols=['variant', 'rsid', 'alt', 'ref'],
                       dtype={'variant': 'string', 'rsid': 'string', 'alt': 'string', 'ref': 'string'})
GWAS_variants = pd.merge(variants, GWAS, on='variant')
print("Success! Now getting effect allele frequency, renaming columns, dropping duplicates and X chromosome variants")
GWAS_variants['effect_allele_frequency'] = GWAS_variants.apply(
    lambda row: get_effect_frequency(row['minor_allele'], row['alt'], row['minor_AF']), axis=1)
GWAS_variants['OR'] = np.exp(GWAS_variants['beta'])
GWAS_variants.rename({'rsid': 'SNP', 'alt': 'effect_allele', 'ref': 'other_allele', 'pval': 'P'}, axis=1, inplace=True)
GWAS_variants.drop_duplicates(subset=['SNP'], inplace=True)
GWAS_variants = GWAS_variants.loc[GWAS_variants['CHR'] != 'X']
print(f"Saving file as summary_stats/sumstats_{args.phenotype}.csv.gz")
GWAS_variants.to_csv(f"summary_stats/sumstats_{args.phenotype}.csv.gz", sep=" ", index=False, compression='gzip')
