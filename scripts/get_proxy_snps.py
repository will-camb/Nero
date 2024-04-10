#!/usr/bin/env python

import argparse
import pandas as pd
import subprocess
import os

LD_IDS = "../../GBR-FIN-TSI_ids.txt"
path_sites = "../../ancestral_paths_merged_filtered.sites.gz"
if not os.path.exists("temp"):
    os.makedirs("temp")


def main(args):
    results_list = []
    # read in the SNP data
    snps = pd.read_csv(args.snp_file, sep=" ")
    print(f"Looking at all SNPs for chromosome {args.chrom} in file {args.snp_file}")
    print("**********")
    paths = pd.read_csv(path_sites, sep="\t")
    paths = paths.loc[paths['# [1]CHROM'] == int(args.chrom)]
    # create the filename for the VCF file
    vcf_file = "{}/chr{}.1kg.phase3.v5a.vcf.gz".format(args.genomes_dir, args.chrom)
    # iterate over the SNPs on the chromosome in the snps dataframe
    for index, row in snps[snps["CHR"] == int(args.chrom)].iterrows():
        SNP = row['SNP']
        print("--------")
        print("Examining SNP {}".format(row["SNP"]))
        if row['BP'] in paths['[2]POS'].tolist():
            print(f"SNP {SNP} is in the paths painting, so no need to find proxies")
            R2 = 1
            results_list.append([row['SNP'], row['CHR'], row['BP'], row['other_allele'], row['effect_allele'],
                                 row['effect_allele_frequency'], row['P'], row['OR'], row['SNP'], row['BP'],
                                 row['effect_allele'], row['effect_allele_frequency'], R2])
            continue
        command = ['plink',
                   '--vcf', vcf_file,
                   '--keep', LD_IDS,
                   '--r2',
                   '--ld-snp', SNP,
                   '--ld-window-r2', '0.7',  # Minimum r2 value
                   '--ld-window', '99999',
                   '--ld-window-kb', '1000',
                   '--out', f'temp/output_{SNP}'
                   ]
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"SNP {SNP} not in specified VCF, so can't find proxies. Skipping this SNP. Error message printed "
                  f"below:")
            print(e)
            continue
        output_ld = pd.read_csv(f"temp/output_{SNP}.ld", delim_whitespace=True)
        output_ld = pd.merge(output_ld, paths['[2]POS'], left_on='BP_B', right_on='[2]POS').sort_values("R2",
                                                                                                        ascending=False)
        best_proxy_SNP = output_ld.iloc[0]['SNP_B']
        best_proxy_BP = output_ld.iloc[0]['BP_B']
        R2 = output_ld.iloc[0]['R2']
        with open(f'temp/{best_proxy_SNP}.txt', 'w') as f:
            f.write(f'{best_proxy_SNP}\n')
        command = ['plink',
                   '--vcf', vcf_file,
                   '--keep', LD_IDS,
                   '--extract', f'temp/{best_proxy_SNP}.txt',
                   '--freq',
                   '-out', f'temp/{best_proxy_SNP}_freq']
        subprocess.run(command, check=True)
        freq = pd.read_csv(f'temp/{best_proxy_SNP}_freq.frq', delim_whitespace=True)
        if row['effect_allele_frequency'] <= 0.5:
            if freq.iloc[0]['MAF'] < 0.5:
                proxy_effect_allele = freq.iloc[0]['A1']
                proxy_effect_allele_frequency = freq.iloc[0]['MAF']
            else:
                proxy_effect_allele = freq.iloc[0]['A2']
                proxy_effect_allele_frequency = 1 - freq.iloc[0]['MAF']
        elif row['effect_allele_frequency'] > 0.5:
            if freq.iloc[0]['MAF'] > 0.5:
                proxy_effect_allele = freq.iloc[0]['A1']
                proxy_effect_allele_frequency = freq.iloc[0]['MAF']
            else:
                proxy_effect_allele = freq.iloc[0]['A2']
                proxy_effect_allele_frequency = 1 - freq.iloc[0]['MAF']
        results_list.append([row['SNP'], row['CHR'], row['BP'], row['other_allele'], row['effect_allele'],
                             row['effect_allele_frequency'], row['P'], row['OR'], best_proxy_SNP, best_proxy_BP,
                             proxy_effect_allele, proxy_effect_allele_frequency, R2])
        print(f"Done getting proxy for {SNP}!")
    if not os.path.exists("proxy_snps.out"):
        pd.DataFrame(columns=['SNP', 'CHR', 'BP', 'other_allele', 'effect_allele', 'effect_allele_frequency', 'P',
                                   'OR', 'proxy_SNP', 'proxy_BP', 'proxy_effect_allele',
                                   'proxy_effect_allele_frequency', 'R2']).to_csv("proxy_snps.out", sep=" ", index=False)
    results = pd.DataFrame.from_records(results_list)
    results.to_csv("proxy_snps.out", mode='a', index=False, header=False, sep=" ")
    print(f"Done getting proxies for all SNPs for chromosome {args.chrom}. Saving results at proxy_snps.out")
    print("**********")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--chrom", type=str, required=True, help="Chromosome number")
    parser.add_argument("--snp_file", type=str, required=True, help="Path to SNP file")
    parser.add_argument("--genomes_dir", type=str, required=True, help="Path to 1000 genomes directory")
    args = parser.parse_args()
    main(args)
