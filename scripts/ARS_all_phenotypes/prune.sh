#!/bin/bash
set -o errexit

CHR_LIST=$(seq 1 22)
ANC_LIST="CHG WHG EHG Farmer Yamnaya African EastAsian"
VCF_PATH="/willerslev/datasets/1000_Genomes_phase3_v5a/individual_chromosomes" # Named as chr"$chr".1kg.phase3.v5a.vcf.gz
LD_IDS="../../GBR-FIN-TSI_ids.txt"
imputed_calls="../../ukb_mfi_chrALL_v3.txt.gz"
path_sites="../../ancestral_paths_merged_filtered.sites.gz"
CLUMP_P1=5e-8
CLUMP_R2=0.05
CLUMP_KB=250

variants = pd.read_csv(os.environ['PYTHON_ARG1'], dtype={'SNP': 'string', 'CHR': 'int8', 'BP': 'int32', 'other_allele': 'string',
                              'effect_allele': 'string', 'effect_allele_frequency': 'float', 'P': 'float', 'OR': 'float'},
                              sep=" ")

# Filter for path sites only
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
  for i in range(1, 23):
      temp = variants.loc[variants['CHR'] == i]
      temp['SNP'].to_csv(str(i) + ".rsids", header=False, index=False)
  END