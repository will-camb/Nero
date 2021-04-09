#!/bin/bash
source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "2" ] ; then
    echo "Usage: get_clinvar_haps.sh <rsID> <chr>"
    echo "<rsID>: the SNP of interest"
    echo "<chr>: the chromosome that the SNP is on"
    exit 0
fi

rsID="$1"
chr="$2"

bgenix -g /willerslev/ukbiobank/haplotype_bgen/ukb_hap_chr"$chr"_v2.bgen -incl-rsids "$rsID" -vcf > "$rsID".vcf
bgzip "$rsID".vcf
tabix "$rsID".vcf.gz
pbwt -readVcfGT "$rsID".vcf.gz -writeVcf "$rsID".GT.vcf
bgzip "$rsID".GT.vcf
tabix "$rsID".GT.vcf.gz
pbwt -readVcfGT "$rsID".GT.vcf.gz -writeImputeHapsG "$rsID".haps
impute2chromopainter.pl "$rsID".haps "$rsID".phase
python3 << END
import pandas as pd
phase=pd.read_csv('$rsID.phase', skiprows=3, header=None)
phase[1] = [val for val in pd.read_csv("UKBB_samples", header=None)[0].unique().tolist() for _ in (0, 1)]
haps = list()
for h in range(int(phase.shape[0] / 2)):
    haps.extend([1, 2])
phase[2] = haps
phase = phase[phase[0]==0]
mapping = pd.read_csv("name2id_UKBB", sep=" ", header=None)
phase[1] = phase[1].map(mapping.set_index(0)[2].to_dict())
phase.dropna(inplace=True)
phase[[1,2]].to_csv('output_files/$rsID.haps', header=False, index=False)
END
rm "$rsID"*