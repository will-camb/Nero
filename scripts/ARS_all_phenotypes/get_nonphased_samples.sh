#!/bin/bash
source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "2" ] ; then
    echo "Usage: get_nonphased_samples.sh <rsID> <chr>"
    echo "<rsID>: the SNP of interest"
    echo "<chr>: the chromosome that the SNP is on"
    exit 0
fi

rsID="$1"
chr="$2"

bgenix -g /willerslev/ukbiobank/imputation_bgen/ukb_imp_chr"$chr"_v3.bgen  -incl-rsids "$rsID" -vcf > "$rsID".vcf
bgzip "$rsID".vcf
tabix "$rsID".vcf.gz
pbwt -readVcfGT "$rsID".vcf.gz -writeVcf "$rsID".GT.vcf
bgzip "$rsID".GT.vcf
tabix "$rsID".GT.vcf.gz
bcftools view -g hom "$rsID".GT.vcf.gz -o "$rsID".GT.hom.vcf.gz
pbwt -readVcfGT "$rsID".GT.hom.vcf.gz -writeImputeHapsG "$rsID".haps
impute2chromopainter.pl "$rsID".haps "$rsID".phase
python3 << END
import pandas as pd
from pathlib import Path

Path("output_files/").mkdir(parents=True, exist_ok=True)
phase=pd.read_csv("$rsID.phase", skiprows=3, header=None)
phase[1] = [val for val in pd.read_csv("UKBB_samples", header=None)[0].unique().tolist() for _ in (0, 1)]
haps = list()
for h in range(int(phase.shape[0] / 2)):
    haps.extend([1, 2])
phase[2] = haps

phase_0 = phase[phase[0]==0]
mapping = pd.read_csv("name2id_UKBB", sep=" ", header=None)
phase_0[1] = phase_0[1].map(mapping.set_index(0)[2].to_dict())
samples = phase_0[1].tolist()
vc = pd.Series(samples).value_counts()
hom_samples = vc[vc > 1].index.tolist()
pd.DataFrame(hom_samples).to_csv("output_files/$rsID.hom.0.samples", header=False, index=False)

phase_1 = phase[phase[0]==1]
mapping = pd.read_csv("name2id_UKBB", sep=" ", header=None)
phase_1[1] = phase_1[1].map(mapping.set_index(0)[2].to_dict())
samples = phase_1[1].tolist()
vc = pd.Series(samples).value_counts()
hom_samples = vc[vc > 1].index.tolist()
pd.DataFrame(hom_samples).to_csv("output_files/$rsID.hom.1.samples", header=False, index=False)
END
rm "$rsID"*

