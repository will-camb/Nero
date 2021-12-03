#!/bin/bash
source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "10" ] ; then
    echo "Usage: run_PRS_calculator.sh <phenotype_file> <copyprobs_directory> <phasefile_directory> <idfile> <bootstrap> <reverse_cols> <LD_prune>"
    echo "<phenotype_file>: file containing filenames of phenotype GWAS results, 1 per row eg E4_OBESITY.gwas.imputed_v3.both_sexes.tsv.bgz"
    echo "<copyprobs_directory_imputed>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<copyprobs_directory_original>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<output_files>: output for homozygous inds"
    echo "<phasefile_directory>: Directory of phasefiles; should be named in form chr#.merged.phase"
    echo "<idfile>: Path to idfile: ordered_all_pop_ids_mapped ('individualID popID 1'); only used for phasefile index"
    echo "<variants>: GWAS variants file containing SNP details - chr, pos, rsID, OR/BETA"
    echo "<bootstrap True/False>: Whether to bootstrap or not. Can take value either True or False"
    echo "<reverse_cols True/False>: Whether to reverse cols or not. Can take value either True or False"
    echo "<LD_prune True/False>: Whether to LD prune or not. Can take value either True or False"
    exit 0
fi

phenotype_file="$1"
copyprobs_directory_imputed="$2"
copyprobs_directory_original="$3"
output_files="$4"
phasefile_directory="$5"
idfile="$6"
variants="$7"
bootstrap="$8"
reverse_cols="$9"
LD_prune="${10}"
chrlist=$(seq 1 22)
ancestries="CHG EHG Farmer African EastAsian WHG Yamnaya"

for chr in $chrlist; do
  for anc in $ancestries; do
    echo "python3 PRS_calculator_GWAScat.py -copyprobs_file_original $copyprobs_directory_original/$anc.$chr.master_all_copyprobsperlocus.txt.gz -copyprobs_file_imputed $copyprobs_directory_imputed/temp.$anc.$chr.master_all_copyprobsperlocus.txt.gz -phasefile $phasefile_directory/transformed.$chr.merged.phase.gz -idfile $idfile -variants $variants -phenotypes $phenotype_file -output_files $output_files -chr $chr -anc $anc -bootstrap $bootstrap -reverse_cols $reverse_cols -LD_prune $LD_prune" >> PRS_calculator_commands_temp
  done
done

echo "Now running commands in PRS_calculator_commands in parallel"
shuf PRS_calculator_commands_temp > PRS_calculator_commands_"$phenotype_file"
rm PRS_calculator_commands_temp
cat PRS_calculator_commands_"$phenotype_file" | parallel  -j 50 # Change as appropriate

echo "*** All done, results are in PRS_calculations_v3 ***"