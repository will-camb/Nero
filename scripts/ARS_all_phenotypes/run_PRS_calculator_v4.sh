#!/bin/bash
#cd $PBS_O_WORKDIR
#module load anaconda3/4.4.0
source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "5" ] ; then
    echo "Usage: run_PRS_calculator.sh <copyprobs_directory_imputed> <idfile> <output_files> <variants> <bootstrap True/False> <reverse_cols True/False>"
    echo "<copyprobs_directory_imputed>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<output_files>: output for homozygous inds"
    echo "<variants>: pruned variants containing: SNP CHR BP other_allele effect_allele effect_allele_frequency P OR"
    echo "<bootstrap True/False>: Whether to bootstrap or not. Can take value either True or False"
    echo "<reverse_cols True/False>: Whether to reverse cols or not. Can take value either True or False"
    exit 0
fi

copyprobs_directory_imputed="$1"
output_files="$2"
variants="$3"
bootstrap="$4"
reverse_cols="$5"
chrlist=$(seq 1 22)
ancestries="CHG EHG Farmer African EastAsian WHG Yamnaya Steppe"

# For running on a subset of painted individuals in UKB
for chr in $chrlist; do
  for anc in $ancestries; do
    echo "python3 PRS_calculator_v4.py -copyprobs_file_imputed $copyprobs_directory_imputed/temp.$anc.$chr.master_all_copyprobsperlocus.txt.gz -output_files $output_files -variants $variants -chr $chr -anc $anc -bootstrap $bootstrap -reverse_cols $reverse_cols" >> PRS_calculator_commands_temp
  done
done
j=50
echo "Now running commands in PRS_calculator_commands in parallel, across $j cores"
shuf PRS_calculator_commands_temp > PRS_calculator_commands_shuffled
rm PRS_calculator_commands_temp
cat PRS_calculator_commands_shuffled | parallel  -j $j # Change as appropriate

echo "*** All done, results are in PRS_calculations_v4 ***"

rm PRS_calculator_commands_shuffled