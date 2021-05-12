#!/bin/bash
#cd $PBS_O_WORKDIR
#module load anaconda3/4.4.0
source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "4" ] ; then
    echo "Usage: run_PRS_calculator.sh <phenotype_file> <copyprobs_directory> <phasefile_directory> <idfile>"
    echo "<phenotype_file>: file containing filenames of phenotype GWAS results, 1 per row eg E4_OBESITY.gwas.imputed_v3.both_sexes.tsv.bgz"
    echo "<copyprobs_directory>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<phasefile_directory>: Directory of phasefiles; should be named in form chr#.merged.phase"
    echo "<idfile>: Path to idfile: ordered_all_pop_ids_mapped ('individualID popID 1')"
    exit 0
fi

phenotype_file="$1"
copyprobs_directory="$2"
phasefile_directory="$3"
idfile="$4"
pruned_dir="/willerslev/ukbiobank/LD-pruning/"
chrlist=$(seq 1 22)
ancestries="CHG EHG Farmer African EastAsian WHG Yamnaya"

for chr in $chrlist; do
  for anc in $ancestries; do
    echo "python3 PRS_calculator.py -copyprobs_file $copyprobs_directory/temp.$anc.$chr.master_all_copyprobsperlocus.txt.gz -phasefile $phasefile_directory/transformed.$chr.merged.phase.gz -idfile $idfile -phenotypes $phenotype_file -pruned_dir $pruned_dir -chr $chr -anc $anc" >> PRS_calculator_commands_temp
  done
done

echo "Now running commands in PRS_calculator_commands in parallel"
shuf PRS_calculator_commands_temp > PRS_calculator_commands
rm PRS_calculator_commands_temp
cat PRS_calculator_commands | parallel  -j 40

echo "*** All done, results are in PRS_calculations ***"

#rm PRS_calculator_commands