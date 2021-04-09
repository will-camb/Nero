#!/bin/bash
#cd $PBS_O_WORKDIR
#module load anaconda3/4.4.0
#source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "5" ] ; then
    echo "Usage: run_PRS_calculator.sh <file_name> <copyprobs_directory> <phasefile_directory>"
    echo "<phenotype_file>: file containing filenames of phenotype GWAS results, 1 per row eg E4_OBESITY.gwas.imputed_v3.both_sexes.tsv.bgz"
    echo "<copyprobs_file>: Directory of all_copyprobsperlocus files; should be named in form n.master_all_copyprobsperlocus.txt.gz"
    echo "<phasefile>: Directory of phasefiles; should be named in form chr#.merged.phase"
    echo "<idfile>: Path to idfile: ordered_all_pop_ids_mapped ('individualID popID 1')"
    exit 0
fi

phenotype_file="$1"
copyprobs_file="$2"
phasefile="$3"
idfile="$4"
pruned_dir="/willerslev/ukbiobank/LD-pruning/"
chrlist=$(seq 1 22)
ancestries="CHG EHG Farmer African EastAsian WHG Yamnaya"

for chr in $chrlist; do
  for anc in $ancestries; do
    echo "python3 PRS_calculator.py -copyprobs_directory $copyprobs_directory -phasefile_directory $phasefile_directory -idfile $idfile -phenotypes $phenotype_file -pruned_dir $pruned_dir -chr $chr -anc $anc" >> PRS_calculator_commands
  done
done

echo "Now running commands in PRS_calculator_commands in parallel"
cat PRS_calculator_commands | parallel

echo "*** All done, results are in PRS_calculations ***"

#rm PRS_calculator_commands