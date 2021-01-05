#!/bin/bash
cd $PBS_O_WORKDIR
module load anaconda3/4.4.0

if [ "$#" -ne "5" ] ; then
    echo "Usage: run_PRS_calculator.sh <manifest> <file_name> <copyprobs_directory> <phasefile_directory>"
    echo "<manifest>: the manifest csv file"
    echo "<file_name>: filename from manifest to be calculated"
    echo "<copyprobs_directory>: Directory of all_copyprobsperlocus files; should be named in form n.all_copyprobsperlocus.txt.gz"
    echo "<phasefile_directory>: Directory of phasefiles; should be named in form chr#.merged.phase"
    echo "<idfile>: Path to idfile: ordered_all_pop_ids_mapped ('individualID popID 1')"
    exit 0
fi

manifest="$1"
file_name="$2"
copyprobs_directory="$3"
phasefile_directory="$4"
idfile="$5"
echo "*** Extracting results file from Dropbox ***"
`awk -F ',' -v file_name="$file_name" '$5 == file_name { print $6 }' $manifest`
MY_FILE=$file_name
NEW_EXT=${MY_FILE/bgz/gz}
mv $file_name $NEW_EXT
echo "*** Success! Now running PRS calculator: ***"
echo "python3 PRS_calculator.py -copyprobs_directory $copyprobs_directory -phasefile_directory $phasefile_directory -idfile $idfile -file_name $NEW_EXT"
`python3 PRS_calculator.py -copyprobs_directory $copyprobs_directory -phasefile_directory $phasefile_directory -idfile $idfile -file_name $NEW_EXT`
echo "*** All done, now tidying up temp files ***"
rm $NEW_EXT