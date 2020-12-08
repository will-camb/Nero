#!/bin/bash
cd $PBS_O_WORKDIR
module load anaconda3/4.4.0

if [ "$#" -ne "3" ] ; then
    echo "Usage: run_PRS_calculator.sh <manifest> <copyprobs_directory> <phasefile_directory>"
    echo "<manifest>: the manifest csv file"
    echo "<file_name>: filename from manifest to be calculated"
    echo "<copyprobs_directory>: Directory of all_copyprobsperlocus files; should be named in form n.all_copyprobsperlocus.txt.gz"
    echo "<phasefile_directory>: Directory of phasefiles; should be named in form chr#.merged.phase"
    exit 0
fi

manifest="$1"
file_name="$2"
copyprobs_directory="$3"
phasefile_directory="$4"
export manifest
export file_name
export copyprobs_directory
export phasefile_directory

python3 << END
import pandas as pd
import os
manifest = pd.read_csv("os.environ["manifest"]")
manifest[manifest['UK Biobank Data Showcase Link'] != 'N/A']
wget_command = manifest[manifest['File']==os.environ["file_name"]]['wget command]
END
echo $wget_command
