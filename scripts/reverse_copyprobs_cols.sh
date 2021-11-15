#!/bin/bash
source /willerslev/software/venv_python3.6/bin/activate

if [ "$#" -ne "2" ] ; then
    echo "Usage: reverse_copyprobs_cols.sh <copyprobs_file> <copyprobs_file_gzip>"
    echo "<copyprobs_file>: The name of the copyprobs file without .gz"
    echo "<copyprobs_file_gzip>: The name of the copyprobs file with .gz"
    exit 0
fi

copyprobs_file="$1"
copyprobs_file_gzip="$2"

if [ -f "reversed_cols/$copyprobs_file_gzip" ]; then
    echo "File  already exists!"
    exit 0
fi

python3 << END
import pandas as pd
import os

col_df = pd.read_csv("$copyprobs_file_gzip",  sep=" ", nrows=0)
col_names = col_df.columns.tolist()
del col_names[0]
reversed_col_names = col_names[::-1]
reversed_col_names.insert(0, 'ID')
col_df.columns = reversed_col_names
col_df.to_csv("reversed_cols/$copyprobs_file.cols.csv", sep=" ", index=False)
END

gunzip "$copyprobs_file_gzip"
awk '{if (NR!=1) {print}}' "$copyprobs_file" >> reversed_cols/"$copyprobs_file".cols.csv
mv reversed_cols/"$copyprobs_file".cols.csv reversed_cols/"$copyprobs_file"
gzip "$copyprobs_file"
gzip reversed_cols/"$copyprobs_file"