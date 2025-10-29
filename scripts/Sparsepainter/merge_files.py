import gzip
import os
import sys

chrom = '6'
total_files = 400


# Function to merge files
def merge_probfiles(file_type, output_file_name):
    with gzip.open(output_file_name, 'wt') as output_file:
        for i in range(1, total_files + 1):
            print(f"Processing file {i} for {file_type}")
            full_path = f"chr{chrom}/chr{chrom}_target{i}_{file_type}.txt.gz"

            if os.path.isfile(full_path):
                with gzip.open(full_path, 'rt') as input_file:
                    if i != 1:
                        next(input_file, None)
                        next(input_file, None)
                    for line in input_file:
                        output_file.write(line)
            else:
                print(f"File not found: {full_path}")


def merge_chunkfiles(file_type, output_file_name):
    with gzip.open(output_file_name, 'wt') as output_file:
        for i in range(1, total_files + 1):
            print(f"Processing file {i} for {file_type}")
            full_path = f"chr{chrom}/chr{chrom}_target{i}_{file_type}.txt.gz"

            if os.path.isfile(full_path):
                with gzip.open(full_path, 'rt') as input_file:
                    if i != 1:
                        next(input_file, None)

                    for line in input_file:
                        output_file.write(line)
            else:
                print(f"File not found: {full_path}")


# Merge 'prob' files
merge_probfiles('prob', f"chr{chrom}_target_prob.txt.gz")

# Merge 'chunklength' files
merge_chunkfiles('chunklength', f"chr{chrom}_target_chunklength.txt.gz")

# Merge 'chunkcount' files
# merge_chunkfiles('chunkcount', f"chr{chr}_target_chunkcount.txt.gz")
