#!/bin/bash

# Input and output files
input_file="chr6_target_prob.txt.gz"
output_file="chr6_target_prob_archaic.txt"

# Unzip, process, and rezip
gzip -dc "$input_file" | awk -v OFS="\t" '
BEGIN { header_skipped = 0 }
NR == 1 {
    # Skip the first row (metadata or other non-header content)
    header_skipped = 1
    next
}
NR == 2 {
    # Process the second row as the header
    printf $1  # Print the index column name
    for (i = 2; i <= NF; i++) {
        printf OFS $i "_hap1" OFS $i "_hap2"
    }
    print ""
    next
}
{
    # Process the data rows
    row = $1
    for (i = 2; i <= NF; i++) {
        split($i, hap, "|")
        split(hap[1], archaic1, ",")
        split(hap[2], archaic2, ",")
        row = row OFS archaic1[2] OFS archaic2[2]
    }
    print row
}' > "$output_file"

# Compress the output back if needed
gzip "$output_file"

echo "Processing complete. Output saved to $output_file.gz"