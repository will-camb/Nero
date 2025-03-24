#!/bin/bash
module load python/3.9.9

input_dir="genetic_data"
output_dir="genetic_maps"

mkdir -p temp
mkdir -p "$output_dir"

for chr in {1..22}
do
 	bcftools query -f '%POS\n' $input_dir/ref.$chr.vcf.gz > temp/sites.$chr
done

export OUTPUT_DIR="$output_dir"

python3 << EOF
import pandas as pd
import numpy as np
import os

output_dir = os.getenv('OUTPUT_DIR')

for i in range(1, 23):
    map = pd.read_csv(f"/projects/lundbeck/data/genetic_map_hg37-hg38/genetic_map_GRCh37_chr{i}.txt", sep="\t")
    sites = pd.read_csv(f"temp/sites.{i}", header=None)
    map = pd.merge(map, sites, left_on='Position(bp)', right_on=0, how='right')
    for index, row in map.iterrows():
        if np.isnan(row['Map(cM)']):  # Fill in missing values with average of SNP on either side
            try:
                map.at[index, 'Map(cM)'] = (map.iloc[index+1]['Map(cM)'] + map.iloc[index-1]['Map(cM)']) / 2
            except IndexError:
                try:
                    map.at[index, 'Map(cM)'] = map.iloc[index+1]['Map(cM)']
                except:
                	continue
    map.ffill(inplace=True)  # forward fill any missing values
    map['Map(M)'] = map['Map(cM)'] / 100
    map[['Position(bp)', 'Map(cM)']].to_csv(f"{output_dir}/{i}.map", header=False, index=False, sep=" ")
EOF

rm -r temp
