#!/bin/bash
#
# Step 3: Create genetic maps for the SNP set
#

set -e

LAYER=$1
CHR=$2

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"

# Parse config
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
GENETIC_MAPS_TEMPLATE=$(grep "genetic_maps:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
LAYER_DIR="$OUTPUT_BASE/data/layer${LAYER}"

echo "Step 3: Creating genetic map for chr$CHR"

OUTPUT_MAP="$LAYER_DIR/genetic_maps/chr${CHR}.map"

if [ -f "$OUTPUT_MAP" ]; then
    echo "  âœ“ Already exists: $OUTPUT_MAP"
    exit 0
fi

# Expand template
GENETIC_MAP_FILE="${GENETIC_MAPS_TEMPLATE//\{i\}/$CHR}"

if [ ! -f "$GENETIC_MAP_FILE" ]; then
    echo "  Error: Genetic map file not found: $GENETIC_MAP_FILE"
    exit 1
fi

# Get positions from reference phase file
REF_PHASE="$LAYER_DIR/ref/chr${CHR}.ref.phase"
echo "  Extracting SNP positions..."
tail -n +4 $REF_PHASE | awk '{print $2}' > $LAYER_DIR/genetic_maps/chr${CHR}_snps.txt

# Create map using Python
echo "  Generating genetic map..."
python3 << 'EOF'
import pandas as pd
import numpy as np
import sys
import os

chr_num = int(sys.argv[1])
genetic_map_file = sys.argv[2]
snps_file = sys.argv[3]
output_file = sys.argv[4]

# Load genetic map
print(f"  Loading genetic map: {genetic_map_file}")
gmap = pd.read_csv(genetic_map_file, sep="\t")

# Load SNP positions
print(f"  Loading SNP positions: {snps_file}")
snps = pd.read_csv(snps_file, header=None, names=['Position(bp)'])

# Merge
print(f"  Merging...")
merged = pd.merge(gmap, snps, on='Position(bp)', how='right')

# Fill missing values by interpolation
print(f"  Interpolating missing values...")
for index, row in merged.iterrows():
    if pd.isna(row['Map(cM)']):
        try:
            # Average of SNPs on either side
            merged.at[index, 'Map(cM)'] = (
                merged.iloc[index+1]['Map(cM)'] + merged.iloc[index-1]['Map(cM)']
            ) / 2
        except (IndexError, KeyError):
            try:
                # Use next SNP
                merged.at[index, 'Map(cM)'] = merged.iloc[index+1]['Map(cM)']
            except:
                pass

# Forward fill any remaining
merged.ffill(inplace=True)

# Save in format expected by SparsePainter (position cM)
print(f"  Writing output: {output_file}")
merged[['Position(bp)', 'Map(cM)']].to_csv(output_file, header=False, index=False, sep=" ")

print(f"  âœ“ Created genetic map with {len(merged)} positions")
EOF

python3 -c "$(cat)" $CHR "$GENETIC_MAP_FILE" "$LAYER_DIR/genetic_maps/chr${CHR}_snps.txt" "$OUTPUT_MAP"

# Clean up temp file
rm $LAYER_DIR/genetic_maps/chr${CHR}_snps.txt

echo "  âœ“ Step 3 complete for chr$CHR"
