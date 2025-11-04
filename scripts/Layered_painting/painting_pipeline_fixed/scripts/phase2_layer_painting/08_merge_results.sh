#!/bin/bash
#
# Step 8: Merge painting results from all batches
#

set -e

LAYER=$1

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"

# Parse config
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
CHROMOSOMES=$(grep "chromosomes:" $CONFIG_FILE | grep -oP '\[\K[0-9, ]+' | tr -d '[]')

echo "=================================================="
echo "Merging results for Layer $LAYER"
echo "=================================================="

LAYER_DIR="$OUTPUT_BASE/data/layer${LAYER}"
RESULTS_DIR="$OUTPUT_BASE/results/layer${LAYER}"

for CHR in $(echo $CHROMOSOMES | tr ',' ' '); do
    echo ""
    echo "Processing chromosome $CHR..."
    
    BATCH_DIR="$LAYER_DIR/target/batches/chr${CHR}"
    N_BATCHES=$(cat $BATCH_DIR/n_batches.txt)
    
    # Merge probability files
    echo "  Merging probability files..."
    PROB_OUTPUT="$RESULTS_DIR/raw_probabilities/chr${CHR}_prob.txt.gz"

    if [ -f "$PROB_OUTPUT" ]; then
        echo "    âœ" Already exists: $PROB_OUTPUT"
    else
        mkdir -p $RESULTS_DIR/raw_probabilities

        # Create temporary Python script to merge files
        TEMP_SCRIPT=$(mktemp)
        cat > $TEMP_SCRIPT << 'EOF'
import gzip
import sys
from pathlib import Path

chr_num = sys.argv[1]
batch_dir = Path(sys.argv[2])
output_file = sys.argv[3]
n_batches = int(sys.argv[4])

print(f"    Merging {n_batches} probability files...")

with gzip.open(output_file, 'wt') as output:
    for batch_idx in range(1, n_batches + 1):
        batch_file = batch_dir / f"chr{chr_num}_batch{batch_idx}_prob.txt.gz"

        if not batch_file.exists():
            print(f"      Warning: Missing {batch_file}")
            continue

        with gzip.open(batch_file, 'rt') as input_file:
            # Skip header lines for all but first batch
            if batch_idx != 1:
                next(input_file, None)  # Skip first header line
                next(input_file, None)  # Skip second header line

            for line in input_file:
                output.write(line)

        print(f"      Merged batch {batch_idx}/{n_batches}")

print(f"    âœ" Created {output_file}")
EOF

        python3 $TEMP_SCRIPT $CHR $BATCH_DIR $PROB_OUTPUT $N_BATCHES
        rm $TEMP_SCRIPT
    fi
    
    # Merge chunklength files
    echo "  Merging chunklength files..."
    CHUNK_OUTPUT="$RESULTS_DIR/chunklengths/chr${CHR}_chunklength.txt.gz"

    if [ -f "$CHUNK_OUTPUT" ]; then
        echo "    âœ" Already exists: $CHUNK_OUTPUT"
    else
        mkdir -p $RESULTS_DIR/chunklengths

        # Create temporary Python script to merge files
        TEMP_SCRIPT=$(mktemp)
        cat > $TEMP_SCRIPT << 'EOF'
import gzip
import sys
from pathlib import Path

chr_num = sys.argv[1]
batch_dir = Path(sys.argv[2])
output_file = sys.argv[3]
n_batches = int(sys.argv[4])

print(f"    Merging {n_batches} chunklength files...")

with gzip.open(output_file, 'wt') as output:
    for batch_idx in range(1, n_batches + 1):
        batch_file = batch_dir / f"chr{chr_num}_batch{batch_idx}_chunklength.txt.gz"

        if not batch_file.exists():
            print(f"      Warning: Missing {batch_file}")
            continue

        with gzip.open(batch_file, 'rt') as input_file:
            # Skip header line for all but first batch
            if batch_idx != 1:
                next(input_file, None)

            for line in input_file:
                output.write(line)

        print(f"      Merged batch {batch_idx}/{n_batches}")

print(f"    âœ" Created {output_file}")
EOF

        python3 $TEMP_SCRIPT $CHR $BATCH_DIR $CHUNK_OUTPUT $N_BATCHES
        rm $TEMP_SCRIPT
    fi
    
    echo "  âœ“ Chromosome $CHR complete"
done

echo ""
echo "=================================================="
echo "âœ“ Layer $LAYER results merged successfully"
echo "=================================================="
echo ""
echo "Results located in: $RESULTS_DIR"
echo ""
echo "Next steps:"
echo "  - Review results"
echo "  - Run Phase 3 post-processing (once all layers complete)"
