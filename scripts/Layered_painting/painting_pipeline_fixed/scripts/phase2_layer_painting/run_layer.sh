#!/bin/bash
#
# Phase 2: Paint specific layer
# Usage: run_layer.sh <layer_number>
#

set -e

LAYER=$1

if [ -z "$LAYER" ]; then
    echo "Error: Layer number not specified"
    echo "Usage: $0 <layer_number>"
    exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"

# Source configuration
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"
LAYER_CONFIG="$BASE_DIR/config/layers/layer${LAYER}_config.yaml"

if [ ! -f "$LAYER_CONFIG" ]; then
    echo "Error: Layer configuration not found: $LAYER_CONFIG"
    exit 1
fi

# Parse config
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
CHROMOSOMES=$(grep "chromosomes:" $CONFIG_FILE | grep -oP '\[\K[0-9, ]+' | tr -d '[]')
SPARSEPAINTER=$(grep "sparsepainter:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
PBWT=$(grep "pbwt:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
GENETIC_MAPS=$(grep "genetic_maps:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
REF_POPFILE="$BASE_DIR/config/samples/layer1/ref_popfile.txt"
BATCH_SIZE=$(grep "target_batch_size:" $CONFIG_FILE | awk '{print $2}')

echo "=================================================="
echo "Phase 2: Painting Layer $LAYER"
echo "=================================================="
echo "Output: $OUTPUT_BASE"
echo "Chromosomes: $CHROMOSOMES"
echo ""

# Load required modules
echo "Loading required modules..."
module load python/3.9.5
module load bcftools/1.20
module load gsl/2.5
module load perl/5.38.0

# Create layer output directories
LAYER_DIR="$OUTPUT_BASE/data/layer${LAYER}"
mkdir -p $LAYER_DIR/{ref,target,genetic_maps,lambda,ref_vs_ref}
mkdir -p $OUTPUT_BASE/results/layer${LAYER}/{raw_probabilities,chunklengths}

# Run each step
for CHR in $(echo $CHROMOSOMES | tr ',' ' '); do
    echo ""
    echo "=================================================="
    echo "Layer $LAYER - Chromosome $CHR"
    echo "=================================================="
    
    # Step 1: Prepare reference and target VCFs
    bash $SCRIPT_DIR/01_prepare_ref_target.sh $LAYER $CHR
    
    # Step 2: Convert VCF to phase format
    bash $SCRIPT_DIR/02_vcf_to_phase.sh $LAYER $CHR
    
    # Step 3: Create genetic maps
    bash $SCRIPT_DIR/03_make_genetic_maps.sh $LAYER $CHR
    
    # Step 4: Split targets into batches
    python3 $SCRIPT_DIR/04_split_target_batches.py $LAYER $CHR $BATCH_SIZE
    
    # Step 5: Estimate lambda (only once, use chr 6 or first chr)
    if [ ! -f "$LAYER_DIR/lambda/chr${CHR}_lambda.txt" ]; then
        bash $SCRIPT_DIR/05_estimate_lambda.sh $LAYER $CHR
    fi
    
    # Step 6: Paint ref vs ref (get weights)
    if [ ! -f "$LAYER_DIR/ref_vs_ref/chr${CHR}_refvsref_done.flag" ]; then
        bash $SCRIPT_DIR/06_paint_ref_vs_ref.sh $LAYER $CHR
    fi
    
    # Step 7: Submit target painting jobs (array job)
    bash $SCRIPT_DIR/07_submit_paint_targets.sh $LAYER $CHR
    
done

echo ""
echo "=================================================="
echo "Layer $LAYER painting jobs submitted"
echo "=================================================="
echo ""
echo "Monitor with: squeue -u \$USER"
echo ""
echo "Once complete, merge results with:"
echo "  bash $SCRIPT_DIR/08_merge_results.sh $LAYER"
