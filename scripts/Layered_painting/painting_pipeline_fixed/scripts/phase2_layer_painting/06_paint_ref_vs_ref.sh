#!/bin/bash
#
# Step 6: Paint reference vs reference to generate weights/priors
#

set -e

LAYER=$1
CHR=$2

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"

# Parse config
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
SPARSEPAINTER=$(grep "sparsepainter:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
LAYER_DIR="$OUTPUT_BASE/data/layer${LAYER}"

echo "Step 6: Painting reference vs reference for chr$CHR"

DONE_FLAG="$LAYER_DIR/ref_vs_ref/chr${CHR}_refvsref_done.flag"

if [ -f "$DONE_FLAG" ]; then
    echo "  âœ“ Already complete: ref-vs-ref painting for chr$CHR"
    exit 0
fi

# Inputs
REF_PHASE="$LAYER_DIR/ref/chr${CHR}.ref.phase"
GENETIC_MAP="$LAYER_DIR/genetic_maps/chr${CHR}.map"
REF_POPFILE="$BASE_DIR/config/samples/layer1/ref_popfile.txt"
REF_SAMPLES="$BASE_DIR/config/samples/layer1/ref_samples.txt"
LAMBDA=$(cat $LAYER_DIR/lambda/chr${CHR}_lambda.txt)

# Output directory
REFVSREF_DIR="$LAYER_DIR/ref_vs_ref"
mkdir -p $REFVSREF_DIR

# Create namefile (sample IDs, one per line)
NAMEFILE="$REFVSREF_DIR/chr${CHR}_ref_namefile.txt"
cp $REF_SAMPLES $NAMEFILE

echo "  Running SparsePainter (ref vs ref)..."
cd $REFVSREF_DIR

$SPARSEPAINTER \
    -reffile $REF_PHASE \
    -targetfile $REF_PHASE \
    -mapfile $GENETIC_MAP \
    -popfile $REF_POPFILE \
    -namefile $NAMEFILE \
    -haplambda $LAMBDA \
    -prob \
    -chunklength \
    -probstore raw \
    -out chr${CHR}_refvsref

# Gzip outputs
gzip -f chr${CHR}_refvsref_prob.txt 2>/dev/null || true
gzip -f chr${CHR}_refvsref_chunklength.txt 2>/dev/null || true

# Create done flag
touch $DONE_FLAG

echo "  âœ“ Ref-vs-ref painting complete for chr$CHR"
echo "  âœ“ Step 6 complete"
