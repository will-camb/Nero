#!/bin/bash
#
# Step 2: Convert VCF to SparsePainter phase format
#

set -e

LAYER=$1
CHR=$2

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"

# Parse config
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
PBWT=$(grep "pbwt:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
LAYER_DIR="$OUTPUT_BASE/data/layer${LAYER}"

echo "Step 2: Converting VCF to phase format for chr$CHR"

REF_VCF="$LAYER_DIR/ref/chr${CHR}.ref.vcf.gz"
TARGET_VCF="$LAYER_DIR/target/chr${CHR}.target.vcf.gz"

REF_PHASE="$LAYER_DIR/ref/chr${CHR}.ref.phase"
TARGET_PHASE="$LAYER_DIR/target/chr${CHR}.target.phase"

# Convert reference VCF
if [ ! -f "$REF_PHASE" ]; then
    echo "  Converting reference VCF to phase format..."
    $PBWT -readVcfGT $REF_VCF -writePhase $REF_PHASE
    echo "  âœ“ Created $REF_PHASE"
else
    echo "  âœ“ Already exists: $REF_PHASE"
fi

# Convert target VCF
if [ ! -f "$TARGET_PHASE" ]; then
    echo "  Converting target VCF to phase format..."
    $PBWT -readVcfGT $TARGET_VCF -writePhase $TARGET_PHASE
    echo "  âœ“ Created $TARGET_PHASE"
else
    echo "  âœ“ Already exists: $TARGET_PHASE"
fi

echo "  âœ“ Step 2 complete for chr$CHR"
