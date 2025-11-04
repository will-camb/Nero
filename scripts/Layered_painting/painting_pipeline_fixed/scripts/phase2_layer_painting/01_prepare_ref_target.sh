#!/bin/bash
#
# Step 1: Prepare reference and target VCFs
# Extracts reference samples from ancient data and target samples from UKB
#

set -e

LAYER=$1
CHR=$2

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"

# Parse config
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
LAYER_DIR="$OUTPUT_BASE/data/layer${LAYER}"

# Input files
ANCIENT_FILTERED="$OUTPUT_BASE/data/processed/ancient/chr${CHR}.ancient.filtered.vcf.gz"
UKB_PROCESSED="$OUTPUT_BASE/data/processed/ukb/chr${CHR}.ukb.vcf.gz"
COMMON_SITES="$OUTPUT_BASE/data/processed/common_sites/chr${CHR}_common_sites.txt"
REF_SAMPLES="$BASE_DIR/config/samples/layer1/ref_samples.txt"
UKB_SUBSET="$BASE_DIR/config/samples/ukb_subset.txt"

echo "Step 1: Preparing ref and target VCFs for chr$CHR"

# Check if already done
REF_VCF="$LAYER_DIR/ref/chr${CHR}.ref.vcf.gz"
TARGET_VCF="$LAYER_DIR/target/chr${CHR}.target.vcf.gz"

if [ -f "$REF_VCF" ] && [ -f "$TARGET_VCF" ]; then
    echo "  âœ“ Already exists: $REF_VCF and $TARGET_VCF"
    exit 0
fi

# Extract reference samples at common sites from ancient data
echo "  Extracting reference samples from ancient data..."
bcftools view \
    -S $REF_SAMPLES \
    -R $COMMON_SITES \
    $ANCIENT_FILTERED \
    -Oz -o $REF_VCF

tabix -p vcf $REF_VCF

N_REF=$(bcftools query -l $REF_VCF | wc -l)
N_SITES=$(bcftools view -H $REF_VCF | wc -l)
echo "  âœ“ Reference VCF: $N_REF samples, $N_SITES sites"

# Extract target samples (UKB subset) at common sites
echo "  Extracting target samples from UKB..."

if [ -f "$UKB_SUBSET" ]; then
    bcftools view \
        -S $UKB_SUBSET \
        -R $COMMON_SITES \
        $UKB_PROCESSED \
        -Oz -o $TARGET_VCF
else
    echo "  Warning: UKB subset file not found: $UKB_SUBSET"
    echo "  Using all UKB samples (this may be slow)"
    bcftools view \
        -R $COMMON_SITES \
        $UKB_PROCESSED \
        -Oz -o $TARGET_VCF
fi

tabix -p vcf $TARGET_VCF

N_TARGET=$(bcftools query -l $TARGET_VCF | wc -l)
echo "  âœ“ Target VCF: $N_TARGET samples, $N_SITES sites"

echo "  âœ“ Step 1 complete for chr$CHR"
