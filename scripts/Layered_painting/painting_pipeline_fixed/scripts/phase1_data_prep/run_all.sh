#!/bin/bash
#
# Phase 1: Data Preparation
# Prepares UKB and ancient VCF data for painting
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASE_DIR="$(dirname $(dirname $SCRIPT_DIR))"

# Source configuration
CONFIG_FILE="$BASE_DIR/config/master_config.yaml"

# Parse YAML config (simple extraction)
OUTPUT_BASE=$(grep "output_base:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
ANCIENT_VCF=$(grep "ancient_vcf:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
UKB_HAPLOTYPE=$(grep "ukb_haplotype:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
UKB_CHR_MAP=$(grep "ukb_chr_rename_map:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
UKB_SAMPLES=$(grep "ukb_samples_file:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')
CHROMOSOMES=$(grep "chromosomes:" $CONFIG_FILE | grep -oP '\[\K[0-9, ]+' | tr -d '[]')
PBWT=$(grep "pbwt:" $CONFIG_FILE | head -1 | awk '{print $2}' | tr -d '"')

INFO_THRESH=$(grep "info_threshold:" $CONFIG_FILE | awk '{print $2}')
MAF_THRESH=$(grep "maf_global:" $CONFIG_FILE | awk '{print $2}')
MAX_MISSING=$(grep "max_missing_global:" $CONFIG_FILE | awk '{print $2}')

echo "=================================================="
echo "Phase 1: Data Preparation"
echo "=================================================="
echo "Output directory: $OUTPUT_BASE"
echo "Chromosomes: $CHROMOSOMES"
echo ""

# Load required modules
echo "Loading required modules..."
module load python/3.9.5
module load bcftools/1.20
module load gsl/2.5

# Create directory structure
mkdir -p $OUTPUT_BASE/data/{raw,processed/ukb,processed/ancient,processed/common_sites}
mkdir -p $OUTPUT_BASE/logs

# Extract sample lists
echo "Step 0: Extracting sample lists from reference popfile..."
python3 $SCRIPT_DIR/../utils/extract_sample_lists.py \
    $BASE_DIR/config/samples/layer1/ref_popfile.txt \
    $BASE_DIR/config/samples/layer1/

# Process each chromosome
for CHR in $(echo $CHROMOSOMES | tr ',' ' '); do
    echo ""
    echo "=================================================="
    echo "Processing Chromosome $CHR"
    echo "=================================================="
    
    # Expand chromosome in paths
    ANCIENT_VCF_CHR="${ANCIENT_VCF//\{chrom\}/$CHR}"
    UKB_VCF_CHR="${UKB_HAPLOTYPE//\{chrom\}/$CHR}"
    UKB_MAP_CHR="${UKB_CHR_MAP//\{chrom\}/$CHR}"
    
    # Step 1: Preprocess UKB VCF (extract hard calls, add chr info and sample names)
    echo "Step 1: Preprocessing UKB VCF for chr$CHR..."

    UKB_PROCESSED="$OUTPUT_BASE/data/processed/ukb/chr${CHR}.ukb.vcf.gz"

    if [ ! -f "$UKB_PROCESSED" ]; then
        echo "  Extracting hard calls (GT) with pbwt..."
        $PBWT -readVcfGT $UKB_VCF_CHR \
            -writeVcf $OUTPUT_BASE/data/processed/ukb/chr${CHR}.GT.vcf.gz

        echo "  Adding chromosome info..."
        bcftools annotate $OUTPUT_BASE/data/processed/ukb/chr${CHR}.GT.vcf.gz \
            --rename-chrs $UKB_MAP_CHR \
            -Oz -o $OUTPUT_BASE/data/processed/ukb/chr${CHR}.GT.chrinfo.vcf.gz

        echo "  Adding sample names and finalizing..."
        bcftools reheader \
            --samples $UKB_SAMPLES \
            $OUTPUT_BASE/data/processed/ukb/chr${CHR}.GT.chrinfo.vcf.gz \
            -o $UKB_PROCESSED

        tabix -p vcf $UKB_PROCESSED

        # Clean up intermediate files
        rm $OUTPUT_BASE/data/processed/ukb/chr${CHR}.GT.vcf.gz
        rm $OUTPUT_BASE/data/processed/ukb/chr${CHR}.GT.chrinfo.vcf.gz
        echo "  ✓ Created $UKB_PROCESSED"
    else
        echo "  ✓ Already exists: $UKB_PROCESSED"
    fi
    
    # Step 2: Filter ancient VCF for quality
    echo "Step 2: Filtering ancient VCF for quality..."
    
    ANCIENT_FILTERED="$OUTPUT_BASE/data/processed/ancient/chr${CHR}.ancient.filtered.vcf.gz"
    
    if [ ! -f "$ANCIENT_FILTERED" ]; then
        echo "  Applying filters: INFO>$INFO_THRESH, MAF>$MAF_THRESH, F_MISSING<$MAX_MISSING"
        
        # FIXED: Using single bcftools view with combined filter expression
        bcftools view \
            $ANCIENT_VCF_CHR \
            -i "INFO/INFO>$INFO_THRESH && MAF>$MAF_THRESH && F_MISSING<$MAX_MISSING" \
            -Oz -o $ANCIENT_FILTERED
        
        tabix -p vcf $ANCIENT_FILTERED
        
        N_VARIANTS=$(bcftools view -H $ANCIENT_FILTERED | wc -l)
        echo "  ✓ Created $ANCIENT_FILTERED ($N_VARIANTS variants)"
    else
        N_VARIANTS=$(bcftools view -H $ANCIENT_FILTERED | wc -l)
        echo "  ✓ Already exists: $ANCIENT_FILTERED ($N_VARIANTS variants)"
    fi
    
    # Step 3: Find common sites between UKB and ancient data
    echo "Step 3: Finding common sites..."
    
    COMMON_SITES="$OUTPUT_BASE/data/processed/common_sites/chr${CHR}_common_sites.txt"
    
    if [ ! -f "$COMMON_SITES" ]; then
        # Extract positions from both files and find intersection
        bcftools query -f '%CHROM\t%POS\n' $UKB_PROCESSED | sort -k1,1 -k2,2n > $OUTPUT_BASE/data/processed/common_sites/chr${CHR}_ukb_sites.txt
        bcftools query -f '%CHROM\t%POS\n' $ANCIENT_FILTERED | sort -k1,1 -k2,2n > $OUTPUT_BASE/data/processed/common_sites/chr${CHR}_ancient_sites.txt
        
        comm -12 \
            $OUTPUT_BASE/data/processed/common_sites/chr${CHR}_ukb_sites.txt \
            $OUTPUT_BASE/data/processed/common_sites/chr${CHR}_ancient_sites.txt \
            > $COMMON_SITES
        
        N_COMMON=$(wc -l < $COMMON_SITES)
        echo "  ✓ Found $N_COMMON common sites"
        
        # Clean up temp files
        rm $OUTPUT_BASE/data/processed/common_sites/chr${CHR}_ukb_sites.txt
        rm $OUTPUT_BASE/data/processed/common_sites/chr${CHR}_ancient_sites.txt
    else
        N_COMMON=$(wc -l < $COMMON_SITES)
        echo "  ✓ Already exists: $COMMON_SITES ($N_COMMON sites)"
    fi
    
    echo "Chromosome $CHR complete!"
done

echo ""
echo "=================================================="
echo "Phase 1 Complete!"
echo "=================================================="
echo ""
echo "Next steps:"
echo "  1. Review the common sites counts"
echo "  2. Create a UKB subset file: config/samples/ukb_subset.txt"
echo "  3. Run Phase 2: ./run_pipeline.sh layer 1"
