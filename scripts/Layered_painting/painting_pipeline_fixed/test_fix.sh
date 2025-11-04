#!/bin/bash
#
# Quick test to verify the bcftools fix is applied
#

echo "========================================"
echo "Testing Fixed Pipeline"
echo "========================================"
echo ""

# Check if the fixed syntax is present
echo "1. Checking for fixed bcftools syntax..."
if grep -q "F_MISSING<" scripts/phase1_data_prep/run_all.sh; then
    echo "   ✓ Fixed syntax found: F_MISSING field is used"
else
    echo "   ✗ Error: Fixed syntax not found!"
    exit 1
fi

# Check that old broken syntax is NOT present
echo "2. Checking old broken syntax is removed..."
if grep -q "max-missing-count" scripts/phase1_data_prep/run_all.sh; then
    echo "   ✗ Error: Old broken syntax still present!"
    exit 1
else
    echo "   ✓ Old broken syntax removed"
fi

# Check that combined filter is used (efficiency improvement)
echo "3. Checking for efficient combined filter..."
if grep -q "INFO/INFO>.*&&.*MAF>.*&&.*F_MISSING<" scripts/phase1_data_prep/run_all.sh; then
    echo "   ✓ Using efficient single-pass filter"
else
    echo "   ✗ Warning: Combined filter may not be optimal"
fi

# Check all required files are present
echo "4. Checking required files..."
REQUIRED_FILES=(
    "run_pipeline.sh"
    "scripts/phase1_data_prep/run_all.sh"
    "scripts/phase2_layer_painting/run_layer.sh"
    "scripts/phase0_setup/validate_config.py"
    "config/master_config.yaml"
    "config/layers/layer1_config.yaml"
    "README.md"
)

ALL_PRESENT=true
for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "   ✓ $file"
    else
        echo "   ✗ Missing: $file"
        ALL_PRESENT=false
    fi
done

if [ "$ALL_PRESENT" = false ]; then
    exit 1
fi

echo ""
echo "========================================"
echo "✓ All Tests Passed!"
echo "========================================"
echo ""
echo "The pipeline is ready to use."
echo ""
echo "Next steps:"
echo "  1. Edit config/master_config.yaml"
echo "  2. Create config/samples/ukb_subset.txt"
echo "  3. Run: ./run_pipeline.sh validate"
echo "  4. Run: ./run_pipeline.sh prep"
