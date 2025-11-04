#!/bin/bash
#
# Master orchestrator for Temporal Ancestry Painting Pipeline
# Usage: ./run_pipeline.sh {prep|layer <N>|all-layers|postprocess|full}
#

set -e  # Exit on error

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if config is valid
check_config() {
    print_info "Validating configuration..."
    if python3 scripts/phase0_setup/validate_config.py config/; then
        print_info "âœ“ Configuration valid"
        return 0
    else
        print_error "âœ— Configuration validation failed"
        return 1
    fi
}

# Function to run Phase 1: Data Preparation
run_prep() {
    print_info "========================================="
    print_info "Phase 1: Data Preparation"
    print_info "========================================="
    
    if [ ! -f "scripts/phase1_data_prep/run_all.sh" ]; then
        print_error "Phase 1 script not found: scripts/phase1_data_prep/run_all.sh"
        return 1
    fi
    
    bash scripts/phase1_data_prep/run_all.sh
    
    if [ $? -eq 0 ]; then
        print_info "âœ“ Phase 1 completed successfully"
        return 0
    else
        print_error "âœ— Phase 1 failed"
        return 1
    fi
}

# Function to run Phase 2: Paint specific layer
run_layer() {
    LAYER=$1
    
    if [ -z "$LAYER" ]; then
        print_error "Layer number not specified"
        echo "Usage: $0 layer <N>"
        return 1
    fi
    
    # Check if layer config exists
    LAYER_CONFIG="config/layers/layer${LAYER}_config.yaml"
    if [ ! -f "$LAYER_CONFIG" ]; then
        print_error "Layer configuration not found: $LAYER_CONFIG"
        return 1
    fi
    
    print_info "========================================="
    print_info "Phase 2: Painting Layer $LAYER"
    print_info "========================================="
    
    if [ ! -f "scripts/phase2_layer_painting/run_layer.sh" ]; then
        print_error "Phase 2 script not found: scripts/phase2_layer_painting/run_layer.sh"
        return 1
    fi
    
    bash scripts/phase2_layer_painting/run_layer.sh "$LAYER"
    
    if [ $? -eq 0 ]; then
        print_info "âœ“ Layer $LAYER completed successfully"
        return 0
    else
        print_error "âœ— Layer $LAYER failed"
        return 1
    fi
}

# Function to run all layers in parallel
run_all_layers() {
    print_info "========================================="
    print_info "Phase 2: Painting All Layers (Parallel)"
    print_info "========================================="
    
    # Find all layer config files
    LAYER_CONFIGS=(config/layers/layer*_config.yaml)
    
    if [ ${#LAYER_CONFIGS[@]} -eq 0 ]; then
        print_error "No layer configurations found in config/layers/"
        return 1
    fi
    
    print_info "Found ${#LAYER_CONFIGS[@]} layer configurations"
    
    # Submit each layer as a separate job (or run sequentially for testing)
    for config in "${LAYER_CONFIGS[@]}"; do
        # Extract layer number from filename
        LAYER=$(echo "$config" | grep -oP 'layer\K\d+')
        
        print_info "Submitting Layer $LAYER..."
        
        # Option 1: Submit as SLURM job (uncomment when ready)
        # sbatch --job-name="paint_layer${LAYER}" \
        #        --output="logs/layer${LAYER}_%j.out" \
        #        --error="logs/layer${LAYER}_%j.err" \
        #        scripts/phase2_layer_painting/run_layer.sh "$LAYER"
        
        # Option 2: Run sequentially (for testing)
        bash scripts/phase2_layer_painting/run_layer.sh "$LAYER" &
    done
    
    # Wait for all background jobs to complete
    print_info "Waiting for all layers to complete..."
    wait
    
    print_info "âœ“ All layers completed"
    return 0
}

# Function to run Phase 3: Post-processing
run_postprocess() {
    print_info "========================================="
    print_info "Phase 3: Post-Processing"
    print_info "========================================="

    if [ ! -f "scripts/phase3_postprocess/combine_all_layers.sh" ]; then
        print_warning "Phase 3 is not yet implemented."
        print_warning "Script not found: scripts/phase3_postprocess/combine_all_layers.sh"
        print_info ""
        print_info "This feature is planned for a future release."
        print_info "It will include:"
        print_info "  - Combining results from all layers"
        print_info "  - Creating SQLite database"
        print_info "  - VCF annotation with ancestry information"
        return 1
    fi

    # Check if any layers have been completed
    if [ ! -d "data/results" ] || [ -z "$(ls -A data/results/)" ]; then
        print_warning "No layer results found in data/results/"
        print_warning "Have you run any painting layers yet?"
        return 1
    fi

    bash scripts/phase3_postprocess/combine_all_layers.sh

    if [ $? -eq 0 ]; then
        print_info "Phase 3 completed successfully"
        return 0
    else
        print_error "âœ— Phase 3 failed"
        return 1
    fi
}

# Function to run complete pipeline
run_full() {
    print_info "========================================="
    print_info "Running Full Pipeline"
    print_info "========================================="

    # Step 1: Validate config
    if ! check_config; then
        return 1
    fi

    # Step 2: Data preparation
    print_info "\nStep 1/2: Data Preparation"
    if ! run_prep; then
        print_error "Pipeline stopped due to Phase 1 failure"
        return 1
    fi

    # Step 3: Paint all layers
    print_info "\nStep 2/2: Painting All Layers"
    if ! run_all_layers; then
        print_error "Pipeline stopped due to Phase 2 failure"
        return 1
    fi

    print_info "========================================="
    print_info "Full pipeline completed successfully!"
    print_info "========================================="
    print_info ""
    print_info "Note: Phase 3 (post-processing) is not yet implemented."
    print_info "You can run it manually when available with: ./run_pipeline.sh postprocess"
    return 0
}

# Function to show usage
show_usage() {
    cat << EOF
Temporal Ancestry Painting Pipeline

Usage: $0 <command> [options]

Commands:
    prep              Run Phase 1: Data Preparation
    layer <N>         Run Phase 2: Paint specific layer (1-5)
    all-layers        Run Phase 2: Paint all layers in parallel
    postprocess       Run Phase 3: Post-processing and database creation
    full              Run complete pipeline (prep + all layers + postprocess)
    validate          Validate configuration files
    help              Show this help message

Examples:
    $0 prep                    # Prepare data
    $0 layer 1                 # Paint layer 1
    $0 all-layers              # Paint all layers in parallel
    $0 postprocess             # Build database from results
    $0 full                    # Run everything

Configuration:
    Edit config/master_config.yaml and config/layers/layerN_config.yaml
    before running the pipeline.

For more information, see README.md
EOF
}

# Main script logic
COMMAND=$1
LAYER=$2

case $COMMAND in
    prep)
        check_config || exit 1
        run_prep
        ;;
        
    layer)
        check_config || exit 1
        run_layer "$LAYER"
        ;;
        
    all-layers)
        check_config || exit 1
        run_all_layers
        ;;
        
    postprocess)
        run_postprocess
        ;;
        
    full)
        run_full
        ;;
        
    validate)
        check_config
        ;;
        
    help|--help|-h)
        show_usage
        exit 0
        ;;
        
    *)
        print_error "Unknown command: $COMMAND"
        echo ""
        show_usage
        exit 1
        ;;
esac

exit $?
