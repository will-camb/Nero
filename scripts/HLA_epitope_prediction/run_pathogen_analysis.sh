#!/bin/bash
# Batch runner for HLA epitope prediction across multiple pathogens
# Uses improved incremental processing pipeline

set -e  # Exit on error

# Configuration
EMAIL="wb275@cam.ac.uk"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/hla_epitope_predictor.py"

# Tool paths (modify if needed)
NETMHCPAN_PATH="/home/jdn321/software/netMHCpan-4.2/netMHCpan"
NETMHCIIPAN_PATH="/home/jdn321/software/netMHCIIpan-4.3/netMHCIIpan"
SIGNALP_MODELS="/home/jdn321/software/signalp6_fast/signalp-6-package/models/"

# Allele files (use defaults if not specified)
CLASS_I_ALLELES="${SCRIPT_DIR}/alleles_class_i_default.txt"
CLASS_II_ALLELES="${SCRIPT_DIR}/alleles_class_ii_default.txt"

# Number of parallel jobs (reduced from 80 to prevent memory issues)
# For large proteomes (>3000 proteins), use 20-30
# For small proteomes (<500 proteins), can use 40-60
N_JOBS=30

# Peptide generation parameters (conservative defaults for ~10x speedup)
# Class I: 9-mers only (instead of 8,9,10,11)
# Step size: 3 (instead of 1 = fully overlapping)
CLASS_I_LENGTHS="9"
CLASS_I_STEP=3
CLASS_II_LENGTH=15
CLASS_II_STEP=3

echo "============================================"
echo "HLA Epitope Prediction Pipeline"
echo "============================================"
echo "Email: $EMAIL"
echo "Parallel jobs: $N_JOBS"
echo "Class I: lengths=$CLASS_I_LENGTHS, step=$CLASS_I_STEP"
echo "Class II: length=$CLASS_II_LENGTH, step=$CLASS_II_STEP"
echo "Class I alleles: $CLASS_I_ALLELES"
echo "Class II alleles: $CLASS_II_ALLELES"
echo ""

# Function to run analysis for one organism
run_analysis() {
    local organism_query="$1"
    local output_dir="$2"
    local organism_name="$3"
    local n_jobs="${4:-$N_JOBS}"  # Use default if not provided

    echo ""
    echo "========================================"
    echo "Analyzing: $organism_name"
    echo "Output: $output_dir"
    echo "Parallel jobs: $n_jobs"
    echo "========================================"

    python3 "$SCRIPT" \
        --email "$EMAIL" \
        --organism "$organism_query" \
        --output-dir "$output_dir" \
        --name "$organism_name" \
        --class-i-alleles "$CLASS_I_ALLELES" \
        --class-ii-alleles "$CLASS_II_ALLELES" \
        --class-i-lengths $CLASS_I_LENGTHS \
        --class-i-step "$CLASS_I_STEP" \
        --class-ii-length "$CLASS_II_LENGTH" \
        --class-ii-step "$CLASS_II_STEP" \
        --netmhcpan-path "$NETMHCPAN_PATH" \
        --netmhciipan-path "$NETMHCIIPAN_PATH" \
        --signalp-model-dir "$SIGNALP_MODELS" \
        --n-jobs "$n_jobs"

    if [ $? -eq 0 ]; then
        echo "✓ $organism_name completed successfully"
    else
        echo "✗ $organism_name failed"
        return 1
    fi
}

# Run analyses for multiple pathogens
# Syntax: run_analysis "NCBI query" "output_dir" "name" [n_jobs]

# Small proteomes (can use more parallel jobs)
run_analysis "Variola virus[Organism] AND RefSeq[Filter]" \
             "./smallpox_analysis" \
             "smallpox" \
             40

run_analysis "HIV[Organism] AND RefSeq[Filter]" \
             "./HIV_analysis" \
             "HIV" \
             40

run_analysis "measles[Organism] AND RefSeq[Filter]" \
             "./measles_analysis" \
             "measles" \
             40

# Medium proteomes
run_analysis "Yersinia pestis CO92[Organism] AND RefSeq[Filter]" \
             "./plague_analysis" \
             "plague" \
             30

run_analysis "Mycobacterium leprae[Organism] AND RefSeq[Filter]" \
             "./leprosy_analysis" \
             "leprosy" \
             30

run_analysis "Borrelia recurrentis[Organism] AND RefSeq[Filter]" \
             "./B_recurrentis_analysis" \
             "B_recurrentis" \
             30

run_analysis "Treponema pallidum[Organism] AND RefSeq[Filter]" \
             "./T_pallidum_analysis" \
             "T_pallidum" \
             30

# Large proteomes (reduce parallel jobs to avoid memory issues)
run_analysis "Leptospira interrogans Copenhageni[Organism] AND RefSeq[Filter]" \
             "./L_interrogans_analysis" \
             "L_interrogans" \
             20

run_analysis "Mycobacterium tuberculosis H37Rv[Organism] AND RefSeq[Filter]" \
             "./M_tuberculosis_analysis" \
             "M_tuberculosis" \
             20

run_analysis "Plasmodium vivax[Organism] AND RefSeq[Filter]" \
             "./P_vivax_analysis" \
             "P_vivax" \
             20

echo ""
echo "============================================"
echo "All analyses complete!"
echo "============================================"
echo ""
echo "Next steps:"
echo "1. Review results in each *_analysis/ directory"
echo "2. Check *_hla_predictions.csv for binding predictions"
echo "3. Run downstream analysis (e.g., aggregate by allele)"
echo ""
