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

# Run analyses for multiple pathogens using REFERENCE STRAINS
# Syntax: run_analysis "NCBI query" "output_dir" "name" [n_jobs]
#
# Using reference strains provides:
# - Higher quality annotations
# - Reduced redundancy
# - Literature-comparable results
# - Faster analysis (5-50x fewer proteins)

echo "=========================================="
echo "VIRAL PATHOGENS"
echo "=========================================="

# Small viral proteomes (can use more parallel jobs)
run_analysis "Variola virus strain Bangladesh-1975[Organism]" \
             "./smallpox_analysis" \
             "smallpox" \
             40

run_analysis "Human immunodeficiency virus 1 HXB2[Organism]" \
             "./HIV_analysis" \
             "HIV" \
             40

run_analysis "Measles morbillivirus strain Edmonston[Organism]" \
             "./measles_analysis" \
             "measles" \
             40

run_analysis "Hepatitis B virus subtype adw2[Organism]" \
             "./hepatitis_b_analysis" \
             "hepatitis_b" \
             40

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS - MAJOR HISTORICAL"
echo "=========================================="

# Plague - Black Death, Justinian Plague
run_analysis "Yersinia pestis CO92[Organism]" \
             "./plague_analysis" \
             "plague" \
             30

# Tuberculosis - major historical killer
run_analysis "Mycobacterium tuberculosis H37Rv[Organism]" \
             "./tuberculosis_analysis" \
             "tuberculosis" \
             30

# Leprosy - medieval stigmatized disease
run_analysis "Mycobacterium leprae TN[Organism]" \
             "./leprosy_analysis" \
             "leprosy" \
             30

# Typhoid - major waterborne killer
run_analysis "Salmonella enterica subsp. enterica serovar Typhi str. CT18[Organism]" \
             "./typhoid_analysis" \
             "typhoid" \
             30

# Syphilis - post-Columbian epidemic in Europe
run_analysis "Treponema pallidum subsp. pallidum str. Nichols[Organism]" \
             "./syphilis_analysis" \
             "syphilis" \
             30

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS - RESPIRATORY"
echo "=========================================="

# Pneumonia - leading historical cause of death
run_analysis "Streptococcus pneumoniae TIGR4[Organism]" \
             "./pneumonia_analysis" \
             "pneumonia" \
             30

# Diphtheria - major childhood killer
run_analysis "Corynebacterium diphtheriae NCTC 13129[Organism]" \
             "./diphtheria_analysis" \
             "diphtheria" \
             30

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS - STREPTOCOCCAL"
echo "=========================================="

# Scarlet fever, puerperal fever, necrotizing fasciitis
run_analysis "Streptococcus pyogenes M1 GAS[Organism]" \
             "./scarlet_fever_analysis" \
             "scarlet_fever" \
             30

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS - VECTOR-BORNE"
echo "=========================================="

# Epidemic typhus - war and famine disease
run_analysis "Rickettsia prowazekii str. Madrid E[Organism]" \
             "./typhus_analysis" \
             "typhus" \
             30

# Relapsing fever - louse-borne epidemic disease
run_analysis "Borrelia recurrentis A1[Organism]" \
             "./relapsing_fever_analysis" \
             "relapsing_fever" \
             30

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS - ZOONOTIC"
echo "=========================================="

# Anthrax - occupational disease of farmers/herders
run_analysis "Bacillus anthracis str. Ames[Organism]" \
             "./anthrax_analysis" \
             "anthrax" \
             30

# Brucellosis - from livestock, undulant fever
run_analysis "Brucella melitensis 16M[Organism]" \
             "./brucellosis_analysis" \
             "brucellosis" \
             30

# Listeriosis - foodborne, important in pregnant women
run_analysis "Listeria monocytogenes EGD-e[Organism]" \
             "./listeriosis_analysis" \
             "listeriosis" \
             30

echo ""
echo "=========================================="
echo "BACTERIAL PATHOGENS - WATERBORNE"
echo "=========================================="

# Cholera - 19th century pandemics in Europe
run_analysis "Vibrio cholerae O1 biovar El Tor str. N16961[Organism]" \
             "./cholera_analysis" \
             "cholera" \
             30

# Leptospirosis - waterborne, occupational
run_analysis "Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130[Organism]" \
             "./leptospirosis_analysis" \
             "leptospirosis" \
             20

echo ""
echo "=========================================="
echo "PARASITIC PATHOGENS"
echo "=========================================="

# Malaria - major selective pressure in Mediterranean/southern Europe
run_analysis "Plasmodium vivax Sal-1[Organism]" \
             "./malaria_vivax_analysis" \
             "malaria_vivax" \
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
